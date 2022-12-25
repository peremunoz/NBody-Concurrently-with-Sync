/* ---------------------------------------------------------------
Práctica 1.
Código fuente : nbfast.c
Grau Informàtica
Pere Muñoz Figuerol - 48252062V
--------------------------------------------------------------- */

// Usage: NBody* <particle number> <iterations number> [particle file] [0: console mode / other: graphic mode] [threads number]
// You can leave any of the [optional] parameters empty with the - character

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<stdio.h>
#include<unistd.h>
#include<pthread.h>
#include<semaphore.h>
#include <stdbool.h>

#ifdef D_GLFW_SUPPORT
    #include<GLFW/glfw3.h>
#endif

// Macros to make code a little bit easier to understand because for speedup reasons, I'll use only 1D arrays
#define PX(i) (3*i+1)
#define PY(i) (3*i+2)
#define MASS(i) (3*i+3)

#define VX(i) (4*i+0)
#define VY(i) (4*i+1)
#define AX(i) (4*i+2)
#define AY(i) (4*i+3)

double G=0.0001;
double dt=0.005;
double rcutoff=0.35;
double rlimit=0.03;

struct Node{
    struct Node *children[4];
    int external;

    double CMX;
    double CMY;
    double mass;
    double TRX;
    double TRY;

    double LLX;
    double LLY;

    double GCX;
    double GCY;
};

struct GraphicInterfaceStruct{
    struct Node* tree;
    double* sharedBuff;
    int* indexes;
    int nShared;
    int nLocal;
    double* localBuff;
    double* radius;
    int count;
    int steps;
    int nOriginal;
    int M;
};

struct BuildTreeStruct{
    struct Node* tree;
    double* sharedBuff;
    int* indexes;
    int nShared;
    int remainingThreads;
};

struct CalculateForceStruct{
    struct Node* tree;
    double* sharedBuff;
    double* localBuff;
    int index;
    int remainingThreads;
};

struct GlobalStruct{
    int nLocal;
    double* localBuff;
    int *indexes;
    double *sharedBuff;
    struct Node *tree;
    int nShared;
    bool lastIteration;
};

// Global synchronization variables
sem_t calculateForceSem;
pthread_mutex_t remainingParticlesMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t calculateForceBarrier;
sem_t endedCalculatingForce;

// Global variables
pthread_t *threads;
int numberOfThreads;
int remainingParticles;
struct GlobalStruct globalStruct;

void cancelThreads(pthread_t* tids, int nThreads){
    for(int i=0;i<nThreads;i++){
        pthread_cancel(tids[i]);
    }
    exit(-1);
}

void destroySyncVariables(){

}

void buildTree(struct Node* node, double* shrdBuff, int *indexes, int n){
    if(n==1){ //This is an external node!
        node->external=1;
        node->CMX=shrdBuff[PX(indexes[0])];
        node->CMY=shrdBuff[PY(indexes[0])];
        node->mass=shrdBuff[MASS(indexes[0])];
    } else {
        node->external=0;
        //Arrays of indexes of particles per quartile
        int *NEi = (int *) malloc(sizeof(int)*n);
        int *NWi = (int *) malloc(sizeof(int)*n);
        int *SWi = (int *) malloc(sizeof(int)*n);
        int *SEi = (int *) malloc(sizeof(int)*n);
        int NWc=0, SWc=0,SEc=0, NEc=0;

        int i;
        /** For each particle we will check where is it located relative to the geometric center,
            to sort them into the 4 children nodes**/
        for(i=0;i<n;i++){
            if(shrdBuff[PY(indexes[i])] < node->GCY ){ //South half
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    SWi[SWc]=indexes[i];
                    SWc++;
                } else {
                    SEi[SEc]=indexes[i];
                    SEc++;
                }
            } else { //North half
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    NWi[NWc]=indexes[i];
                    NWc++;
                } else {
                    NEi[NEc]=indexes[i];
                    NEc++;
                }
            }
        }
        //If there are particles in the NorthWest quarter
        if(NEc>0){
            //This instruction declares a new node in the position 0
            node->children[0]= malloc(sizeof *node->children[0]);
            //We give the values of the Low Left and Top Right corner, and also the geometric center.
            node->children[0]->TRX=node->TRX;
            node->children[0]->TRY=node->TRY;
            node->children[0]->LLX=node->GCX;
            node->children[0]->LLY=node->GCY;
            node->children[0]->GCX=(node->GCX+node->TRX)/2;
            node->children[0]->GCY=(node->GCY+node->TRY)/2;
            //We build a tree in the new node, with the particles that are inside
            buildTree(node->children[0],shrdBuff,NEi,NEc);
        } else {
            //If not, we set the children to null
            node->children[0]=NULL;
        }
        //The next three blocks are exactly the same thing but for the other three nodes
        if(NWc>0){
            node->children[1]= malloc(sizeof *node->children[1]);
            node->children[1]->TRX=node->GCX;
            node->children[1]->TRY=node->TRY;
            node->children[1]->LLX=node->LLX;
            node->children[1]->LLY=node->GCY;
            node->children[1]->GCX=(node->LLX+node->GCX)/2;
            node->children[1]->GCY=(node->GCY+node->TRY)/2;
            buildTree(node->children[1],shrdBuff,NWi,NWc);
        } else {
            node->children[1]=NULL;
        }
        if(SWc>0){
            node->children[2]= malloc(sizeof *node->children[2]);
            node->children[2]->TRX=node->GCX;
            node->children[2]->TRY=node->GCY;
            node->children[2]->LLX=node->LLX;
            node->children[2]->LLY=node->LLY;
            node->children[2]->GCX=(node->LLX+node->GCX)/2;
            node->children[2]->GCY=(node->LLY+node->GCY)/2;
            buildTree(node->children[2],shrdBuff,SWi,SWc);
        } else {
            node->children[2]=NULL;
        }
        if(SEc>0){
            node->children[3]= malloc(sizeof *node->children[3]);
            node->children[3]->TRX=node->TRX;
            node->children[3]->TRY=node->GCY;
            node->children[3]->LLX=node->GCX;
            node->children[3]->LLY=node->LLY;
            node->children[3]->GCX=(node->GCX+node->TRX)/2;
            node->children[3]->GCY=(node->LLY+node->GCY)/2;
            buildTree(node->children[3],shrdBuff,SEi,SEc);
        } else {
            node->children[3]=NULL;
        }
        node->mass=0;
        node->CMX=0;
        node->CMY=0;
        //Now that we have finished building the 4 trees beneath this node, we calculate the Center of Mass
        //based on the center of mass of the children
        for(i=0;i<4;i++){
            if(node->children[i]!=NULL){
                node->mass+=node->children[i]->mass;
                node->CMX+=node->children[i]->CMX*node->children[i]->mass;
                node->CMY+=node->children[i]->CMY*node->children[i]->mass;
            }
        }
        node->CMX=node->CMX/node->mass;
        node->CMY=node->CMY/node->mass;
        //And tadaaa
    }
}

void calculateForce(struct Node *tree, double *shrdBuff, double *localBuff, int index){
    double distance = sqrt((tree->CMX-shrdBuff[PX(index)])*(tree->CMX-shrdBuff[PX(index)])+
                           (tree->CMY-shrdBuff[PY(index)])*(tree->CMY-shrdBuff[PY(index)]));
    //First we check if the node is not actually the same particle we are calculating
    if(distance>0){
        //Now, we know it is not because the is some distance between the Center of Mass and the particle
        //If the node is external (only contains one particle) or is far away enough, we calculate the force with the center of mass
        if(distance>rcutoff || tree->external){
            double f;
            if(distance<rlimit){
                f=G*tree->mass/(rlimit*rlimit*distance);
            } else {
                f=G*tree->mass/(distance*distance*distance);
            }
            localBuff[AX(index)]+=f*(tree->CMX-shrdBuff[PX(index)]);
            localBuff[AY(index)]+=f*(tree->CMY-shrdBuff[PY(index)]);
        } else {
            //If not, we recursively call the calculateForce() function in the children that are not empty.
            int i;
            for(i=0;i<4;i++){
                if(tree->children[i]!=NULL){
                    calculateForce(tree->children[i],shrdBuff,localBuff,index);
                }
            }
        }
    }
}

void buildTreeThread(struct BuildTreeStruct* data){
    // Unpack the variables from the data struct
    struct Node* node = data->tree;
    double* shrdBuff = data->sharedBuff;
    int* indexes = data->indexes;
    int n = data->nShared;
    int remainingThreads = data->remainingThreads;

    // Create variables for controlling the concurrency
    int possibleSubThreads = 0;
    int remainingThreadsPerSubThread = 0;
    int extraThreadsPerSubThread = 0;

    pthread_t NE_Tid = 0, NW_Tid = 0, SW_Tid = 0, SE_Tid = 0;
    struct BuildTreeStruct *NE_data=0, *NW_data=0, *SW_data=0, *SE_data=0;

    if(n==1){ //This is an external node!
        node->external=1;
        node->CMX=shrdBuff[PX(indexes[0])];
        node->CMY=shrdBuff[PY(indexes[0])];
        node->mass=shrdBuff[MASS(indexes[0])];
    } else {
        node->external=0;
		//Arrays of indexes of particles per quartile
        int *NEi = (int *) malloc(sizeof(int)*n);
        int *NWi = (int *) malloc(sizeof(int)*n);
        int *SWi = (int *) malloc(sizeof(int)*n);
        int *SEi = (int *) malloc(sizeof(int)*n);
        int NWc=0, SWc=0,SEc=0, NEc=0;

        int i;
		/** For each particle we will check where is it located relative to the geometric center,
			to sort them into the 4 children nodes**/
        for(i=0;i<n;i++){
            if(shrdBuff[PY(indexes[i])] < node->GCY ){ //South half
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    SWi[SWc]=indexes[i];
                    SWc++;
                } else {
                    SEi[SEc]=indexes[i];
                    SEc++;
                }
            } else { //North half
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    NWi[NWc]=indexes[i];
                    NWc++;
                } else {
                    NEi[NEc]=indexes[i];
                    NEc++;
                }
            }
        }

        // Pre-calculate how many possible subThreads we could create
        if(NEc>0) possibleSubThreads++;
        if(NWc>0) possibleSubThreads++;
        if(SWc>0) possibleSubThreads++;
        if(SEc>0) possibleSubThreads++;

        // How many threads can we assign to each subThread
        remainingThreadsPerSubThread = remainingThreads / possibleSubThreads;
        extraThreadsPerSubThread = remainingThreads % possibleSubThreads;

		//If there are particles in the NorthWest quarter
        if(NEc>0){
			//This instruction declares a new node in the position 0
            node->children[0]= malloc(sizeof *node->children[0]);
			//We give the values of the Low Left and Top Right corner, and also the geometric center.
            node->children[0]->TRX=node->TRX;
            node->children[0]->TRY=node->TRY;
            node->children[0]->LLX=node->GCX;
            node->children[0]->LLY=node->GCY;
            node->children[0]->GCX=(node->GCX+node->TRX)/2;
            node->children[0]->GCY=(node->GCY+node->TRY)/2;
			//We build a tree in the new node, with the particles that are inside
            // If we have remaining threads to create, we build the sub-tree concurrently
            if (remainingThreads > 0) {
                // Create the build tree struct only if we have to create a thread
                NE_data = malloc(sizeof(struct BuildTreeStruct));
                NE_data->tree = node->children[0];
                NE_data->sharedBuff = shrdBuff;
                NE_data->indexes = NEi;
                NE_data->nShared = NEc;

                int assignedThreadsToSubThread = remainingThreadsPerSubThread;

                if (extraThreadsPerSubThread > 0) {
                    assignedThreadsToSubThread++;
                    extraThreadsPerSubThread--;
                }

                remainingThreads = remainingThreads - assignedThreadsToSubThread;

                NE_data->remainingThreads = assignedThreadsToSubThread;
                if(pthread_create(&NE_Tid, NULL, (void *(*)(void *)) buildTreeThread, NE_data)){
                    perror("Error creating the buildTreeThread [NE children] thread: ");
                    exit(-1);
                }
            } else { // If we don't have threads remaining, we execute the normal function for saving space and time
                buildTree(node->children[0], shrdBuff, NEi, NEc);
            }
        } else {
			//If not, we set the children to null
            node->children[0]=NULL;
        }
		//The next three blocks are exactly the same thing but for the other three nodes
        if(NWc>0){
            node->children[1]= malloc(sizeof *node->children[1]);
            node->children[1]->TRX=node->GCX;
            node->children[1]->TRY=node->TRY;
            node->children[1]->LLX=node->LLX;
            node->children[1]->LLY=node->GCY;
            node->children[1]->GCX=(node->LLX+node->GCX)/2;
            node->children[1]->GCY=(node->GCY+node->TRY)/2;
            // If we have remaining threads to create, we build the sub-tree concurrently
            if (remainingThreads > 0) {
                // Create the build tree struct
                NW_data = malloc(sizeof(struct BuildTreeStruct));
                NW_data->tree = node->children[1];
                NW_data->sharedBuff = shrdBuff;
                NW_data->indexes = NWi;
                NW_data->nShared = NWc;
                NW_data->remainingThreads = 0;

                int assignedThreadsToSubThread = remainingThreadsPerSubThread;

                if (extraThreadsPerSubThread > 0) {
                    assignedThreadsToSubThread++;
                    extraThreadsPerSubThread--;
                }

                remainingThreads = remainingThreads - assignedThreadsToSubThread;

                NW_data->remainingThreads = assignedThreadsToSubThread;
                if(pthread_create(&NW_Tid, NULL, (void *(*)(void *)) buildTreeThread, NW_data)){
                    perror("Error creating the buildTreeThread [NW children] thread: ");
                    exit(-1);
                }
            } else {
                buildTree(node->children[1], shrdBuff, NWi, NWc);
            }
        } else {
            node->children[1]=NULL;
        }
        if(SWc>0){
            node->children[2]= malloc(sizeof *node->children[2]);
            node->children[2]->TRX=node->GCX;
            node->children[2]->TRY=node->GCY;
            node->children[2]->LLX=node->LLX;
            node->children[2]->LLY=node->LLY;
            node->children[2]->GCX=(node->LLX+node->GCX)/2;
            node->children[2]->GCY=(node->LLY+node->GCY)/2;
            // If we have remaining threads to create, we build the sub-tree concurrently
            if (remainingThreads > 0) {
                // Create the build tree struct
                SW_data = malloc(sizeof(struct BuildTreeStruct));
                SW_data->tree = node->children[2];
                SW_data->sharedBuff = shrdBuff;
                SW_data->indexes = SWi;
                SW_data->nShared = SWc;
                SW_data->remainingThreads = 0;

                int assignedThreadsToSubThread = remainingThreadsPerSubThread;

                if (extraThreadsPerSubThread > 0) {
                    assignedThreadsToSubThread++;
                    extraThreadsPerSubThread--;
                }

                remainingThreads = remainingThreads - assignedThreadsToSubThread;

                SW_data->remainingThreads = assignedThreadsToSubThread;
                if(pthread_create(&SW_Tid, NULL, (void *(*)(void *)) buildTreeThread, SW_data)){
                    perror("Error creating the buildTreeThread [SW children] thread: ");
                    exit(-1);
                }
            } else {
                buildTree(node->children[2], shrdBuff, SWi, SWc);
            }
        } else {
            node->children[2]=NULL;
        }
        if(SEc>0){
            node->children[3]= malloc(sizeof *node->children[3]);
            node->children[3]->TRX=node->TRX;
            node->children[3]->TRY=node->GCY;
            node->children[3]->LLX=node->GCX;
            node->children[3]->LLY=node->LLY;
            node->children[3]->GCX=(node->GCX+node->TRX)/2;
            node->children[3]->GCY=(node->LLY+node->GCY)/2;
            // If we have remaining threads to create, we build the sub-tree concurrently
            if (remainingThreads > 0) {
                // Create the build tree struct
                SE_data = malloc(sizeof(struct BuildTreeStruct));
                SE_data->tree = node->children[3];
                SE_data->sharedBuff = shrdBuff;
                SE_data->indexes = SEi;
                SE_data->nShared = SEc;
                SE_data->remainingThreads = 0;

                int assignedThreadsToSubThread = remainingThreadsPerSubThread;

                if (extraThreadsPerSubThread > 0) {
                    assignedThreadsToSubThread++;
                    extraThreadsPerSubThread--;
                }

                remainingThreads = remainingThreads - assignedThreadsToSubThread;

                SE_data->remainingThreads = assignedThreadsToSubThread;
                if(pthread_create(&SE_Tid, NULL, (void *(*)(void *)) buildTreeThread, SE_data)){
                    perror("Error creating the buildTreeThread [SE children] thread: ");
                    exit(-1);
                }
            } else {
                buildTree(node->children[3], shrdBuff, SEi, SEc);
            }
        } else {
            node->children[3]=NULL;
        }
        node->mass=0;
        node->CMX=0;
        node->CMY=0;
        // Wait for the 4 node trees creation to finish
        if (NE_Tid != 0) {
            if (pthread_join(NE_Tid, NULL)) {
                perror("Error joining the NE children thread: ");
                exit(-1);
            }
            free(NE_data);
        }
        if (NW_Tid != 0) {
            if (pthread_join(NW_Tid, NULL)) {
                perror("Error joining the NW children thread: ");
                exit(-1);
            }
            free(NW_data);
        }
        if (SW_Tid != 0) {
            if (pthread_join(SW_Tid, NULL)) {
                perror("Error joining the SW children thread: ");
                exit(-1);
            }
            free(SW_data);
        }
        if (SE_Tid != 0) {
            if (pthread_join(SE_Tid, NULL)) {
                perror("Error joining the SE children thread: ");
                exit(-1);
            }
            free(SE_data);
        }

		//Now that we have finished building the 4 trees beneath this node, we calculate the Center of Mass
		//based on the center of mass of the children
        for(i=0;i<4;i++){
            if(node->children[i]!=NULL){
                node->mass+=node->children[i]->mass;
                node->CMX+=node->children[i]->CMX*node->children[i]->mass;
                node->CMY+=node->children[i]->CMY*node->children[i]->mass;
            }
        }
        node->CMX=node->CMX/node->mass;
        node->CMY=node->CMY/node->mass;
		//And tadaaa
    }
}

void calculateForceThread(struct CalculateForceStruct* data){
    // Unpack the information from the data struct
    struct Node *tree = data->tree;
    double* shrdBuff = data->sharedBuff;
    double* localBuff = data->localBuff;
    int index = data->index;
    int remainingThreads = data->remainingThreads;

    double distance = sqrt((tree->CMX-shrdBuff[PX(index)])*(tree->CMX-shrdBuff[PX(index)])+
                           (tree->CMY-shrdBuff[PY(index)])*(tree->CMY-shrdBuff[PY(index)]));
    // Necessary variables for concurrency
    int tids_index = 0;
    int possibleSubThreads = 0;
    int remainingThreadsPerSubThread = 0;
    int extraThreadsPerSubThread = 0;

	//First we check if the node is not actually the same particle we are calculating
    if(distance>0){
		//Now, we know it is not because the is some distance between the Center of Mass and the particle
		//If the node is external (only contains one particle) or is far away enough, we calculate the force with the center of mass
        if(distance>rcutoff || tree->external){
            double f;
            if(distance<rlimit){
                f=G*tree->mass/(rlimit*rlimit*distance);
            } else {
                f=G*tree->mass/(distance*distance*distance);
            }
            localBuff[AX(index)]+=f*(tree->CMX-shrdBuff[PX(index)]);
            localBuff[AY(index)]+=f*(tree->CMY-shrdBuff[PY(index)]);
        } else {
			//If not, we recursively call the calculateForceThread() function in the children that are not empty.
            int i;
            //First we calculate the valid children of the node for assigning threads correctly.
            for(i=0;i<4;i++){
                if(tree->children[i]!=NULL){
                    possibleSubThreads++;
                }
            }
            //We calculate how many threads per children have to be assigned.
            remainingThreadsPerSubThread = remainingThreads / possibleSubThreads;
            extraThreadsPerSubThread = remainingThreads % possibleSubThreads;

            pthread_t *tids = malloc(sizeof(pthread_t) * possibleSubThreads);
            struct CalculateForceStruct *forceData = malloc(sizeof(struct CalculateForceStruct) * possibleSubThreads);

            for(i=0;i<4;i++){
                if(tree->children[i]!=NULL){
                    // If there are free threads to be created, we execute the next recursive call concurrently.
                    if(remainingThreads > 0) {
                        // Create the CalculateForce struct only if we have threads available.
                        forceData[tids_index].tree = tree->children[i];
                        forceData[tids_index].sharedBuff = shrdBuff;
                        forceData[tids_index].localBuff = localBuff;
                        forceData[tids_index].index = index;
                        forceData[tids_index].remainingThreads = 0;

                        int assignedThreadsToSubThread = remainingThreadsPerSubThread;

                        if (extraThreadsPerSubThread > 0) {
                            assignedThreadsToSubThread++;
                            extraThreadsPerSubThread--;
                        }

                        forceData[tids_index].remainingThreads = assignedThreadsToSubThread;

                        remainingThreads = remainingThreads - assignedThreadsToSubThread;
                        if(pthread_create(&tids[tids_index], NULL, (void *(*)(void *)) calculateForceThread, &forceData[tids_index])){
                            perror("Error creating the calculateForceThread thread: ");
                            cancelThreads(tids, tids_index);
                        }
                        tids_index++;
                    } else { // If we don't have threads, we execute the sequential function without doing any malloc
                        calculateForce(tree->children[i], shrdBuff, localBuff, index);
                    }
                }
            }
            for(i=0;i<tids_index;i++) {
                if(pthread_join(tids[i], NULL)){
                    perror("Error joining the calculateForceThread thread: ");
                    cancelThreads(tids, tids_index);
                }
            }
            free(forceData);
            free(tids);
            tids_index = 0;
        }
    }
}

void moveParticle(double *shrdBuff, double *localBuff, int index){
    //Unprecise but fast euler method for solving the time differential equation
    double oldX=shrdBuff[PX(index)];
    double oldY=shrdBuff[PY(index)];
    shrdBuff[PX(index)]+=localBuff[VX(index)]*dt+localBuff[AX(index)]*dt*dt*0.5;
    shrdBuff[PY(index)]+=localBuff[VY(index)]*dt+localBuff[AY(index)]*dt*dt*0.5;
    localBuff[VX(index)]=(shrdBuff[PX(index)]-oldX)/dt;
    localBuff[VY(index)]=(shrdBuff[PY(index)]-oldY)/dt;
}

void threadFunction(int id) {
    while (1) {
        printf("Thread %d is waiting for a job\n", id);
        // Wait for the main thread to signal that the tree is built and the particles are initialized.
        sem_wait(&calculateForceSem);
        printf("Thread %d is calculating the force\n", id);


        // Calculate index of the particles to be calculated.
        int particlesPerThread = globalStruct.nLocal / numberOfThreads;

        int start = id * particlesPerThread + 1;
        int end = start + particlesPerThread - 1;

        pthread_mutex_lock(&remainingParticlesMutex);
        if (remainingParticles > 0 && id == numberOfThreads - 1) {
            end += remainingParticles;
            remainingParticles = 0;
        }
        pthread_mutex_unlock(&remainingParticlesMutex);

        if (end > globalStruct.nLocal) {
            end = globalStruct.nLocal;
        }

        double *localBuff = globalStruct.localBuff;
        int *indexes = globalStruct.indexes;
        double *sharedBuff = globalStruct.sharedBuff;
        struct Node* tree = globalStruct.tree;

        printf("Thread %d is calculating the force for particles %d to %d\n", id, start, end);
        // Then, we can calculate the forces.
        int i;
        for(i=start; i < end; i++){
            //Set initial accelerations to zero
            localBuff[AX(indexes[i])]=0;
            localBuff[AY(indexes[i])]=0;
            int s;

            for(s=0;s<4;s++){
                //Recursively calculate accelerations
                if(tree->children[s]!=NULL){
                    // If there are free threads to be created, we execute the next recursive call concurrently.
                    calculateForce(tree->children[s], sharedBuff, localBuff, indexes[i]);
                }
            }
            //Calculate new position
            moveParticle(sharedBuff,localBuff,indexes[i]);
        }
        int ret = pthread_barrier_wait(&calculateForceBarrier);
        // Block the semaphore again.
        if (ret != 0 && ret != PTHREAD_BARRIER_SERIAL_THREAD) {
            perror("Error in the barrier: ");
            cancelThreads(threads, numberOfThreads-1);
        }
        if (ret == PTHREAD_BARRIER_SERIAL_THREAD) {
            printf("Thread %d is unlocking the main thread for doing another iteration\n", id);
            // Barrier released, we can signal the main thread that the calculation is done.
            sem_post(&endedCalculatingForce);
        }
        // If it is the last iteration, we can end the thread.
        if (globalStruct.lastIteration) {
            printf("Thread %d is ending\n", id);
            pthread_exit(NULL);
        }
    }
}


#ifdef D_GLFW_SUPPORT
void drawParticle(double *shrdBuff, double *radius, int index){
    glBegin(GL_TRIANGLE_FAN);
    int k;
    glVertex2f(shrdBuff[PX(index)],shrdBuff[PY(index)]);
    for(k=0;k<20;k++){
        float angle=(float) (k)/19*2*3.141592;
        glVertex2f(shrdBuff[PX(index)]+radius[index]*cos(angle),shrdBuff[PY(index)]+radius[index]*sin(angle));
    }
    glEnd();
}

void drawBarnesHutDivisions(struct Node *rootNode){
    if(!rootNode->external){
        glBegin(GL_LINES);
        glVertex2f(rootNode->GCX,rootNode->LLY);
        glVertex2f(rootNode->GCX,rootNode->TRY);
        glVertex2f(rootNode->LLX,rootNode->GCY);
        glVertex2f(rootNode->TRX,rootNode->GCY);
        glEnd();
        int i;
        for(i=0;i<4;i++){
            if(rootNode->children[i]!=NULL){
                drawBarnesHutDivisions(rootNode->children[i]);
            }
        }
    }
}
#endif

void SaveGalaxy(int count, int nShared, int *indexes, double *sharedBuff);
void SaveGalaxyFile(char *filename, int nShared, int *indexes, double *sharedBuff);

void SaveGalaxy(int count, int nShared, int *indexes, double *sharedBuff)
{
    char filename[100];
    sprintf(filename,"./res/galaxy_%dB_%di.out",nShared,count);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);
}

void SaveGalaxyFile(char *filename, int nShared, int *indexes, double *sharedBuff)
{
    int i;
    FILE *res = fopen(filename,"w");

    fprintf(res,"%d\n",nShared);
    for(i=0;i<nShared;i++){
        fprintf(res,"%d\t%e\t%e\t%e\n",indexes[i],sharedBuff[PX(indexes[i])],sharedBuff[PY(indexes[i])],sharedBuff[MASS(indexes[i])]);
    }
    fclose(res);
}

void ReadGalaxyFile(char *filename, int *nShared, int **indexes, double **sharedBuff)
{
    int i,ind;
    FILE *input;

    printf("Reading file %s\n",filename);

    input = fopen(filename,"r");
    if (input==NULL) {
        printf("Error opening file.\n");
        exit(1);
    }
    // Read number of bodies.
    if (fscanf(input,"%d\n",nShared)<1){
        printf("Error reading number of particles.\n");
        exit(1);
    }


    printf("Reading %d bodies\n",*nShared);

    // Reserve memory for indexes and particles.
    *indexes = (int*) malloc(sizeof(int)*(*nShared));
    *sharedBuff = (double *) malloc(sizeof(double)*(3*(*nShared)+1));

    for(i=0;i<(*nShared);i++){
        (*indexes)[i]=i;
    }

    for(i=0;i<(*nShared);i++){
        if (fscanf(input,"%d\t%le\t%le\t%le\n", &ind,&((*sharedBuff)[PX((*indexes)[i])]),&((*sharedBuff)[PY((*indexes)[i])]),&((*sharedBuff)[MASS((*indexes)[i])]))<4){
            printf("Error reading number of particles.\n");
            exit(1);
        }
        //printf("Body %d: (%le,%le) %le\n", ind,(*sharedBuff)[PX((*indexes)[i])],(*sharedBuff)[PY((*indexes)[i])],(*sharedBuff)[MASS((*indexes)[i])]);
    }

    fclose(input);
}

#define DSaveIntermediateState 1
#define DIntervalIntermediateState 100
#define DShowStatistics 1
#define DIntervalStatistics 1

clock_t StartTime, EndTime;
double TimeSpent;

void ShowWritePartialResults(int count,int nOriginal, int nShared, int *indexes, double *sharedBuff)
{
    if (DSaveIntermediateState && !(count % DIntervalIntermediateState))
        SaveGalaxy(count, nOriginal, indexes, sharedBuff);

    if (DShowStatistics && !(count % DIntervalStatistics))
    {
        int i=0;
        double CurrentTime;
        CurrentTime = clock();
        TimeSpent = (double)(CurrentTime - StartTime) / CLOCKS_PER_SEC;
        //Mins = (int)TimeSpent/60;
        //Secs = (TimeSpent-(Mins*60));
        printf("[%.3f] Iteration %d => %d Bodies (%d) \t(Body %d: (%le, %le) %le).\n",TimeSpent, count, nShared, nOriginal, i, sharedBuff[PX(indexes[i])],sharedBuff[PY(indexes[i])],sharedBuff[MASS(indexes[i])]);
    }
}

int main(int argc, char *argv[]){
    int nShared=500;
	int steps=100;
    int M=4;
    int iterationsToPrint=25;
    double *sharedBuff;
    double *radius;
    int *indexes, i;
    char filename[100];

    printf("NBody with %d arguments. ",argc);
    if(atoi(argv[4]) == 0) {
        printf("Graphics OFF. ");
    } else {
        printf("Graphics ON. ");
    }
    StartTime = clock();


	if(argc>1){
		nShared=atoi(argv[1]);
		if(argc>2){
		  steps=atoi(argv[2]);
		} if(argc>5){
            M=atoi(argv[5]);
        } if(argc>6){
            iterationsToPrint=atoi(argv[6]);
        }
	}
    M++; // Count the main thread

    printf("Execution with %d threads.\n",M);

    // Initialize global variables
    numberOfThreads = M;
    globalStruct.lastIteration = false;
    // Init the calculateForcesBarrier
    pthread_barrier_init(&calculateForceBarrier, NULL, M-1);
    // Init semaphores
    sem_init(&calculateForceSem, 0, 0);
    sem_init(&endedCalculatingForce, 0, 0);

    // Start the threads
    threads = malloc(sizeof(pthread_t) * M-1);
    int t;
    for(t=0;t<M-1;t++) {
        if(pthread_create(&threads[t], NULL, (void *(*)(void *)) threadFunction, (void *) (t+1))) {
            perror("Error creating the thread for calculating the particles' forces: ");
            cancelThreads(threads, t);
        }
    }

    if(argc>3 && access(argv[3], F_OK) == 0)
    {
        printf("Read file..\n");
        /* Read bodies initial state from file */
        ReadGalaxyFile(argv[3], &nShared, &indexes, &sharedBuff);
        argc--;
    }
    else
    {   /* Inicialize the bodies randomly */

        //Buffers to hold the position of the particles and their mass
        sharedBuff = (double *) malloc(sizeof(double) * (3 * nShared + 1));

        srand(time(NULL));
        for (i = 0; i < nShared; i++) {
            //I start with an almost random distribution of particles
            sharedBuff[PX(i)] = (float) (i) / (nShared - 1) * 0.8 + 0.1;
            sharedBuff[PY(i)] = (float) (rand() % 4096) / 4095 * 0.8 + 0.1;
            //With a random Mass between 1 and 3
            sharedBuff[MASS(i)]=(double) (rand()%2048)/2047*2+1;
        }

        //Index array, to speed up the creation of the tree (faster than passing the 3 floats per particle of x,y and mass)
        indexes = (int*) malloc(sizeof(int)*nShared);
        for(i=0;i<nShared;i++){
            indexes[i]=i;
        }
    }

    int nLocal=nShared;
    int nOriginal=nShared;
    globalStruct.nLocal = nLocal;
    //Buffer to hold velocity in x and y, and acceleration in x and y also
    double* localBuff = (double *) malloc(sizeof(double)*(4*nLocal));
    //This is for opengl
    radius = (double *) malloc(sizeof(double)*(nShared));

    for(i=0;i<nShared;i++){
        // init bodies mass
        radius[i]=sqrt(sharedBuff[MASS(i)])*0.0025;

        //With zero speed, and zero acceleration
        localBuff[VX(i)]=0;
        localBuff[VY(i)]=0;
        localBuff[AX(i)]=0;
        localBuff[AY(i)]=0;
    }

    //This is the main node, the one that holds the first four children nodes that make the calculation zone
    struct Node* tree = malloc(sizeof *tree);
	//LLX is the x coordinate of the Low Left corner
    tree->LLX=0;
	//This is the y coordinate..
    tree->LLY=0;

	//Now the same but for the top right corner
    tree->TRX=1;
    tree->TRY=1;
	//The coordinates of the geometric center of the node in x and y
    tree->GCX=0.5;
    tree->GCY=0.5;

    // Save initial state.
    sprintf(filename,"./res/galaxy_%dB_initial.out",nOriginal);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);

    int count=1;
	//If we need to visualize
#ifdef D_GLFW_SUPPORT
	if(atoi(argv[4])!=0){
        //GRAPHIC MODE
	} else {
#endif
		//This is the pure algorithm, without visualization
		//system("mkdir res");
    	while(count<=steps){
			//First we build the tree
            // Create a BuildTreeStruct for passing the params to the buildTreeThread function
            struct BuildTreeStruct* data = malloc(sizeof(struct BuildTreeStruct));
            data->tree = tree;
            data->sharedBuff = sharedBuff;
            data->indexes = indexes;
            data->nShared = nShared;
            data->remainingThreads = M;
            buildTreeThread(data);

            // Calculate the indexes of the particles to calculate its forces by the main thread
            int totalThreads = M;
            int particlesPerThread = nLocal / totalThreads;

            pthread_mutex_lock(&remainingParticlesMutex);
            remainingParticles = nLocal % totalThreads;

            if (remainingParticles > 0) {
                particlesPerThread++;
                remainingParticles--;
            }
            pthread_mutex_unlock(&remainingParticlesMutex);

            int start = 0;
            int end = particlesPerThread;

            printf("Main thread calculating forces for particles %d to %d\n", start, end);

            // Assign values to the global struct
            globalStruct.localBuff = localBuff;
            globalStruct.indexes = indexes;
            globalStruct.sharedBuff = sharedBuff;
            globalStruct.tree = tree;
            globalStruct.nShared = nShared;

            if (count==steps) {
                globalStruct.lastIteration = true;
            }

            // Unlock the semaphore to start the threads calculating the forces
            for (i = 0; i < M-1; i++) {
                sem_post(&calculateForceSem);
            }

        	for(i=start;i<end;i++){
				//Set initial accelerations to zero
            	localBuff[AX(indexes[i])]=0;
            	localBuff[AY(indexes[i])]=0;
            	int s;

            	for(s=0;s<4;s++){
					//Recursively calculate accelerations
                	if(tree->children[s]!=NULL){
                        // If there are free threads to be created, we execute the next recursive call concurrently.
                        calculateForce(tree->children[s], sharedBuff, localBuff, indexes[i]);
                    }
            	}
				//Calculate new position
            	moveParticle(sharedBuff,localBuff,indexes[i]);
        	}

            // Wait for all the threads to finish its work
            sem_wait(&endedCalculatingForce);
            // Then check if there are particles outside the boundaries, and remove them
            printf("All threads finished calculating forces, main thread will process the eliminated particles and continue to the next iteration\n");
            for (i = 0; i < globalStruct.nLocal; ++i) {
                //Kick out particle if it went out of the box (0,1)x(0,1)
                if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                    int r;
                    globalStruct.nLocal--;
                    nShared--;
                    for(r=i;r<globalStruct.nLocal;r++){
                        indexes[r]=indexes[r+1];
                    }
                    i--;
                }
            }

			//To be able to store the positions of the particles
            ShowWritePartialResults(count,nOriginal, nShared, indexes, sharedBuff);
            //We advance one step
			count++;
            free(data);
		}
#ifdef D_GLFW_SUPPORT
	}
#endif
    // Join the threads
    for (i = 0; i < M-1; i++) {
        pthread_join(threads[i], NULL);
    }

    // Free the memory and destroy the mutexes and semaphores
    free(threads);
    destroySyncVariables();

    EndTime = clock();
    TimeSpent = (double)(EndTime - StartTime) / CLOCKS_PER_SEC;
    printf("NBody Simulation took %.3f seconds.\n",TimeSpent);

    // Save initial state.
    sprintf(filename,"./res/galaxy_%dB_%di_final.out",nOriginal, count-1);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);

	free(sharedBuff);
	free(localBuff);
	free(radius);
	free(indexes);

    return 0;
}
