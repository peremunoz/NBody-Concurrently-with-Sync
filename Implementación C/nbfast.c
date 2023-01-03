/* ---------------------------------------------------------------
Práctica 3.
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

struct BuildTreeStruct{
    struct Node* tree;
    double* sharedBuff;
    int* indexes;
    int nShared;
    int remainingThreads;
};

struct GlobalStruct{
    int nLocal;
    double* localBuff;
    int *indexes;
    double *sharedBuff;
    struct Node *tree;
    int actualIteration;
    bool lastIteration;
};

struct UIStruct{
    struct Node *tree;
    int nShared;
    double *sharedBuff;
    double *radius;
    int *indexes;
};

struct Statistics{
    double timePerIteration;
    double timePerMIterations;
    int evaluatedParticles;
    int removedParticles;
    int numberOfSimplifications;
};

// Global synchronization variables
sem_t calculateForceSem;
pthread_mutex_t remainingParticlesMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t calculateForceBarrier;
sem_t endedCalculatingForce;
pthread_mutex_t statisticsMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t visualizationCond = PTHREAD_COND_INITIALIZER;
pthread_mutex_t visualizationMutex = PTHREAD_MUTEX_INITIALIZER;

// Global variables
pthread_t *threads;
pthread_t uiThread;
int numberOfThreads;
int iterationsToPrint;
int remainingParticles;
struct GlobalStruct globalStruct;
struct Statistics* threadsStatistics;
struct Statistics globalStatistics;
struct UIStruct uiStruct;


void cancelThreads(pthread_t* tids, int nThreads){
    for(int i=0;i<nThreads;i++){
        pthread_cancel(tids[i]);
    }
    if(uiThread != 0){
        pthread_cancel(uiThread);
    }
    exit(-1);
}

void destroySyncVariables(){
    sem_destroy(&calculateForceSem);
    pthread_mutex_destroy(&remainingParticlesMutex);
    pthread_barrier_destroy(&calculateForceBarrier);
    sem_destroy(&endedCalculatingForce);
    pthread_mutex_destroy(&statisticsMutex);
    pthread_cond_destroy(&visualizationCond);
    pthread_mutex_destroy(&visualizationMutex);
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

int calculateForce(struct Node *tree, double *shrdBuff, double *localBuff, int index){
    double distance = sqrt((tree->CMX-shrdBuff[PX(index)])*(tree->CMX-shrdBuff[PX(index)])+
                           (tree->CMY-shrdBuff[PY(index)])*(tree->CMY-shrdBuff[PY(index)]));
    int simplifications=0;
    //First we check if the node is not actually the same particle we are calculating
    if(distance>0){
        //Now, we know it is not because the is some distance between the Center of Mass and the particle
        //If the node is external (only contains one particle) or is far away enough, we calculate the force with the center of mass
        if(distance>rcutoff || tree->external){
            double f;
            if(distance<rlimit){
                f=G*tree->mass/(rlimit*rlimit*distance);
            } else {
                simplifications++;
                f=G*tree->mass/(distance*distance*distance);
            }
            localBuff[AX(index)]+=f*(tree->CMX-shrdBuff[PX(index)]);
            localBuff[AY(index)]+=f*(tree->CMY-shrdBuff[PY(index)]);
        } else {
            //If not, we recursively call the calculateForce() function in the children that are not empty.
            int i;
            for(i=0;i<4;i++){
                if(tree->children[i]!=NULL){
                    return simplifications + calculateForce(tree->children[i],shrdBuff,localBuff,index);
                }
            }
        }
    }
    return 0;
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
        // Wait for the main thread to signal that the tree is built and the particles are initialized.
        sem_wait(&calculateForceSem);

        if (globalStruct.actualIteration > 0 && globalStruct.actualIteration % iterationsToPrint == 1) {
            pthread_mutex_lock(&statisticsMutex);
            threadsStatistics[id].timePerMIterations = 0;
            pthread_mutex_unlock(&statisticsMutex);
        }

        // When unlocked, init the timer for the calculation of the statistics
        struct timespec startTime, endTime;
        clock_gettime(CLOCK_MONOTONIC, &startTime);

        // Calculate index of the particles to be calculated.
        int particlesPerThread = globalStruct.nLocal / numberOfThreads;

        int start = id * particlesPerThread;
        int end = start + particlesPerThread;

        if (id == numberOfThreads - 1 && remainingParticles > 0) {
            end += remainingParticles;
            remainingParticles = 0;
        }

        if (end > globalStruct.nLocal) {
            end = globalStruct.nLocal;
        }

        printf("Thread %d: start = %d, end = %d\n", id, start, end);

        double *localBuff = globalStruct.localBuff;
        int *indexes = globalStruct.indexes;
        double *sharedBuff = globalStruct.sharedBuff;
        struct Node* tree = globalStruct.tree;

        // Then, we can calculate the forces.
        int i;
        int removedParticles = 0;
        int simplifications = 0;
        for(i=start; i < end; i++){
            //Set initial accelerations to zero
            localBuff[AX(indexes[i])]=0;
            localBuff[AY(indexes[i])]=0;
            int s;

            for(s=0;s<4;s++){
                //Recursively calculate accelerations
                if(tree->children[s]!=NULL){
                    // If there are free threads to be created, we execute the next recursive call concurrently.
                    simplifications += calculateForce(tree->children[s], sharedBuff, localBuff, indexes[i]);
                }
            }
            //Calculate new position
            moveParticle(sharedBuff,localBuff,indexes[i]);

            if (sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1) {
                // If the particle is out of the limits, we count it for the statistics.
                removedParticles++;
            }
        }

        // Now that the calculation is done, we register the time it took to calculate the forces.
        clock_gettime(CLOCK_MONOTONIC, &endTime);
        double timeElapsed = (endTime.tv_sec - startTime.tv_sec);
        timeElapsed += (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;

        // Update the statistics
        pthread_mutex_lock(&statisticsMutex);
        threadsStatistics[id].timePerIteration += timeElapsed;
        threadsStatistics[id].timePerMIterations += timeElapsed;
        threadsStatistics[id].evaluatedParticles += end - start;
        threadsStatistics[id].removedParticles += removedParticles;
        threadsStatistics[id].numberOfSimplifications += simplifications;
        globalStatistics.numberOfSimplifications += simplifications;
        pthread_mutex_unlock(&statisticsMutex);

        // Barrier for waiting for all threads to finish the calculation of the forces.
        int ret = pthread_barrier_wait(&calculateForceBarrier);

        // Once all threads have finished its calculation and updated its statistics, we can print the partial results if needed.
        if (globalStruct.actualIteration > 0 && globalStruct.actualIteration % iterationsToPrint == 0) {
            double averageTimePerMIterations = 0;
            pthread_mutex_lock(&statisticsMutex);
            for (i = 0; i < numberOfThreads; i++) {
                averageTimePerMIterations += threadsStatistics[i].timePerMIterations;
            }
            averageTimePerMIterations /= numberOfThreads;
            pthread_mutex_unlock(&statisticsMutex);
            printf("Thread [%d]: Average of %f seconds per iteration, %f seconds of unbalance %i evaluated particles, %i removed particles, %i of simplifications done.\n", id, threadsStatistics[id].timePerIteration/globalStruct.actualIteration, threadsStatistics[id].timePerMIterations-averageTimePerMIterations, threadsStatistics[id].evaluatedParticles, threadsStatistics[id].removedParticles, threadsStatistics[id].numberOfSimplifications);
        }

        if (ret != 0 && ret != PTHREAD_BARRIER_SERIAL_THREAD) {
            perror("Error in the barrier: ");
            cancelThreads(threads, numberOfThreads-1);
        }
        if (ret == PTHREAD_BARRIER_SERIAL_THREAD) {
            // Barrier released, we can signal the main thread that the calculation is done.
            sem_post(&endedCalculatingForce);
        }
        // If it is the last iteration, we can end the thread.
        if (globalStruct.lastIteration) {
            printf("[%d] is ending\n", id);
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

#ifdef D_GLFW_SUPPORT
void uiThreadFunction() {
    if(!glfwInit()) {
        printf("Failed to start GLFW\n");
        return;
    }
    GLFWwindow *window = glfwCreateWindow(2000,2000,"Simulation",NULL,NULL);
    if(!window){
        printf("Failed to open window\n");
        return;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,1,0,1,0,1);
    glMatrixMode(GL_MODELVIEW);

    // Wait for receiving the condition to refresh the window.
    while(1) {
        pthread_mutex_lock(&visualizationMutex);
        pthread_cond_wait(&visualizationCond, &visualizationMutex);
        pthread_mutex_unlock(&visualizationMutex);
        glClear(GL_COLOR_BUFFER_BIT);
        double t=glfwGetTime();

        // If the condition is met, we can refresh the window.
        drawBarnesHutDivisions(uiStruct.tree);
        int k;
        for(k=0;k<uiStruct.nShared;k++){
            drawParticle(uiStruct.sharedBuff,uiStruct.radius,uiStruct.indexes[k]);
        }
        t=glfwGetTime()-t;
        if(t<0.013){
            usleep(1000*1000*(0.013-t));
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (glfwWindowShouldClose(window) || globalStruct.lastIteration) {
            break;
        }
    }
    glfwTerminate();
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
    iterationsToPrint=25;
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

    printf("Execution with %d threads (%d + the main one).\n",M, M-1);

    // Initialize global variables
    numberOfThreads = M;
    globalStruct.lastIteration = false;

    // Init the threads statistics array
    threadsStatistics = (struct Statistics *) malloc(sizeof(struct Statistics) * M);

    // Init the calculateForcesBarrier
    pthread_barrier_init(&calculateForceBarrier, NULL, M-1);

    // Init semaphores
    sem_init(&calculateForceSem, 0, 0);
    sem_init(&endedCalculatingForce, 0, 0);

    // Start the threads
    threads = malloc(sizeof(pthread_t) * M-1);
    int t;
    for(t=0;t<M-1;t++) {
        if(pthread_create(&threads[t], NULL, (void *(*)(void *)) threadFunction, (void *) (size_t) (t+1))) {
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
	//This is the y coordinate.
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
        // Start the visualization thread
        if(pthread_create(&uiThread, NULL, (void *(*)(void *)) uiThreadFunction, NULL)) {
            perror("Error creating the thread for visualizing the particles: ");
            cancelThreads(threads, numberOfThreads-1);
        }

        // Alloc memory for the ui struct pointers
        uiStruct.sharedBuff = (double *) malloc(sizeof(double) * (3 * nShared + 1));
        uiStruct.indexes = (int *) malloc(sizeof(int) * nShared);
        uiStruct.radius = (double *) malloc(sizeof(double) * nShared);
        uiStruct.tree = malloc(sizeof(struct Node));

        while(count<=steps){
            // First we build the tree
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
            int particlesPerThread = globalStruct.nLocal / totalThreads;

            // For preventing miscalculation of the average time for all threads, we update this time per M iterations variable here, in the M'th + 1 iteration
            if (globalStruct.actualIteration > 0 && globalStruct.actualIteration % iterationsToPrint == 1) {
                pthread_mutex_lock(&statisticsMutex);
                threadsStatistics[0].timePerMIterations = 0;
                pthread_mutex_unlock(&statisticsMutex);
            }

            int start = 0;
            int end = particlesPerThread;

            printf("[Main thread] calculating forces for particles %d to %d\n", start, end);

            // Assign values to the global struct
            globalStruct.localBuff = localBuff;
            globalStruct.indexes = indexes;
            globalStruct.sharedBuff = sharedBuff;
            globalStruct.tree = tree;

            if (count==steps) {
                globalStruct.lastIteration = true;
            }
            globalStruct.actualIteration = count;

            globalStatistics.evaluatedParticles += globalStruct.nLocal;

            // Unlock the semaphore to start the threads calculating the forces
            for (i = 0; i < M-1; i++) {
                sem_post(&calculateForceSem);
            }

            // Start the main thread clock
            struct timespec startTime, endTime;
            clock_gettime(CLOCK_MONOTONIC, &startTime);

            int removedParticles = 0;

            int simplifications = 0;

            for(i=start;i<end;i++){
                //Set initial accelerations to zero
                localBuff[AX(indexes[i])]=0;
                localBuff[AY(indexes[i])]=0;
                int s;

                for(s=0;s<4;s++){
                    //Recursively calculate accelerations
                    if(tree->children[s]!=NULL){
                        // If there are free threads to be created, we execute the next recursive call concurrently.
                        simplifications += calculateForce(tree->children[s], sharedBuff, localBuff, indexes[i]);
                    }
                }
                //Calculate new position
                moveParticle(sharedBuff,localBuff,indexes[i]);

                if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                    removedParticles++;
                }
            }

            // Wait for all the threads to finish its work
            sem_wait(&endedCalculatingForce);

            // Then check if there are particles outside the boundaries, and remove them
            for (i = 0; i < globalStruct.nLocal; ++i) {
                //Kick out particle if it went out of the box (0,1)x(0,1)
                if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                    globalStatistics.removedParticles++;
                    int r;
                    globalStruct.nLocal--;
                    nShared--;
                    for(r=i;r<globalStruct.nLocal;r++){
                        indexes[r]=indexes[r+1];
                    }
                    i--;
                }
            }

            SaveGalaxy(count, nShared, indexes, sharedBuff);

            // Update the ui struct copying the actual iteration results to it
            uiStruct.nShared = nShared;
            memcpy(uiStruct.sharedBuff, sharedBuff, sizeof(double) * (3 * nShared + 1));
            memcpy(uiStruct.indexes, indexes, sizeof(int) * nShared);
            memcpy(uiStruct.radius, radius, sizeof(double) * nShared);
            *uiStruct.tree = *tree;

            // Communicate the visualization thread that the particles have been updated and it can draw them
            pthread_cond_signal(&visualizationCond);

            // Stop the main thread clock
            clock_gettime(CLOCK_MONOTONIC, &endTime);
            double timeElapsed = (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;

            // Update the main thread statistics
            pthread_mutex_lock(&statisticsMutex);
            threadsStatistics[0].timePerIteration += timeElapsed;
            threadsStatistics[0].timePerMIterations += timeElapsed;
            threadsStatistics[0].evaluatedParticles += end - start;
            threadsStatistics[0].removedParticles += removedParticles;
            threadsStatistics[0].numberOfSimplifications += simplifications;
            globalStatistics.numberOfSimplifications += simplifications;
            pthread_mutex_unlock(&statisticsMutex);

            // When all the threads have finished, we can print the main thread statistics
            if (globalStruct.actualIteration > 0 && globalStruct.actualIteration % iterationsToPrint == 0) {
                double averageTimePerMIterations = 0;
                pthread_mutex_lock(&statisticsMutex);
                for (i = 0; i < numberOfThreads; i++) {
                    averageTimePerMIterations += threadsStatistics[i].timePerMIterations;
                }
                averageTimePerMIterations /= numberOfThreads;
                pthread_mutex_unlock(&statisticsMutex);
                printf("Thread [%d]: Average of %f seconds per iteration, %f seconds of unbalance, %i evaluated particles, %i removed particles, %i of simplifications done.\n", 0, threadsStatistics[0].timePerIteration/count, threadsStatistics[0].timePerMIterations - averageTimePerMIterations, threadsStatistics[0].evaluatedParticles, threadsStatistics[0].removedParticles, threadsStatistics[0].numberOfSimplifications);

                // We print also the global statistics
                printf("Global statistics: %i evaluated particles, %i removed particles, %i of simplifications done.\n", globalStatistics.evaluatedParticles, globalStatistics.removedParticles, globalStatistics.numberOfSimplifications);
            }

            //To be able to store the positions of the particles
            ShowWritePartialResults(count,nOriginal, nShared, indexes, sharedBuff);
            //We advance one step
            count++;
            free(data);
        }
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
            int particlesPerThread = globalStruct.nLocal / totalThreads;

            // For preventing miscalculation of the average time for all threads, we update this time per M iterations variable here, in the M'th + 1 iteration
            if (globalStruct.actualIteration > 0 && globalStruct.actualIteration % iterationsToPrint == 1) {
                pthread_mutex_lock(&statisticsMutex);
                threadsStatistics[0].timePerMIterations = 0;
                pthread_mutex_unlock(&statisticsMutex);
            }

            remainingParticles = globalStruct.nLocal % totalThreads;

            int start = 0;
            int end = particlesPerThread;

            printf("[Main thread] calculating forces for particles %d to %d\n", start, end);

            // Assign values to the global struct
            globalStruct.localBuff = localBuff;
            globalStruct.indexes = indexes;
            globalStruct.sharedBuff = sharedBuff;
            globalStruct.tree = tree;

            if (count==steps) {
                globalStruct.lastIteration = true;
            }
            globalStruct.actualIteration = count;

            globalStatistics.evaluatedParticles += globalStruct.nLocal;

            // Unlock the semaphore to start the threads calculating the forces
            for (i = 0; i < M-1; i++) {
                sem_post(&calculateForceSem);
            }

            // Start the main thread clock
            struct timespec startTime, endTime;
            clock_gettime(CLOCK_MONOTONIC, &startTime);

            int removedParticles = 0;

            int simplifications = 0;

        	for(i=start;i<end;i++){
				//Set initial accelerations to zero
            	localBuff[AX(indexes[i])]=0;
            	localBuff[AY(indexes[i])]=0;
            	int s;

            	for(s=0;s<4;s++){
					//Recursively calculate accelerations
                	if(tree->children[s]!=NULL){
                        // If there are free threads to be created, we execute the next recursive call concurrently.
                        simplifications += calculateForce(tree->children[s], sharedBuff, localBuff, indexes[i]);
                    }
            	}
				//Calculate new position
            	moveParticle(sharedBuff,localBuff,indexes[i]);

                if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                    removedParticles++;
                }
        	}

            // Wait for all the threads to finish its work
            sem_wait(&endedCalculatingForce);

            // Then check if there are particles outside the boundaries, and remove them
            for (i = 0; i < globalStruct.nLocal; ++i) {
                //Kick out particle if it went out of the box (0,1)x(0,1)
                if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                    globalStatistics.removedParticles++;
                    int r;
                    globalStruct.nLocal--;
                    nShared--;
                    for(r=i;r<globalStruct.nLocal;r++){
                        indexes[r]=indexes[r+1];
                    }
                    i--;
                }
            }

            // Stop the main thread clock
            clock_gettime(CLOCK_MONOTONIC, &endTime);
            double timeElapsed = (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;

            // Update the main thread statistics
            pthread_mutex_lock(&statisticsMutex);
            threadsStatistics[0].timePerIteration += timeElapsed;
            threadsStatistics[0].timePerMIterations += timeElapsed;
            threadsStatistics[0].evaluatedParticles += end - start;
            threadsStatistics[0].removedParticles += removedParticles;
            threadsStatistics[0].numberOfSimplifications += simplifications;
            globalStatistics.numberOfSimplifications += simplifications;
            pthread_mutex_unlock(&statisticsMutex);

            // When all the threads have finished, we can print the main thread statistics
            if (globalStruct.actualIteration > 0 && globalStruct.actualIteration % iterationsToPrint == 0) {
                double averageTimePerMIterations = 0;
                pthread_mutex_lock(&statisticsMutex);
                for (i = 0; i < numberOfThreads; i++) {
                    averageTimePerMIterations += threadsStatistics[i].timePerMIterations;
                }
                averageTimePerMIterations /= numberOfThreads;
                pthread_mutex_unlock(&statisticsMutex);
                printf("Thread [%d]: Average of %f seconds per iteration, %f seconds of unbalance, %i evaluated particles, %i removed particles, %i of simplifications done.\n", 0, threadsStatistics[0].timePerIteration/count, threadsStatistics[0].timePerMIterations - averageTimePerMIterations, threadsStatistics[0].evaluatedParticles, threadsStatistics[0].removedParticles, threadsStatistics[0].numberOfSimplifications);

                // We print also the global statistics
                printf("Global statistics: %i evaluated particles, %i removed particles, %i of simplifications done.\n", globalStatistics.evaluatedParticles, globalStatistics.removedParticles, globalStatistics.numberOfSimplifications);
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
        printf("[%d] joined\n", i+1);
    }
    if(atoi(argv[4])!=0){
        pthread_join(uiThread, NULL);
        printf("[UI Thread] joined\n");
    }

    EndTime = clock();
    TimeSpent = (double)(EndTime - StartTime) / CLOCKS_PER_SEC;

    // Print the global statistics
    printf("Global statistics:\n");
    printf("Total simulation time: %f seconds, %i evaluated particles, %i removed particles, %i of simplifications done.\n", TimeSpent, globalStatistics.evaluatedParticles, globalStatistics.removedParticles, globalStatistics.numberOfSimplifications);


    // Print the threads statistics
    printf("Threads statistics:\n");

    double averageTimeOfAllThreads = 0;

    for (i = 0; i < numberOfThreads; i++) {
        averageTimeOfAllThreads += threadsStatistics[i].timePerIteration;
    }
    averageTimeOfAllThreads /= numberOfThreads;
    for (i = 0; i < M; i++) {
        printf("Thread [%d]: %f seconds of total work, %f seconds per iteration, %f seconds of total unbalance, %i evaluated particles, %i removed particles, %i of simplifications done.\n", i, threadsStatistics[i].timePerIteration, threadsStatistics[i].timePerIteration/steps, threadsStatistics[i].timePerIteration-averageTimeOfAllThreads,threadsStatistics[i].evaluatedParticles, threadsStatistics[i].removedParticles, threadsStatistics[i].numberOfSimplifications);
    }

    // Free the memory and destroy the mutexes and semaphores
    free(threadsStatistics);
    free(threads);
    destroySyncVariables();

    // Save initial state.
    sprintf(filename,"./res/galaxy_%dB_%di_final.out",nOriginal, count-1);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);

	free(sharedBuff);
	free(localBuff);
	free(radius);
	free(indexes);

    return 0;
}