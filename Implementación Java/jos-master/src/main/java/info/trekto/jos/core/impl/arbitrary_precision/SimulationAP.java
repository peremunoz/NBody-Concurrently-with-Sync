/*
Práctica 4.
Código fuente: SimulationAP.java
Grau Informàtica
48252062V - Pere Muñoz Figuerol
*/

package info.trekto.jos.core.impl.arbitrary_precision;

import info.trekto.jos.core.Controller;
import info.trekto.jos.core.ForceCalculator;
import info.trekto.jos.core.Simulation;
import info.trekto.jos.core.exceptions.SimulationException;
import info.trekto.jos.core.impl.Iteration;
import info.trekto.jos.core.impl.SimulationProperties;
import info.trekto.jos.core.model.ImmutableSimulationObject;
import info.trekto.jos.core.model.SimulationObject;
import info.trekto.jos.core.model.impl.SimulationObjectImpl;
import info.trekto.jos.core.numbers.Number;
import info.trekto.jos.util.Utils;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.time.StopWatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.time.Clock;
import java.util.*;

import static info.trekto.jos.core.Controller.C;
// import static info.trekto.jos.core.impl.arbitrary_precision.SimulationRecursiveAction.threshold;
import static info.trekto.jos.core.Controller.main;
import static info.trekto.jos.util.Utils.*;

/**
 * This implementation uses fork/join Java framework introduced in Java 7.
 *
 * @author Trayan Momkov
 * 2017-May-18
 */
public class SimulationAP implements Simulation {
    private static final Logger logger = LoggerFactory.getLogger(SimulationAP.class);
    public static final int PAUSE_SLEEP_MILLISECONDS = 100;
    public static final int SHOW_REMAINING_INTERVAL_SECONDS = 2;

    private final SimulationLogicAP simulationLogic;
    protected SimulationProperties properties;
    private ForceCalculator forceCalculator;
    protected long iterationCounter;

    private List<SimulationObject> objects;
    private List<SimulationObject> auxiliaryObjects;
    
    private volatile boolean collisionExists;

    // Define a thread list for the threads that will be created
    private static List<Thread> threads;

    // Define a List for storing the statistics of each thread
    private static List<Statistics> threadsStatistics;
    private static Statistics mainThreadStatistics;

    // Define a Statistics object for storing the global statistics
    private static Statistics globalStatistics;

    // Get M from the arguments
    private static int M;

    public SimulationAP(SimulationProperties properties) {
        threads = new ArrayList<>();
        simulationLogic = new SimulationLogicAP(this);
        this.properties = properties;

        // Initialize the threadsStatistics list
        threadsStatistics = new LinkedList<>();
        for (int i = 0; i < properties.getNumberOfThreads()+1; i++) {
            threadsStatistics.add(new Statistics());
        }

        mainThreadStatistics = new Statistics();

        threadsStatistics.set(0, mainThreadStatistics);

        // Initialize the threads and start them
        for (int i = 0; i < properties.getNumberOfThreads(); i++) {
            int finalI = i;
            threads.add(new Thread(() -> {
                try {
                    simulationLogic.threadFunction(finalI + 1);
                } catch (Exception e) {
                    cancelThreads();
                    e.printStackTrace();
                }
            }));
        }

        for (Thread thread : threads) {
            thread.start();
        }

        // Initialize the globalStatistics object
        globalStatistics = new Statistics();

        // Get M from the program arguments
        M = Integer.parseInt(Controller.getArgs()[0]);
    }

    public static void cancelThreads() {
        for (Thread thread : threads) {
            thread.interrupt();
        }
        System.exit(-1);
    }

    public static Statistics getStatistics() {
        return globalStatistics;
    }

    public static Statistics getMainThreadStatistics() {
        return mainThreadStatistics;
    }

    public static List<Statistics> getThreadsStatistics() {
        return threadsStatistics;
    }

    public static int getM() {
        return M;
    }

    public boolean collisionExists(List<SimulationObject> objects) {
        for (SimulationObject o : objects) {
            for (SimulationObject o1 : objects) {
                if (o == o1) {
                    continue;
                }
                // distance between centres
                Number distance = simulationLogic.calculateDistance(o, o1);

                if (distance.compareTo(o.getRadius().add(o1.getRadius())) < 0) {
                    info(logger, String.format("Collision between object A(x:%f, y:%f, r:%f) and B(x:%f, y:%f, r:%f)",
                                               o.getX().doubleValue(), o.getY().doubleValue(), o.getRadius().doubleValue(),
                                               o1.getX().doubleValue(), o1.getY().doubleValue(), o1.getRadius().doubleValue()));
                    return true;
                }
            }
        }
        return false;
    }

    public boolean duplicateIdExists(List<SimulationObject> objects) {
        Set<String> ids = new HashSet<>();
        for (SimulationObject object : objects) {
            if (!ids.add(object.getId())) {
                return true;
            }
        }
        return false;
    }

    public void doIteration(boolean saveCurrentIterationToFile, long iterationCounter) throws InterruptedException {
        auxiliaryObjects = deepCopy(objects);

        // Start the clock for the main thread statistics
        float startTime = System.nanoTime();

        simulationLogic.calculateAllNewValues();

        System.out.println("Iteration " + iterationCounter + " with " + objects.size() + " objects");

        // Update the global statistics
        globalStatistics.addEvaluatedParticles(objects.size());


        // Wait for the threads to finish their calculations
        float startWaitingTime = System.nanoTime();
        simulationLogic.getCalculateNewValuesSemaphore().acquire(properties.getNumberOfThreads());
        float elapsedWaitingTime = System.nanoTime() - startWaitingTime;

        /* Collision */

        // Calculate the main thread collisions' indexes
        int start = 0;
        int end = objects.size() / (properties.getNumberOfThreads()+1);
        CollisionCheckAP collisionCheck = new CollisionCheckAP(start, end, this);
        collisionExists = false;

        System.out.println("[MAIN THREAD] Start: " + start + " End: " + end);

        // Signal the threads to start checking for collisions
        simulationLogic.getCheckCollisionsSemaphore().release(properties.getNumberOfThreads());

        // Start the main thread's collision check
        synchronized (SimulationAP.class) {
            collisionCheck.checkAllCollisions();
        }
        System.out.println("[MAIN THREAD] Ended checking collisions");

        // Wait for the threads to finish checking for collisions
        float startWaitingTime2 = System.nanoTime();
        simulationLogic.getEndedCheckingCollisionsSemaphore().acquire(properties.getNumberOfThreads());

        elapsedWaitingTime += System.nanoTime() - startWaitingTime2;

        /* If collision/s exists execute sequentially on a single thread */
        if (collisionExists) {
            simulationLogic.processCollisions(this);
        }

        objects = auxiliaryObjects;

        // Stop the clock for the main thread statistics
        float elapsedTime = System.nanoTime() - startTime - elapsedWaitingTime;

        // Update the main thread statistics
        mainThreadStatistics.addTime(elapsedTime);

        // Print the main thread statistics if needed
        if (iterationCounter % M == 0 && iterationCounter != 0 && iterationCounter != properties.getNumberOfIterations()) {

            // Now signal the threads that they can continue because all collisions have been checked
            simulationLogic.getContinueToPrintStatisticsSemaphore().release(properties.getNumberOfThreads());

            // Wait for the threads to finish printing their statistics
            simulationLogic.getPrintStatisticsSemaphore().acquire(properties.getNumberOfThreads());

            // Print the statistics
            simulationLogic.getStatisticsLock().lock();
            System.out.println("\n[Thread 0] Statistics for iteration " + iterationCounter + ": ");
            System.out.print("Time per iteration: " + elapsedTime/((float)iterationCounter*(float)1000000) + "ms, ");
            System.out.print("Total execution time: " + mainThreadStatistics.getTotalTime()/((float)1000000000) + "s, ");

            float averageMTime = 0;
            for (Statistics threadStatistic : threadsStatistics) {
                averageMTime += threadStatistic.getTimePerMIterations();
            }
            averageMTime /= threadsStatistics.size();
            System.out.print("Generated unbalance: " + (mainThreadStatistics.getTimePerMIterations()-averageMTime)/(1000000) + "ms, ");
            System.out.print("Evaluated particles: " + mainThreadStatistics.getEvaluatedParticles());
            System.out.println(" Merged particles: " + mainThreadStatistics.getMergedParticles());
            simulationLogic.getStatisticsLock().unlock();

            // Print the global statistics
            System.out.println("\nGlobal statistics for iteration " + iterationCounter + ": ");
            System.out.print("Evaluated particles: " + globalStatistics.getEvaluatedParticles());
            System.out.println(" Merged particles: " + globalStatistics.getMergedParticles());

            mainThreadStatistics.resetTimePerMIterations();
        }

        /* Slow down visualization */
        if (properties.isRealTimeVisualization() && properties.getPlayingSpeed() < 0) {
            Thread.sleep(-properties.getPlayingSpeed());
        }

        if (properties.isSaveToFile() && saveCurrentIterationToFile) {
            C.getReaderWriter().appendObjectsToFile(objects, properties, iterationCounter);
        }
    }

    public void playSimulation(String inputFile) {
        try {
            // Only reset reader pointer. Do not change properties! We want to have the latest changes from the GUI.
            C.getReaderWriter().readPropertiesForPlaying(inputFile);
        } catch (IOException e) {
            error(logger, "Cannot reset input file for playing.", e);
        }
        C.setVisualizer(C.createVisualizer(properties));
        long previousTime = System.nanoTime();
        long previousVisualizationTime = previousTime;
        C.setRunning(true);
        C.setEndText("END.");
        try {
            while (C.getReaderWriter().hasMoreIterations()) {
                if (C.hasToStop()) {
                    doStop();
                    break;
                }
                while (C.isPaused()) {
                    Thread.sleep(PAUSE_SLEEP_MILLISECONDS);
                }
                Iteration iteration = C.getReaderWriter().readNextIteration();
                if (iteration == null) {
                    break;
                }

                if (System.nanoTime() - previousTime >= NANOSECONDS_IN_ONE_SECOND * SHOW_REMAINING_INTERVAL_SECONDS) {
                    previousTime = System.nanoTime();
                    info(logger, "Cycle: " + iteration.getCycle() + ", number of objects: " + iteration.getNumberOfObjects());
                }

                if (properties.getPlayingSpeed() < 0) {
                    /* Slow down */
                    Thread.sleep(-properties.getPlayingSpeed());
                    C.getVisualizer().visualize(iteration);
                    previousVisualizationTime = System.nanoTime();
                } else if ((System.nanoTime() - previousVisualizationTime) / NANOSECONDS_IN_ONE_MILLISECOND >= properties.getPlayingSpeed()) {
                    C.getVisualizer().visualize(iteration);
                    previousVisualizationTime = System.nanoTime();
                }
            }
            info(logger, "End.");
            C.getVisualizer().end();
        } catch (IOException e) {
            error(logger, "Error while reading simulation object.", e);
        } catch (InterruptedException e) {
            error(logger, "Thread interrupted.", e);
        } finally {
            C.setRunning(false);
        }
    }

    protected void doStop() {
        C.setHasToStop(false);
        if (properties.isSaveToFile()) {
            C.getReaderWriter().endFile();
        }
        if (C.getVisualizer() != null) {
            C.setEndText("Stopped!");
            C.getVisualizer().closeWindow();
        }
    }

    @Override
    public void startSimulation() throws SimulationException {
        init(true);

        info(logger, "Start simulation...");
        C.setEndText("END.");
        long startTime = System.nanoTime();
        long previousTime = startTime;
        long previousVisualizationTime = startTime;
        long endTime;

        C.setRunning(true);
        C.setHasToStop(false);
        try {
            for (long i = 0; properties.isInfiniteSimulation() || i < properties.getNumberOfIterations(); i++) {
                try {
                    if (C.hasToStop()) {
                        doStop();
                        break;
                    }
                    while (C.isPaused()) {
                        Thread.sleep(PAUSE_SLEEP_MILLISECONDS);
                    }

                    iterationCounter = i + 1;

                    if (System.nanoTime() - previousTime >= NANOSECONDS_IN_ONE_SECOND * SHOW_REMAINING_INTERVAL_SECONDS) {
                        showRemainingTime(i, startTime, properties.getNumberOfIterations(), objects.size());
                        previousTime = System.nanoTime();
                    }

                    if (properties.isRealTimeVisualization()) {
                        endTime = System.nanoTime();
                        if (properties.getPlayingSpeed() < 0) {
                            /* Slow down */
                            Thread.sleep(-properties.getPlayingSpeed());
                            C.getVisualizer().visualize(objects, endTime-startTime, iterationCounter);
                            previousVisualizationTime = System.nanoTime();
                        } else if ((System.nanoTime() - previousVisualizationTime) / NANOSECONDS_IN_ONE_MILLISECOND >= properties.getPlayingSpeed()) {
                            C.getVisualizer().visualize(objects, endTime-startTime, iterationCounter);
                            previousVisualizationTime = System.nanoTime();
                        }
                    }

                    doIteration(i % properties.getSaveEveryNthIteration() == 0, iterationCounter);
                } catch (InterruptedException e) {
                    error(logger, "Concurrency failure. One of the threads interrupted in cycle " + i, e);
                }
            }

            if (properties.isRealTimeVisualization()) {
                C.getVisualizer().end();
            }
            endTime = System.nanoTime();

            // Print final statistics
            System.out.println("\n* Global statistics *");
            System.out.println("Total time: " + (float)(endTime - startTime)/(float)NANOSECONDS_IN_ONE_SECOND + " s");
            System.out.println("Evaluated particles: " + globalStatistics.getEvaluatedParticles());
            System.out.println("Merged particles: " + globalStatistics.getMergedParticles());

            System.out.println("\n* Thread statistics *");
            for (int i = 0; i < threadsStatistics.size(); i++) {
                System.out.println("Thread " + i + ":");
                System.out.print("Total time: " + threadsStatistics.get(i).getTotalTime() /(float)NANOSECONDS_IN_ONE_SECOND + "s, ");
                System.out.print("Evaluated particles: " + threadsStatistics.get(i).getEvaluatedParticles());
                System.out.println(", Merged particles: " + threadsStatistics.get(i).getMergedParticles());
            }

            // Join threads
            for (Thread thread : threads) {
                thread.join();
            }

        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        } finally {
            C.setRunning(false);
            if (properties.isSaveToFile()) {
                C.getReaderWriter().endFile();
            }
        }

        info(logger, "End of simulation. Time: " + nanoToHumanReadable(endTime - startTime));
    }

    public void init(boolean printInfo) throws SimulationException {
        if (printInfo) {
            logger.info("Initialize simulation...");
        }

        switch (properties.getInteractingLaw()) {
            case NEWTON_LAW_OF_GRAVITATION:
                forceCalculator = new NewtonGravityAP();
                break;
            case COULOMB_LAW_ELECTRICALLY:
                throw new NotImplementedException("COULOMB_LAW_ELECTRICALLY is not implemented");
                // break;
            default:
                forceCalculator = new NewtonGravityAP();
                break;
        }

        objects = new ArrayList<>();

        for (SimulationObject simulationObject : properties.getInitialObjects()) {
            objects.add(new SimulationObjectImpl(simulationObject));
        }

        if (duplicateIdExists(objects)) {
            throw new SimulationException("Objects with duplicate IDs exist!");
        }

        if (collisionExists(objects)) {
            throw new SimulationException("Initial collision exists!");
        }
        if (printInfo) {
            logger.info("Done.\n");
            Utils.printConfiguration(this);
        }
    }

    public void initSwitchingFromGpu(List<SimulationObject> currentObjects) {
        logger.info("Switching to CPU - Initialize simulation...");

        switch (properties.getInteractingLaw()) {
            case NEWTON_LAW_OF_GRAVITATION:
                forceCalculator = new NewtonGravityAP();
                break;
            case COULOMB_LAW_ELECTRICALLY:
                throw new NotImplementedException("COULOMB_LAW_ELECTRICALLY is not implemented");
                // break;
            default:
                forceCalculator = new NewtonGravityAP();
                break;
        }

        objects = currentObjects;
        logger.info("Done.\n");

        Utils.printConfiguration(this);
    }

    @Override
    public List<SimulationObject> getObjects() {
        return objects;
    }

    @Override
    public List<SimulationObject> getAuxiliaryObjects() {
        return auxiliaryObjects;
    }

    @Override
    public long getCurrentIterationNumber() {
        return iterationCounter;
    }

    @Override
    public ForceCalculator getForceCalculator() {
        return forceCalculator;
    }

    public SimulationLogicAP getSimulationLogic() {
        return simulationLogic;
    }

    @Override
    public SimulationProperties getProperties() {
        return properties;
    }

    @Override
    public void setProperties(SimulationProperties properties) {
        this.properties = properties;
    }

    @Override
    public Number calculateDistance(ImmutableSimulationObject object, ImmutableSimulationObject object1) {
        return simulationLogic.calculateDistance(object, object1);
    }

    public boolean isCollisionExists() {
        return collisionExists;
    }
    
    public void upCollisionExists() {
        collisionExists = true;
    }
}
