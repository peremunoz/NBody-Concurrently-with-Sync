package info.trekto.jos.core.impl.arbitrary_precision;

public class Statistics {
    private long totalTime;
    private long timePerMIterations;
    private int evaluatedParticles;
    private int mergedParticles;

    public Statistics() {
        totalTime = 0;
        timePerMIterations = 0;
        evaluatedParticles = 0;
        mergedParticles = 0;
    }

    public void addTime(long time) {
        totalTime += time;
        timePerMIterations += time;
    }

    public void resetTimePerMIterations() {
        timePerMIterations = 0;
    }

    public void addEvaluatedParticles(int particles) {
        evaluatedParticles += particles;
    }

    public void addMergedParticles(int particles) {
        mergedParticles += particles;
    }

    public long getTotalTime() {
        return totalTime;
    }

    public long getTimePerMIterations() {
        return timePerMIterations;
    }

    public int getEvaluatedParticles() {
        return evaluatedParticles;
    }

    public int getMergedParticles() {
        return mergedParticles;
    }
}
