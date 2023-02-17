/*
Práctica 4.
Código fuente: Statistics.java
Grau Informàtica
Pere Muñoz Figuerol
*/

package info.trekto.jos.core.impl.arbitrary_precision;

public class Statistics {
    private float totalTime;
    private float timePerMIterations;
    private int evaluatedParticles;
    private int mergedParticles;

    public Statistics() {
        totalTime = 0;
        timePerMIterations = 0;
        evaluatedParticles = 0;
        mergedParticles = 0;
    }

    public void addTime(float time) {
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

    public float getTotalTime() {
        return totalTime;
    }

    public float getTimePerMIterations() {
        return timePerMIterations;
    }

    public int getEvaluatedParticles() {
        return evaluatedParticles;
    }

    public int getMergedParticles() {
        return mergedParticles;
    }
}
