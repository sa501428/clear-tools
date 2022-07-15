package cli.utils;

import cli.utils.sift.Zscore;

public class Welford {
    private long counts = 0;
    private double mu = 0;
    private double aggSquaredDiffs = 0;

    public void addValue(double x) {
        counts++;
        double nextMu = mu + ((x - mu) / counts);
        aggSquaredDiffs += (x - mu) * (x - nextMu);
        mu = nextMu;
    }

    public double getMean() {
        return mu;
    }

    public double getStdDev() {
        if (counts > 2) {
            return Math.sqrt(aggSquaredDiffs / (counts - 1));
        }
        return 0;
    }

    public Zscore getZscore() {
        return new Zscore(getMean(), getStdDev());
    }

    public long getCounts() {
        return counts;
    }
}



