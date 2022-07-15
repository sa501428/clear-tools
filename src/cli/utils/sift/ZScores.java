package cli.utils.sift;

public class ZScores {

    private final double[] mean;
    private final double[] stdDev;

    public ZScores(double[] mean, double[] stdDev) {
        this.mean = mean;
        this.stdDev = stdDev;
    }

    public float getZscore(int index, double value) {
        return (float) ((value - mean[index]) / stdDev[index]);
    }
}
