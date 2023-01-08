package cli.utils.peaks;

import org.apache.commons.math3.ml.clustering.Clusterable;

public class Point1D implements Clusterable {
    private final long x;

    public Point1D(long x) {
        this.x = x;
    }

    public long getX() {
        return x;
    }

    @Override
    public double[] getPoint() {
        return new double[]{x};
    }
}
