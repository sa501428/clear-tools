package cli.utils.expected;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class LogExpectedSpline extends ExpectedModel {

    private final int n;
    private final PolynomialSplineFunction mu;
    private final double nearDiagonalSignal;

    public LogExpectedSpline(double[] compressedMu) {
        n = compressedMu.length;
        double[] indices = new double[n];
        for (int i = 0; i < n; i++) {
            indices[i] = i;
        }

        SplineInterpolator interpolator = new SplineInterpolator();
        mu = interpolator.interpolate(indices, compressedMu);
        nearDiagonalSignal = Math.expm1(mu.value(0.5));
    }

    public void print() {
        for (double k = 0; k <= n - 1; k += 0.5) {
            System.out.println(k + "  " + mu.value(k));
        }
    }

    @Override
    public double getExpectedFromUncompressedBin(int dist0) {
        double dist = Math.min(logp1(dist0), n - 1);
        return Math.expm1(mu.value(dist));
    }

    @Override
    public double getNearDiagonalSignal() {
        return nearDiagonalSignal;
    }
}
