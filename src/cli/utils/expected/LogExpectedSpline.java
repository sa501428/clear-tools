package cli.utils.expected;

import javastraw.reader.block.ContactRecord;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class LogExpectedSpline extends ExpectedModel {

    private final int n;
    private final PolynomialSplineFunction mu;
    double maxSignal;

    public LogExpectedSpline(double[] mean, double[] stddev) {
        n = mean.length;
        double[] indices = new double[n];
        for (int i = 0; i < n; i++) {
            indices[i] = i;
        }

        SplineInterpolator interpolator = new SplineInterpolator();
        mu = interpolator.interpolate(indices, mean);
        maxSignal = Math.expm1(mu.value(0.5));
    }

    public boolean isReasonablePercentContact(ContactRecord cr, LogExpectedModel model) {
        double percentContact = getPercentContact(cr);
        return percentContact > 0.01;// && percentContact < 0.4;
    }

    public static float getP(double obs, double expected, double superDiagonal) {
        // P = (O - E)/(SD - E)
        return (float) ((obs - expected) / (superDiagonal - expected));
    }

    public float getPercentContact(ContactRecord cr) {
        double dist = (Math.log(1 + ExpectedUtils.getDist(cr))); // floor
        dist = Math.min(dist, n - 1);
        double baseline = Math.expm1(mu.value(dist));
        return getP(cr.getCounts(), baseline, maxSignal);
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
}
