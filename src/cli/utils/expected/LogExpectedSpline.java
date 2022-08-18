package cli.utils.expected;

import javastraw.reader.block.ContactRecord;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class LogExpectedSpline {

    private final int n;
    private final PolynomialSplineFunction mu;

    public LogExpectedSpline(double[] mean, double[] stddev) {
        n = mean.length;
        double[] indices = new double[n];
        for (int i = 0; i < n; i++) {
            indices[i] = i;
        }

        SplineInterpolator interpolator = new SplineInterpolator();
        mu = interpolator.interpolate(indices, mean);
    }

    public static float getP(double obs, double expected, double superDiagonal) {
        // P = (O - E)/(SD - E)
        return (float) ((obs - expected) / (superDiagonal - expected));
    }

    public boolean isReasonablePercentContact(ContactRecord cr, LogExpectedModel model) {
        double percentContact = getPercentContact(cr);
        return percentContact > 0.01 && percentContact < 0.4;
    }

    public float getPercentContact(ContactRecord cr) {
        double dist = Math.log(1 + ExpectedUtils.getDist(cr));
        dist = Math.min(dist, n - 1);
        double baseline = mu.value(dist);
        double maxSignal = mu.value(1);
        return getP(cr.getCounts(), baseline, maxSignal);
    }

    public void print() {
        for (double k = 0; k <= n - 1; k += 0.5) {
            System.out.println(k + "  " + mu.value(k));
        }
    }
}
