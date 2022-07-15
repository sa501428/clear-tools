package cli.utils.expected;

import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;

public class LogExpectedModel {

    private final double[] compressedLogExpected;
    private final double[] compressedExpected;

    public LogExpectedModel(MatrixZoomData zd, NormalizationType norm, int maxBinDist) {
        WelfordStats stats = getSummaryStats(zd, maxBinDist, false, 0, norm);
        compressedLogExpected = stats.getMean();
        compressedExpected = expm1(compressedLogExpected);
    }

    public static WelfordStats getSummaryStats(MatrixZoomData zd, int length, boolean useNone, int minVal,
                                               NormalizationType norm) {
        WelfordStats stats = new WelfordStats(length);

        Iterator<ContactRecord> it;
        if (useNone) {
            it = zd.getDirectIterator();
        } else {
            it = zd.getNormalizedIterator(norm);
        }

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > minVal) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < length) {
                    stats.addValue(dist, logp1(cr.getCounts()));
                }
            }
        }
        return stats;
    }

    public static int logp1i(int x) {
        return (int) Math.floor(Math.log(1 + x));
    }

    public static double logp1(double x) {
        return Math.log(1 + x);
    }

    public double getExpFromBin(int dist) {
        return compressedExpected[logp1i(dist)];
    }

    public double getLogExpFromBin(int dist) {
        return compressedLogExpected[logp1i(dist)];
    }

    private double[] expm1(double[] input) {
        double[] vec = new double[input.length];
        for (int k = 0; k < vec.length; k++) {
            vec[k] = Math.expm1(input[k]);
        }
        return vec;
    }

    public static float getP(double obs, double expected, double superDiagonal) {
        // P = (O - E)/(SD - E)
        return (float) ((obs - expected) / (superDiagonal - expected));
    }

    public float getPercentContact(int dist, float counts) {
        double baseline = getExpFromBin(dist);
        double maxSignal = getExpFromBin(1);
        return getP(counts, baseline, maxSignal);
    }
}
