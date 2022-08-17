package cli.utils.expected;

import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import cli.utils.sift.ExtremePixels;
import cli.utils.sift.ZScores;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;
import java.util.List;

public class LogExpectedModel {

    private final WelfordStats stats;
    private final double[] compressedExpected;

    public LogExpectedModel(MatrixZoomData zd, NormalizationType norm, int maxBinDist, int minVal) {
        stats = getSummaryStats(zd, maxBinDist, minVal, norm);
        compressedExpected = expm1(stats.getMean());
    }

    public LogExpectedModel(List<ContactRecord> records, int maxBinDist) {
        stats = getSummaryStats(records, maxBinDist);
        compressedExpected = expm1(stats.getMean());
    }

    private WelfordStats getSummaryStats(MatrixZoomData zd, int maxBin, int minVal,
                                         NormalizationType norm) {
        Iterator<ContactRecord> it = ExtremePixels.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > minVal) {
                updateStats(cr, maxBin, stats);
            }
        }
        return stats;
    }

    private WelfordStats getSummaryStats(List<ContactRecord> records, int maxBin) {
        WelfordStats stats = new WelfordStats(logp1i(maxBin) + 1);
        for (ContactRecord cr : records) {
            updateStats(cr, maxBin, stats);
        }
        return stats;
    }

    private void updateStats(ContactRecord cr, int maxBin, WelfordStats stats) {
        int dist = ExpectedUtils.getDist(cr);
        if (dist < maxBin) {
            stats.addValue(logp1i(dist), logp1(cr.getCounts()));
        }
    }

    public int logp1i(int x) {
        return (int) Math.floor(Math.log(1 + x));
    }

    public static double logp1(double x) {
        return Math.log(1 + x);
    }

    public double getExpFromUncompressedBin(int dist) {
        return compressedExpected[logp1i(dist)];
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
        double baseline = getExpFromUncompressedBin(dist);
        double maxSignal = getExpFromUncompressedBin(1);
        return getP(counts, baseline, maxSignal);
    }

    public ZScores getZscores() {
        return stats.getZscores();
    }
}
