package cli.utils.expected;

import cli.utils.sift.ExtremePixels;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;
import java.util.List;

public class LogBinnedExpectedModel extends ExpectedModel {

    private final WelfordArray stats;
    private final double[] compressedExpected;

    public LogBinnedExpectedModel(MatrixZoomData zd, NormalizationType norm, int maxBinDist, int minVal) {
        stats = getSummaryStats(zd, maxBinDist, minVal, norm);
        compressedExpected = expm1(stats.getMean());
    }

    public LogBinnedExpectedModel(List<ContactRecord> records, int maxBinDist) {
        stats = getSummaryStats(records, maxBinDist);
        compressedExpected = expm1(stats.getMean());
    }

    private WelfordArray getSummaryStats(MatrixZoomData zd, int maxBin, int minVal,
                                         NormalizationType norm) {
        WelfordArray stats = new WelfordArray(logp1i(maxBin) + 1);
        Iterator<ContactRecord> it = ExtremePixels.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > minVal) {
                updateStats(cr, maxBin, stats);
            }
        }
        return stats;
    }

    private WelfordArray getSummaryStats(List<ContactRecord> records, int maxBin) {
        WelfordArray stats = new WelfordArray(logp1i(maxBin) + 1);
        for (ContactRecord cr : records) {
            updateStats(cr, maxBin, stats);
        }
        return stats;
    }

    private void updateStats(ContactRecord cr, int maxBin, WelfordArray stats) {
        int dist = ExpectedUtils.getDist(cr);
        if (dist < maxBin) {
            stats.addValue(logp1i(dist), logp1(cr.getCounts()));
        }
    }

    public double getExpectedFromUncompressedBin(int dist) {
        return compressedExpected[logp1i(dist)];
    }

    public ZScoreArray getZscores() {
        return stats.getZscores();
    }

    public LogExpectedSpline getSpline() {
        return new LogExpectedSpline(stats.getMean(), stats.getStdDev());
    }
}
