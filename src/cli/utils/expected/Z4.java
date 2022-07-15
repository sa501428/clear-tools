package cli.utils.expected;

import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import cli.utils.sift.Z4Scores;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.Iterator;

public class Z4 {
    private static Z4Scores getZ4scores(MatrixZoomData zd, int length,
                                        double[] nvVC, double[] nvVCSqrt, double[] nvSCALE) {
        WelfordStats rawStats = new WelfordStats(length);
        WelfordStats vcStats = new WelfordStats(length);
        WelfordStats vcSqrtStats = new WelfordStats(length);
        WelfordStats scaleStats = new WelfordStats(length);
        Iterator<ContactRecord> it = zd.getDirectIterator();

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = LogExpectedModel.logp1i(ExpectedUtils.getDist(cr));
                if (dist < length) {
                    rawStats.addValue(dist, LogExpectedModel.logp1(cr.getCounts()));
                    populateNormedStats(cr, nvVC, vcStats, dist);
                    populateNormedStats(cr, nvVCSqrt, vcSqrtStats, dist);
                    populateNormedStats(cr, nvSCALE, scaleStats, dist);
                }
            }
        }
        return new Z4Scores(rawStats, vcStats, vcSqrtStats, scaleStats);
    }

    private static void populateNormedStats(ContactRecord cr, double[] norm, WelfordStats stats, int dist) {
        double denom = norm[cr.getBinX()] * norm[cr.getBinY()];
        if (denom > 0) {
            double val = cr.getCounts() / denom;
            if (val > 1) {
                stats.addValue(dist, LogExpectedModel.logp1(val));
            }
        }
    }
}
