package cli.utils.general;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class VectorCleaner {
    public static void inPlaceClean(double[] vec) {
        for (int k = 0; k < vec.length; k++) {
            if (Double.isInfinite(vec[k]) || vec[k] < 1e-20) {
                vec[k] = Double.NaN;
            }
        }

        double logThreshold = getLogLowerZscoreBound(vec, -2);
        //double logThreshold = getLogLowerBound(vec);
        double lowerBound = Math.min(1, Math.exp(logThreshold));
        //double lowerBound = Math.min(1, getPercentile(vec, 16));
        for (int k = 0; k < vec.length; k++) {
            if (vec[k] < lowerBound) {
                vec[k] = Double.NaN;
            }
        }
    }

    private static double getLogLowerZscoreBound(double[] vec, float z) {
        Zscore stats = createLogZscore(vec);
        return stats.getZscore(z);
    }

    private static double getPercentile(double[] vec, int percentile) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double v : vec) {
            if (v > 0) {
                stats.addValue(v);
            }
        }
        return stats.getPercentile(percentile);
    }

    private static double getLogLowerBound(double[] vec) {
        DescriptiveStatistics stats = createLogStats(vec);
        double p25 = stats.getPercentile(25);
        double iqr = stats.getPercentile(75) - p25;
        return p25 - (iqr * 1.5);
    }

    private static DescriptiveStatistics createLogStats(double[] vec) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double v : vec) {
            if (v > 0) {
                stats.addValue(Math.log(v));
            }
        }
        return stats;
    }

    private static Zscore createLogZscore(double[] vec) {
        Welford welford = new Welford();
        for (double v : vec) {
            if (v > 0) {
                welford.addValue(Math.log(v));
            }
        }
        return welford.getZscore();
    }
}
