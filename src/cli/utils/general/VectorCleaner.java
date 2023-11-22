package cli.utils.general;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class VectorCleaner {

    public static void inPlaceZscore(double[] vec) {
        cleanExtremes(vec);
        Zscore stats = createLogZscore(vec);
        for (int k = 0; k < vec.length; k++) {
            vec[k] = stats.getZscore(Math.log(vec[k]));
        }
    }

    public static void inPlaceClean(double[] vec) {
        cleanExtremes(vec);
        double logThreshold = getLogLowerZscoreBound(vec, -2);
        double lowerBound = Math.min(1, Math.exp(logThreshold));
        for (int k = 0; k < vec.length; k++) {
            if (vec[k] < lowerBound) {
                vec[k] = Double.NaN;
            }
        }
    }

    private static void cleanExtremes(double[] vec) {
        for (int k = 0; k < vec.length; k++) {
            if (Double.isInfinite(vec[k]) || vec[k] < 1e-20) {
                vec[k] = Double.NaN;
            }
        }
    }

    private static double getLogLowerZscoreBound(double[] vec, float z) {
        Zscore stats = createLogZscore(vec);
        return stats.getValForZscore(z);
    }

    public static double getPercentile(double[] vec, int percentile, int minVal) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double v : vec) {
            if (v > minVal) {
                stats.addValue(v);
            }
        }
        return stats.getPercentile(percentile);
    }

    public static int getPercentile(int[] vec, int percentile, int minVal) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int v : vec) {
            if (v > minVal) {
                stats.addValue(v);
            }
        }
        return (int) stats.getPercentile(percentile);
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

    public static float[] logMedianScale(float[] original) {
        float[] data = new float[original.length];
        System.arraycopy(original, 0, data, 0, original.length);
        inPlaceLog(data);
        cleanNanInfTo(data, 0);
        double median = getMedian(data, 0);
        for (int k = 0; k < data.length; k++) {
            data[k] = (float) (data[k] / median);
        }
        return data;
    }

    private static double getMedian(float[] data, int minVal) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (float v : data) {
            if (v > minVal) {
                stats.addValue(v);
            }
        }
        return stats.getPercentile(50);
    }

    private static void inPlaceLog(float[] data) {
        for (int k = 0; k < data.length; k++) {
            data[k] = (float) Math.log(data[k]);
        }
    }

    private static void cleanNanInfTo(float[] data, int value) {
        for (int k = 0; k < data.length; k++) {
            if (Float.isNaN(data[k]) || Float.isInfinite(data[k])) {
                data[k] = value;
            }
        }
    }
}
