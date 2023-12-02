package cli.utils.stripes;

import cli.utils.anchors.Convolution1DTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TTest;

import java.util.BitSet;
import java.util.LinkedList;
import java.util.List;

public class StripeUtils {
    public static List<int[]> findContiguousStretches(float[][] dataSlice, float[][] dataOESlice, int minLength) {
        float[][] smooth = smoothen(dataSlice);
        BitSet enriched = getEnriched(smooth);
        List<int[]> stretches = new LinkedList<>();

        int start = -1;
        int end = -1;
        int currentLength = 0;

        for (int i = 0; i < enriched.size(); i++) {
            if (enriched.get(i)) {
                if (start == -1) {
                    start = i;
                }
                currentLength++;
                end = i;
            } else {
                if (currentLength > minLength) {
                    stretches.add(new int[]{start, end + 1}); // non-inclusive second
                }
                // Reset the counters for the next potential stretch
                start = -1;
                end = -1;
                currentLength = 0;
            }
        }

        // Check if the last stretch is longer than length 10
        if (currentLength > minLength) {
            stretches.add(new int[]{start, end + 1});
        }

        return onlySignificantStretches(stretches, dataSlice, dataOESlice);
    }

    private static float[][] smoothen(float[][] dataSlice) {
        float[][] smooth = new float[dataSlice.length][dataSlice[0].length];
        for (int i = 0; i < dataSlice.length; i++) {
            smooth[i] = Convolution1DTools.smooth5(dataSlice[i]);
        }
        return smooth;
    }

    private static List<int[]> onlySignificantStretches(List<int[]> stretches, float[][] dataSlice, float[][] dataOESlice) {
        List<int[]> result = new LinkedList<>();
        for (int[] stretch : stretches) {

            double geoMean = getGeometricMean(stretch, dataSlice[2], 0);
            double prevMean = getGeometricMean(stretch, dataSlice[1], 1);
            double nextMean = getGeometricMean(stretch, dataSlice[3], 1);

            if (geoMean > 2 && geoMean > 1.5 * prevMean && geoMean > 1.5 * nextMean) {
                double[] potentialStripe = getValuesInBoundsOfRows(stretch, dataSlice[2], 0);
                double[] prevStripe = getValuesInBoundsOfRows(stretch, dataSlice[1], 0);
                double[] nextStripe = getValuesInBoundsOfRows(stretch, dataSlice[3], 0);
                double[] background = getBackgroundValues(stretch, dataSlice, 0);

                if (isMeanSignificantlyHigher(potentialStripe, background)
                        && isMeanSignificantlyHigher(potentialStripe, prevStripe)
                        && isMeanSignificantlyHigher(potentialStripe, nextStripe)
                ) {
                    potentialStripe = getValuesInBoundsOfRows(stretch, dataOESlice[2], 1);
                    prevStripe = getValuesInBoundsOfRows(stretch, dataOESlice[1], 1);
                    nextStripe = getValuesInBoundsOfRows(stretch, dataOESlice[3], 1);
                    background = getBackgroundValues(stretch, dataOESlice, 1);

                    if (isMeanSignificantlyHigher(potentialStripe, background)
                            && isMeanSignificantlyHigher(potentialStripe, prevStripe)
                            && isMeanSignificantlyHigher(potentialStripe, nextStripe)
                    ) {
                        result.add(stretch);
                    }
                }
            }
        }
        return result;
    }

    private static double[] getBackgroundValues(int[] stretch, float[][] dataSlice, int minVal) {
        List<Double> values = new LinkedList<>();
        for (int i = 0; i < dataSlice.length; i++) {
            if (i != 2) {
                for (int j = stretch[0]; j < stretch[1]; j++) {
                    if (dataSlice[i][j] > minVal) {
                        values.add((double) dataSlice[i][j]);
                    }
                }
            }
        }

        return toArray(values);
    }

    private static double[] toArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = values.get(i);
        }
        return result;
    }

    private static double[] getValuesInBoundsOfRows(int[] stretch, float[] dataSlice, int minVal) {
        List<Double> values = new LinkedList<>();
        for (int j = stretch[0]; j < stretch[1]; j++) {
            if (dataSlice[j] > minVal) {
                values.add((double) dataSlice[j]);
            }
        }
        return toArray(values);
    }

    private static double getGeometricMean(int[] stretch, float[] dataSlice, int minVal) {
        double result = 1;
        int count = 0;
        for (int i = 0; i < stretch[1] - stretch[0]; i++) {
            if (dataSlice[i + stretch[0]] > minVal) {
                result *= dataSlice[i + stretch[0]];
                count++;
            }
        }
        return Math.pow(result, 1.0 / count);
    }

    private static BitSet getEnriched(float[][] dataSlice) {
        BitSet enriched = new BitSet(dataSlice[0].length);
        for (int j = 0; j < dataSlice[0].length; j++) {
            enriched.set(j, enrichedComparedToRows(dataSlice, j));
        }
        fixSurroundedFalse(enriched);
        return enriched;
    }

    private static void fixSurroundedFalse(BitSet input) {
        for (int i = 1; i < input.size() - 1; i++) {
            if (!input.get(i) && input.get(i - 1) && input.get(i + 1)) {
                input.set(i);
            }
        }
    }

    private static boolean enrichedComparedToRows(float[][] dataSlice, int j) {
        int count = 0;
        for (int k = 0; k < 5; k++) {
            if (k != 2 && dataSlice[2][j] > 1.25 * dataSlice[k][j]) {
                count++;
            }
        }
        return count > 2;
    }

    public static boolean isMeanSignificantlyHigher(double[] listA, double[] listB) {
        // Perform 1 sided t-test
        if (listA.length > 5 && listB.length > 5) {
            TTest tTest = new TTest();
            double pValue = tTest.tTest(listA, listB) / 2;

            // Check if the p-value is less than the significance level
            return pValue < 0.01; // && perc(listA,25) > 2*perc(listB,75)
        }
        return false;
    }

    private static int perc(double[] numbers, int n) {
        return (int) new DescriptiveStatistics(numbers).getPercentile(n);
    }

    private static double[] removeOutliers(double[] data) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double value : data) {
            stats.addValue(value);
        }

        double q1 = stats.getPercentile(25);
        double q3 = stats.getPercentile(75);

        List<Double> filteredData = new LinkedList<>();
        for (double value : data) {
            if (value >= q1 && value <= q3) {
                filteredData.add(value);
            }
        }

        double[] result = new double[filteredData.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = filteredData.get(i);
        }
        return result;
    }
}
