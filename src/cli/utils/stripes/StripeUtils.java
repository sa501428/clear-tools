package cli.utils.stripes;

import cli.utils.anchors.Convolution1DTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TTest;

import java.util.BitSet;
import java.util.LinkedList;
import java.util.List;

public class StripeUtils {

    public final static int MIDPT = 2, PRE = 1, POST = 3;
    public final static float RELATIVE_ENRICHMENT = 1.25f;
    public static float GEO_ENRICHMENT = 1f;

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

        return stretches;
    }

    private static float[][] smoothen(float[][] dataSlice) {
        float[][] smooth = new float[dataSlice.length][dataSlice[0].length];
        for (int i = 0; i < dataSlice.length; i++) {
            smooth[i] = Convolution1DTools.smooth5(dataSlice[i]);
        }
        return smooth;
    }

    public static List<int[]> onlySignificantStretches(List<int[]> stretches, float[][] dataSlice, float[][] dataOESlice) {
        List<int[]> result = new LinkedList<>();
        for (int[] stretch : stretches) {

            if (hasGeoMeanHigherThanNeighbors(stretch, dataOESlice)) {
                if (isMeanSignificantlyHigherThanNeighbors(stretch, dataSlice, 0)) {
                    if (isMeanSignificantlyHigherThanNeighbors(stretch, dataOESlice, 1)) {
                        result.add(stretch);
                    }
                }
            }
        }
        return result;
    }

    public static boolean hasGeoMeanHigherThanNeighbors(int[] stretch, float[][] dataOESlice) {
        double geoMean = getGeometricMean(stretch, dataOESlice[MIDPT], 1);
        double prevMean = getGeometricMean(stretch, dataOESlice[PRE], 1);
        double nextMean = getGeometricMean(stretch, dataOESlice[POST], 1);
        // RELATIVE_ENRICHMENT *
        return geoMean > 1.5 && geoMean > GEO_ENRICHMENT * prevMean && geoMean > GEO_ENRICHMENT * nextMean;
    }

    public static boolean isMeanSignificantlyHigherThanNeighbors(int[] stretch, float[][] data, int minVal) {
        double[] potentialStripe = getValuesInBoundsOfRows(stretch, data[MIDPT], minVal);

        List<Double> prev = new LinkedList<>();
        addValuesInBoundsOfRows(stretch, data[PRE], minVal, prev);
        addValuesInBoundsOfRows(stretch, data[PRE - 1], minVal, prev);
        double[] prevStripe = toArray(prev);

        List<Double> next = new LinkedList<>();
        addValuesInBoundsOfRows(stretch, data[POST], minVal, next);
        addValuesInBoundsOfRows(stretch, data[POST + 1], minVal, next);
        double[] nextStripe = toArray(next);


        double[] background = getBackgroundValues(stretch, data, minVal);

        boolean status = isMeanSignificantlyHigher(potentialStripe, background);
        if (prevStripe.length > 5) {
            status = status && isMeanSignificantlyHigher(potentialStripe, prevStripe);
        }
        if (nextStripe.length > 5) {
            status = status && isMeanSignificantlyHigher(potentialStripe, nextStripe);
        }
        return status;
    }

    private static double[] getBackgroundValues(int[] stretch, float[][] dataSlice, int minVal) {
        List<Double> values = new LinkedList<>();
        for (int i = 0; i < dataSlice.length; i++) {
            if (i != MIDPT) {
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

    private static void addValuesInBoundsOfRows(int[] stretch, float[] dataSlice, int minVal, List<Double> values) {
        for (int j = stretch[0]; j < stretch[1]; j++) {
            if (dataSlice[j] > minVal) {
                values.add((double) dataSlice[j]);
            }
        }
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
            if (k != MIDPT && dataSlice[MIDPT][j] > RELATIVE_ENRICHMENT * dataSlice[k][j]) {
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
