package cli.utils.stripes;

import cli.utils.anchors.Convolution1DTools;
import org.apache.commons.math3.stat.inference.TTest;

import java.util.BitSet;
import java.util.LinkedList;
import java.util.List;

public class StripeUtils {
    public static List<int[]> findContiguousStretches(float[][] dataSlice, int minLength) {
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
                    stretches.add(new int[]{start, end});
                }
                // Reset the counters for the next potential stretch
                start = -1;
                end = -1;
                currentLength = 0;
            }
        }

        // Check if the last stretch is longer than length 10
        if (currentLength > minLength) {
            stretches.add(new int[]{start, end});
        }

        return onlySignificantStretches(stretches, dataSlice);
    }

    private static float[][] smoothen(float[][] dataSlice) {
        float[][] smooth = new float[dataSlice.length][dataSlice[0].length];
        for (int i = 0; i < dataSlice.length; i++) {
            smooth[i] = Convolution1DTools.smooth5(dataSlice[i]);
        }
        return smooth;
    }

    private static List<int[]> onlySignificantStretches(List<int[]> stretches, float[][] dataSlice) {
        List<int[]> result = new LinkedList<>();
        for (int[] stretch : stretches) {
            double[] potentialStripe = getStripeValues(stretch, dataSlice);
            double[] background = getBackgroundValues(stretch, dataSlice);

            if (isMeanSignificantlyHigher(potentialStripe, background)) {
                result.add(stretch);
            }
        }
        return result;
    }

    private static double[] getBackgroundValues(int[] stretch, float[][] dataSlice) {
        double[][] result = new double[4][stretch[1] - stretch[0] + 1];
        for (int i = 0; i < result.length; i++) {
            result[0][i] = dataSlice[0][i + stretch[0]];
            result[1][i] = dataSlice[1][i + stretch[0]];
            result[2][i] = dataSlice[3][i + stretch[0]];
            result[3][i] = dataSlice[4][i + stretch[0]];
        }
        return flatten(result);
    }

    private static double[] flatten(double[][] result) {
        double[] flat = new double[result.length * result[0].length];
        for (int i = 0; i < result.length; i++) {
            System.arraycopy(result[i], 0, flat, i * result[0].length, result[0].length);
        }
        return flat;
    }

    private static double[] getStripeValues(int[] stretch, float[][] dataSlice) {
        double[] result = new double[stretch[1] - stretch[0] + 1];
        for (int i = 0; i < result.length; i++) {
            result[i] = dataSlice[2][i + stretch[0]];
        }
        return result;
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
        TTest tTest = new TTest();
        double pValue = tTest.tTest(listA, listB) / 2;

        // Check if the p-value is less than the significance level
        return pValue < 0.01;
    }
}
