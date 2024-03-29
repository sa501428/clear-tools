package cli.utils.general;

import org.apache.commons.math3.stat.regression.SimpleRegression;

public class ManhattanDecay {

    public static float[] calculateDecay(float[][] matrix, int midX, int midY, int window) {
        float[] values = new float[window + 1];
        int[] counts = new int[window + 1];

        int startR = Math.max(midX - window, 0);
        int endR = Math.min(midX + window + 1, matrix.length);
        int startC = Math.max(midY - window, 0);
        int endC = Math.min(midY + window + 1, matrix[0].length);

        for (int i = startR; i < endR; i++) {
            int di = Math.abs(i - midX);
            for (int j = startC; j < endC; j++) {
                int dj = Math.abs(j - midY);

                int dist = di + dj;
                if (dist < window) {
                    values[dist] += matrix[i][j];
                    counts[dist]++;
                }
            }
        }

        for (int k = 0; k < values.length; k++) {
            if (counts[k] > 0) {
                values[k] /= counts[k];
            }
        }

        /*float denominator = 0 + values[0];
        if (denominator > 0) {
            for (int k = 0; k < values.length; k++) {
                values[k] /= denominator;
            }
        }*/

        return values;
    }

    public static boolean passesMonotonicDecreasing(float[] decay, int maxIndex) {
        for (int i = 0; i < maxIndex; i++) {
            if (decay[i + 1] > decay[i]) {
                return false;
            }
        }
        return true;
    }


    private static String toString(float[] array) {
        StringBuilder str = new StringBuilder("" + array[0]);
        for (int k = 1; k < array.length; k++) {
            str.append(",").append(array[k]);
        }
        return str.toString();
    }

    private static float getDecaySlope(float[] decay, int startIndex) {
        SimpleRegression regression = new SimpleRegression();
        for (int i = startIndex; i < decay.length; i++) {
            regression.addData(Math.log(1 + i), Math.log(1 + decay[i]));
        }
        return (float) regression.getSlope();
    }
}
