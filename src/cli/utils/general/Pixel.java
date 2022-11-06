package cli.utils.general;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.LinkedList;
import java.util.List;

public class Pixel {
    public final int row;
    public final int col;
    public final float value;
    public final float zScore;

    public Pixel(int row, int col, float value, float zScore) {
        this.row = row;
        this.col = col;
        this.value = value;
        this.zScore = zScore;
    }

    public static boolean contains(Pixel px, int minR, int maxR, int minC, int maxC) {
        return Utils.inBounds(px.row, minR, maxR) && Utils.inBounds(px.col, minC, maxC);
    }

    public static Pixel getMax(List<Pixel> pixels) {
        Pixel max = pixels.get(0);
        for (Pixel px : pixels) {
            if (px.value > max.value) {
                max = px;
            }
        }
        return max;
    }

    public static Pixel getMax(List<Pixel> pixels, float[] rowSignal, float[] colSignal, int peakWidthLimit) {
        int i = getMaxIndex(rowSignal, peakWidthLimit);
        int j = getMaxIndex(colSignal, peakWidthLimit);

        if (i > -1 && j > -1) {
            for (Pixel pixel : pixels) {
                if (pixel.is(i, j)) return pixel;
            }
        }
        return null;
    }

    private static int getMaxIndex(float[] signal, int limit) {
        double weightedSum = 0;
        double weights = 0;
        float cutoff = 3 * getMedian(signal);
        for (int i = 0; i < signal.length; i++) {
            if (signal[i] > 0.5 && signal[i] > cutoff) {
                weightedSum += signal[i] * i;
                weights += signal[i];
            }
        }

        if (weights > 0) {
            return (int) (weightedSum / weights);
        }

        return -1;
    }

    private static int getMaxIndexV1(float[] signal, int limit) {
        List<Integer> goodIndices = new LinkedList<>();
        float cutoff = 3 * getMedian(signal);
        for (int i = 0; i < signal.length; i++) {
            if (signal[i] > 0.5 && signal[i] > cutoff) {
                goodIndices.add(i);
            }
        }

        if (goodIndices.size() > 0) {
            return mean(goodIndices, limit);
        }
        return -1;
    }

    private static int mean(List<Integer> goodIndices, int limit) {
        if (goodIndices.size() == 1) return goodIndices.get(0);
        long total = 0;
        int minVal = goodIndices.get(0);
        int maxVal = goodIndices.get(0);
        for (int val : goodIndices) {
            total += val;
            minVal = Math.min(val, minVal);
            maxVal = Math.max(val, maxVal);
        }
        //if(maxVal - minVal < limit) {
        return (int) (total / goodIndices.size());
        //}
        //return -1;
    }

    private static int getMaxIndexV0(float[] signal) {
        int m = 0;
        for (int i = 0; i < signal.length; i++) {
            if (signal[m] < signal[i]) {
                m = i;
            }
        }

        float median = getMedian(signal);

        if (signal[m] / median > 3) {
            return m;
        }
        return -1;
    }

    private static float getMedian(float[] signal) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (float val : signal) {
            if (val > 0) {
                stats.addValue(val);
            }
        }
        return (float) stats.getPercentile(50);
    }

    private boolean is(int i, int j) {
        return i == row && j == col;
    }
}
