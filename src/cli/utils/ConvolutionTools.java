package cli.utils;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class ConvolutionTools {

    public static float[][] sparseConvolution(int[][] image) {

        float[][] kernel = getManhattanKernel(getApproxBandWidth(image));
        float maxK = ArrayTools.getMax(kernel);

        float[][] result = new float[image.length][image[0].length];

        int halfWidth = kernel.length / 2;
        int maxR = image.length;
        int maxC = image[0].length;

        for (int i = 0; i < maxR; i++) {
            for (int j = 0; j < maxC; j++) {
                if (image[i][j] > 0) {
                    for (int ki = 0; ki < kernel.length; ki++) {
                        int r = i - halfWidth + ki;
                        if (r >= 0 && r < maxR) {
                            for (int kj = 0; kj < kernel[ki].length; kj++) {
                                int c = j - halfWidth + kj;
                                if (c >= 0 && c < maxC) {
                                    result[r][c] += image[i][j] * kernel[ki][kj];
                                }
                            }
                        }
                    }
                }
            }
        }

        ArrayTools.inPlaceDivideArrayBy(result, maxK);

        return result;
    }

    private static float[][] getManhattanKernel(int width) {
        float[][] kernel = new float[width][width];
        int halfWidth = width / 2;
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                int d2 = Math.abs(i - halfWidth) + Math.abs(j - halfWidth);
                if (d2 < halfWidth) {
                    kernel[i][j] = halfWidth - d2;
                }
            }
        }

        // normalize sum to 1
        float sum = ArrayTools.getSum(kernel);
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < width; j++) {
                kernel[i][j] /= sum;
            }
        }
        return kernel;
    }

    private static int getApproxBandWidth(int[][] image) {
        int bw = Math.max(getApproxRowBandWidth(image), getApproxColBandWidth(image));
        return 4 * bw + 1; // empirically sized
    }

    private static int getApproxColBandWidth(int[][] image) {
        int[] colHistogram = new int[image[0].length];
        for (int[] ints : image) {
            for (int c = 0; c < ints.length; c++) {
                colHistogram[c] += ints[c];
            }
        }
        return getApproxBandWidthFromHist(colHistogram);
    }

    private static int getApproxRowBandWidth(int[][] image) {
        int[] rowHistogram = new int[image.length];
        for (int r = 0; r < image.length; r++) {
            for (int c = 0; c < image[r].length; c++) {
                rowHistogram[r] += image[r][c];
            }
        }
        return getApproxBandWidthFromHist(rowHistogram);
    }

    private static int getApproxBandWidthFromHist(int[] histogram) {
        int n = ArrayTools.getSum(histogram);
        DescriptiveStatistics stats = createStatsFromHist(histogram);
        double iqr = stats.getPercentile(75) - stats.getPercentile(25);
        double sigma = stats.getStandardDeviation();
        double h = 0.9 * Math.min(sigma, iqr / 1.34) / Math.pow(n, 0.2);
        return (int) (h) + 1;
    }

    private static DescriptiveStatistics createStatsFromHist(int[] histogram) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < histogram.length; i++) {
            for (int z = 0; z < histogram[i]; z++) {
                stats.addValue(i);
            }
        }
        return stats;
    }
}
