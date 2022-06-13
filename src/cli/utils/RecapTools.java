package cli.utils;

import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.HashMap;
import java.util.Map;

public class RecapTools {

    public static Map<String, String> getStats(Feature2D loop,
                                               float[][] obsMatrix, float[][] eMatrix,
                                               int resolution, int window, ExpectedValueFunction df,
                                               double superDiagonal, float pseudoCount) {
        long binXStart = (loop.getMidPt1() / resolution) - window;
        long binYStart = (loop.getMidPt2() / resolution) - window;

        float[][] oeMatrix = divide(obsMatrix, eMatrix, pseudoCount);

        Map<String, String> loopAttributes = new HashMap<>();
        float obs = obsMatrix[window][window];
        float expected = eMatrix[window][window];

        loopAttributes.put("PRESENCE", String.valueOf(getP(obs, expected, superDiagonal)));
        loopAttributes.put("PRESENCE_INF", String.valueOf(getP(obs, pseudoCount, superDiagonal)));

        addAttributes(loopAttributes, "OBS_", obsMatrix, window);
        addAttributes(loopAttributes, "OE_", oeMatrix, window);

        return loopAttributes;
    }

    private static void addAttributes(Map<String, String> attributes, String stem, float[][] matrix, int window) {
        DescriptiveStatistics stats = makeStats(matrix);

        float val = matrix[window][window];

        attributes.put(stem + "VAL", String.valueOf(val));

        attributes.put(stem + "STD_DEV", String.valueOf(stats.getStandardDeviation()));
        attributes.put(stem + "KURTOSIS", String.valueOf(stats.getKurtosis()));
        attributes.put(stem + "SKEWNESS", String.valueOf(stats.getSkewness()));

        attributes.put(stem + "MEAN_ENRICHMENT", String.valueOf(val / stats.getMean()));
        attributes.put(stem + "MEDIAN_ENRICHMENT", String.valueOf(val / stats.getPercentile(50)));
        attributes.put(stem + "GEO_ENRICHMENT", String.valueOf(val / stats.getGeometricMean()));
        attributes.put(stem + "MAX_ENRICHMENT", String.valueOf(val / stats.getMax()));
        attributes.put(stem + "MIN_ENRICHMENT", String.valueOf(val / stats.getMin()));

        float[] manhattanDecay = calculateDecay(matrix, window);
        addRegressionStats(manhattanDecay, attributes, stem);
    }

    private static void addRegressionStats(float[] decay, Map<String, String> attributes, String stem) {
        SimpleRegression regression = new SimpleRegression();
        for (int i = 0; i < decay.length; i++) {
            regression.addData(i, Math.log(decay[i]));
        }
        // y = Ae^(kx)  for each (xi, yi) --> wi = log(yi)
        // solve w = a + bx
        // e^a = A  ;  b = k
        attributes.put(stem + "DECAY_A", String.valueOf(Math.exp(regression.getIntercept())));
        attributes.put(stem + "DECAY_k", String.valueOf(regression.getSlope()));
    }

    private static float[] calculateDecay(float[][] matrix, int window) {
        float[] values = new float[window + 1];
        int[] counts = new int[window + 1];

        for (int i = 0; i < matrix.length; i++) {
            int di = Math.abs(i - window);
            for (int j = 0; j < matrix[i].length; j++) {
                int dj = Math.abs(j - window);

                int dist = di + dj;
                if (dist < values.length) {
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

        float denom = 0 + values[0];
        for (int k = 0; k < values.length; k++) {
            values[k] /= denom;
        }

        return values;
    }

    private static DescriptiveStatistics makeStats(float[][] matrix) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (float[] row : matrix) {
            for (float val : row) {
                if (val > 0) {
                    stats.addValue(val);
                }
            }
        }
        return stats;
    }

    private static float[][] divide(float[][] obsMatrix, float[][] eMatrix, double pseudoCount) {
        float[][] oeMatrix = new float[obsMatrix.length][obsMatrix[0].length];
        for (int i = 0; i < oeMatrix.length; i++) {
            for (int j = 0; j < oeMatrix[i].length; j++) {
                oeMatrix[i][j] = (float) ((obsMatrix[i][j] + pseudoCount) / (eMatrix[i][j] + pseudoCount));
            }
        }
        return oeMatrix;
    }

    private static float getP(float obs, float expected, double superDiagonal) {
        // P = (O - E)/(SD - E)
        return (float) ((obs - expected) / (superDiagonal - expected));
    }

    private static double getExpected(ContactRecord rec, ExpectedValueFunction df, int chrIndex) {
        int x = rec.getBinX();
        int y = rec.getBinY();
        int dist = Math.abs(x - y);
        return df.getExpectedValue(chrIndex, dist);
    }
}
