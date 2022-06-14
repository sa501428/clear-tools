package cli.utils;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class RecapTools {

    public static List<String> getCategories(boolean useOE) {
        List<String> categories = new ArrayList<>();

        String[] properties = new String[]{"VAL", "STD_DEV", "KURTOSIS", "SKEWNESS", "MEAN_ENRICHMENT",
                "MEDIAN_ENRICHMENT", "GEO_ENRICHMENT", "MAX_ENRICHMENT", "MIN_ENRICHMENT", "DECAY_A", "DECAY_k"};
        String[] types = new String[]{"OBS_"};

        if (useOE) {
            categories.add("PRESENCE");
            categories.add("PRESENCE_INF");
            types = new String[]{"OBS_", "OE_"};
        }
        for (String stem : types) {
            for (String property : properties) {
                categories.add(stem + property);
            }
        }
        return categories;
    }

    public static Map<String, String> getStats(float[][] obsMatrix, float[][] eMatrix,
                                               int window, double superDiagonal, float pseudoCount) {

        Map<String, String> loopAttributes = new HashMap<>();
        float obs = obsMatrix[window][window];
        addAttributes(loopAttributes, "OBS_", obsMatrix, window);

        if (eMatrix != null) {
            float[][] oeMatrix = divide(obsMatrix, eMatrix, pseudoCount);
            float expected = eMatrix[window][window];
            loopAttributes.put("PRESENCE", String.valueOf(getP(obs, expected, superDiagonal)));
            loopAttributes.put("PRESENCE_INF", String.valueOf(getP(obs, pseudoCount, superDiagonal)));
            addAttributes(loopAttributes, "OE_", oeMatrix, window);
        }

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

    public static void exportAllMatrices(Feature2DList refinedLoops, String[] names, File outFolder, boolean useOE) {
        int n = refinedLoops.getNumTotalFeatures();
        int m = names.length;

        List<String> categories = getCategories(useOE);

        List<float[][]> outputs = new ArrayList<>();
        for (int k = 0; k < categories.size(); k++) {
            outputs.add(new float[n][m]);
        }

        AtomicInteger loopIndex = new AtomicInteger(0);
        refinedLoops.processLists((s, list) -> {
            for (Feature2D loop : list) {
                int currIndex = loopIndex.getAndIncrement();
                for (int k = 0; k < categories.size(); k++) {
                    for (int w = 0; w < names.length; w++) {
                        String key = names[w] + "_" + categories.get(k);
                        outputs.get(k)[currIndex][w] = Float.parseFloat(loop.getAttribute(key));
                    }
                }
            }
        });

        for (int k = 0; k < categories.size(); k++) {
            MatrixTools.saveMatrixTextNumpy((new File(outFolder, categories.get(k) + ".npy")).getAbsolutePath(),
                    outputs.get(k));
        }
        outputs.clear();
    }
}
