package cli.utils;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class RecapTools {

    public final static String FULL_DECAY = "DECAY_FULL";
    public final static String DECAY_A = "DECAY_A";
    public final static String DECAY_k = "DECAY_k";
    public final static String ROW_STD = "ROW_STD";
    public final static String COL_STD = "COL_STD";
    public final static String ROW_SUM = "ROW_SUM";
    public final static String COL_SUM = "COL_SUM";

    public static List<String> getCategories(boolean isLoopAnalysis) {
        List<String> categories = new ArrayList<>();
        String[] types = new String[]{"OBS_", "OE_"};
        String[] properties = new String[]{DECAY_A, DECAY_k, FULL_DECAY, ROW_SUM, COL_SUM, ROW_STD, COL_STD};
        //new String[]{"VAL", "STD_DEV", "KURTOSIS", "SKEWNESS", "MEAN_ENRICHMENT",
        //"MEDIAN_ENRICHMENT", "GEO_ENRICHMENT", "MAX_ENRICHMENT", "MIN_ENRICHMENT", "DECAY_A", "DECAY_k"};

        if (!isLoopAnalysis) {
            categories.add("PRESENCE");
            categories.add("PRESENCE_INF");
            return categories;
        }

        for (String stem : types) {
            for (String property : properties) {
                categories.add(stem + property);
            }
        }
        return categories;
    }

    public static Map<String, String> getStats(float[][] obsMatrix, float[][] eMatrix,
                                               int window, double superDiagonal, float pseudoCount, boolean isDeepLoopAnalysis) {

        Map<String, String> loopAttributes = new HashMap<>();
        float obs = obsMatrix[window][window];
        addAttributes(loopAttributes, "OBS_", obsMatrix, window, isDeepLoopAnalysis);

        if (eMatrix != null) {
            float[][] oeMatrix = divide(obsMatrix, eMatrix, pseudoCount);
            if (!isDeepLoopAnalysis) {
                float expected = eMatrix[window][window];
                loopAttributes.put("PRESENCE", String.valueOf(getP(obs, expected, superDiagonal)));
                loopAttributes.put("PRESENCE_INF", String.valueOf(getP(obs, pseudoCount, superDiagonal)));
            }
            addAttributes(loopAttributes, "OE_", oeMatrix, window, isDeepLoopAnalysis);
        }

        return loopAttributes;
    }

    private static void addAttributes(Map<String, String> attributes, String stem, float[][] matrix, int window,
                                      boolean isDeepLoopAnalysis) {
        /*
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
        */

        if (isDeepLoopAnalysis) {
            addMatrixSums(matrix, attributes, stem);

            // calculates the decay vetor
            float[] manhattanDecay = calculateDecay(matrix, window);
            addRegressionStats(manhattanDecay, attributes, stem);
        }
    }

    private static void addMatrixSums(float[][] matrix, Map<String, String> attributes, String stem) {

        // TODO @Justin calculate the row and column sums of matrix
        float[] rowSums = new float[1];
        float[] colSums = new float[1];
        //TODO normalize them by the respective middle values of the vector

        // TODO get the std dev from the row and column sums
        // save them

        attributes.put(stem + ROW_SUM, convertVectorToString(rowSums));
        attributes.put(stem + COL_SUM, convertVectorToString(colSums));
    }

    private static void addRegressionStats(float[] decay, Map<String, String> attributes, String stem) {
        SimpleRegression regression = new SimpleRegression();
        for (int i = 0; i < decay.length; i++) {
            regression.addData(i, Math.log(decay[i]));
        }
        // y = Ae^(kx)  for each (xi, yi) --> wi = log(yi)
        // solve w = a + bx
        // e^a = A  ;  b = k
        attributes.put(stem + DECAY_A, String.valueOf(Math.exp(regression.getIntercept())));
        attributes.put(stem + DECAY_k, String.valueOf(regression.getSlope()));
        attributes.put(stem + FULL_DECAY, convertVectorToString(decay));
    }

    private static String convertVectorToString(float[] vector) {
        StringJoiner joiner = new StringJoiner(",");
        for (float val : vector) {
            joiner.add("" + val);
        }
        return joiner.toString();
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

    public static void exportAllMatrices(Feature2DList refinedLoops, String[] names, File outFolder,
                                         boolean isDeepLoopAnalysis, int window) {
        int n = refinedLoops.getNumTotalFeatures();
        int m = names.length;

        List<String> categories = getCategories(isDeepLoopAnalysis);

        List<float[][]> outputs = new ArrayList<>();
        for (int k = 0; k < categories.size(); k++) {
            if (categories.get(k).contains(FULL_DECAY)) {
                outputs.add(new float[n][m * (window + 1)]);
            } else if (categories.get(k).contains(ROW_SUM) || categories.get(k).contains(COL_SUM)) {
                outputs.add(new float[n][m * (2 * window + 1)]);
            } else {
                outputs.add(new float[n][m]);
            }
        }

        AtomicInteger loopIndex = new AtomicInteger(0);
        refinedLoops.processLists((s, list) -> {
            for (Feature2D loop : list) {
                int currIndex = loopIndex.getAndIncrement();
                for (int k = 0; k < categories.size(); k++) {
                    for (int w = 0; w < names.length; w++) {
                        String key = names[w] + "_" + categories.get(k);
                        if (categories.get(k).contains(FULL_DECAY)) {
                            fillInVector(outputs.get(k), loop.getAttribute(key), currIndex, w, window + 1);
                        } else if (categories.get(k).contains(ROW_SUM) || categories.get(k).contains(COL_SUM)) {
                            fillInVector(outputs.get(k), loop.getAttribute(key), currIndex, w, 2 * window + 1);
                        } else {
                            outputs.get(k)[currIndex][w] = Float.parseFloat(loop.getAttribute(key));
                        }
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

    private static void fillInVector(float[][] matrix, String values, int r, int fileIndex, int vLength) {
        String[] vals = values.split(",");
        for (int z = 0; z < vLength; z++) {
            int realColumn = z + (vLength * fileIndex);
            matrix[r][realColumn] = Float.parseFloat(vals[z]);
        }
    }
}
