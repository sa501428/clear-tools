package cli.utils.recap;

import cli.utils.flags.Utils;
import javastraw.expected.ExpectedModel;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.tools.MatrixTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.io.File;
import java.util.*;

// import static cli.utils.apa.APADataExporter.normalizeBySum;

public class RecapTools {

    public final static String FULL_DECAY = "DECAY_FULL";
    public final static String DECAY_A = "DECAY_A";
    public final static String DECAY_k = "DECAY_k";
    public final static String ROW_SUM = "ROW_SUM";
    public final static String COL_SUM = "COL_SUM";
    public final static String AMP_ROW = "AMP_ROW";
    public final static String AMP_COL = "AMP_COL";
    public final static String SPREAD_ROW = "SPREAD_ROW";
    public final static String SPREAD_COL = "SPREAD_COL";

    public static List<String> getCategories(boolean isLoopAnalysis) {
        List<String> categories = new ArrayList<>();
        String[] types = new String[]{"OBS_", "OE_"};
        String[] properties = new String[]{DECAY_A, DECAY_k, FULL_DECAY, ROW_SUM, COL_SUM, AMP_ROW, AMP_COL,
                SPREAD_ROW, SPREAD_COL};
        //new String[]{"VAL", "STD_DEV", "KURTOSIS", "SKEWNESS", "MEAN_ENRICHMENT",
        //"MEDIAN_ENRICHMENT", "GEO_ENRICHMENT", "MAX_ENRICHMENT", "MIN_ENRICHMENT", "DECAY_A", "DECAY_k"};

        if (!isLoopAnalysis) {
            categories.add("PRESENCE");
            return categories;
        }

        for (String stem : types) {
            for (String property : properties) {
                categories.add(stem + property);
            }
        }
        return categories;
    }

    public static Map<String, String> getStats(float[][] obsMatrix,
                                               int window, float pseudoCount,
                                               boolean isDeepLoopAnalysis, ExpectedModel polynomial,
                                               Feature2D loop, int resolution) {

        Map<String, String> attributes = new HashMap<>();

        if (isDeepLoopAnalysis) {
            addMatrixSums(obsMatrix, attributes, "OBS_");
            float[] manhattanDecay = calculateDecay(obsMatrix, window);
            addRegressionStats(manhattanDecay, attributes, "OBS_");

            float[][] eMatrix = new float[obsMatrix.length][obsMatrix.length];
            Utils.fillInExpectedMatrix(eMatrix, loop, obsMatrix.length, polynomial, resolution, window);

            float[][] oeMatrix = divide(obsMatrix, eMatrix, pseudoCount);
            addMatrixSums(oeMatrix, attributes, "OE_");
            manhattanDecay = calculateDecay(oeMatrix, window);
            addRegressionStats(manhattanDecay, attributes, "OE_");
        } else {
            float obs = obsMatrix[window][window];
            double p = getPresenceFrom(polynomial, loop, resolution, obs);
            attributes.put("PRESENCE", String.valueOf(p));
        }

        return attributes;
    }

    private static double getPresenceFrom(ExpectedModel model, Feature2D loop, int resolution, float counts) {
        int midX = (int) (loop.getMidPt1() / resolution);
        int midY = (int) (loop.getMidPt2() / resolution);
        int dist = Math.abs(midX - midY);
        return model.getPercentContact(dist, counts);
    }

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

    private static void addMatrixSums(float[][] matrix, Map<String, String> attributes, String stem) {

        // initializes variables
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        float[] rowSum = new float[numRows];
        float[] colSum = new float[numCols];

        // calculates row and col sums, and produces normalized row and col sums
        populateRowColSums(matrix, rowSum, colSum);
        float[] normalizedRowSum = getNormalizedSum(rowSum);
        float[] normalizedColSum = getNormalizedSum(colSum);

        // calculate std deviations (amplitude and spread)
        // getStd functions return doubles
        double ampRow = getAmplitudeStdDev(rowSum);
        double ampCol = getAmplitudeStdDev(colSum);
        double spreadRow = getSpreadStdDev(normalizedRowSum);
        double spreadCol = getSpreadStdDev(normalizedColSum);

        // adds these metrics to the attributes object
        attributes.put(stem + ROW_SUM, convertVectorToString(rowSum));
        attributes.put(stem + COL_SUM, convertVectorToString(colSum));
        attributes.put(stem + AMP_ROW, String.valueOf(ampRow));
        attributes.put(stem + AMP_COL, String.valueOf(ampCol));
        attributes.put(stem + SPREAD_ROW, String.valueOf(spreadRow));
        attributes.put(stem + SPREAD_COL, String.valueOf(spreadCol));
    }

    private static double getAmplitudeStdDev(float[] sumVector) {
        // helper fxn
        // moves rowSum elements to rowStats container to use getStandardDeviation method
        // note: amplitude std dev is like the "vertical" std dev
        DescriptiveStatistics rowStats = new DescriptiveStatistics();
        for (float val : sumVector) {
            rowStats.addValue(val);
        }
        return rowStats.getStandardDeviation();
    }

    private static double getSpreadStdDev(float[] pdf) {
        // helper fxn
        // NOTE: assumption here is that the "average" (mu) is the center of the window, may not always be true
        // note 2: spread std dev is like the "horizontal" std dev
        int mu = pdf.length / 2;
        double total = 0;
        for(int k = 0; k < pdf.length; k++){
            total += (k - mu)*(k - mu)*pdf[k];
        }
        return Math.sqrt(total);
    }

    private static float[] getNormalizedSum(float[] vector) {
        // output: normalized vector, each element of vector is divided by the sum of all the elements
        // this produces a vector whose elements all sum to 1
        // aka, a probability distribution
        float sum = getSum(vector);
        if (sum != 0) {
            for (int k = 0; k < vector.length; k++) {
                vector[k] /= sum;
            }
        }
        return vector;
    }

    private static float getSum(float[] vector) {
        // helper fxn for helper fxn getNormalizedSum
        // output: sum of elements in vector
        float sum = 0;
        for (float v : vector) {
            sum += v;
        }
        return sum;
    }


    private static void populateRowColSums(float[][] matrix, float[] rowSum, float[] colSum) {
        for (int row = 0; row < rowSum.length; row++) {
            for (int col = 0; col < colSum.length; col++) {
                rowSum[row] += matrix[row][col];
                colSum[col] += matrix[row][col];
            }
        }
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

    private static double getExpected(ContactRecord rec, ExpectedValueFunction df, int chrIndex) {
        int x = rec.getBinX();
        int y = rec.getBinY();
        int dist = Math.abs(x - y);
        return df.getExpectedValue(chrIndex, dist);
    }

    public static void exportAllMatrices(Chromosome[] chromosomes, Feature2DList refinedLoops, String[] names,
                                         File outFolder, boolean isDeepLoopAnalysis, int window) {
        int n = refinedLoops.getNumTotalFeatures();
        int m = names.length;

        List<String> categories = getCategories(isDeepLoopAnalysis);

        List<float[][]> outputs = new ArrayList<>();
        for (String category : categories) {
            if (category.contains(FULL_DECAY)) {
                outputs.add(new float[n][m * (window + 1)]);
            } else if (category.contains(ROW_SUM) || category.contains(COL_SUM)) {
                outputs.add(new float[n][m * (2 * window + 1)]);
            } else {
                outputs.add(new float[n][m]);
            }
        }

        int currIndex = 0;
        for (Chromosome chrom : chromosomes) {
            List<Feature2D> loops = refinedLoops.get(chrom.getIndex(), chrom.getIndex());
            Collections.sort(loops);
            for (Feature2D loop : loops) {
                //int currIndex = loopIndex.getAndIncrement();
                for (int k = 0; k < categories.size(); k++) {
                    for (int w = 0; w < names.length; w++) {
                        String key = names[w] + "_" + categories.get(k);
                        if (categories.get(k).contains(FULL_DECAY)) {
                            fillInVector(outputs.get(k), loop.getAttribute(key), currIndex, w, window + 1);
                        } else if (categories.get(k).contains(ROW_SUM) || categories.get(k).contains(COL_SUM)) {
                            fillInVector(outputs.get(k), loop.getAttribute(key), currIndex, w, 2 * window + 1);
                        } else {
                            try {
                                outputs.get(k)[currIndex][w] = Float.parseFloat(loop.getAttribute(key));
                            } catch (Exception e) {
                                System.err.println("failed key " + key);
                                e.printStackTrace();
                            }
                        }
                    }
                }
                currIndex++;
            }
        }

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
