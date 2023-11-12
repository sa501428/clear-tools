package cli.utils.apa;

import cli.utils.data.SparseContactMatrixWithMasking;
import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MultiAPAManager {

    private final Map<Integer, float[][]> upStreamAPAMatrices = new HashMap<>();
    private final Map<Integer, float[][]> downStreamAPAMatrices;

    public MultiAPAManager(List<Feature2D> loops, int window, int resolution, int matrixWidth,
                           SparseContactMatrixWithMasking scm, double[] vector, boolean isAggNorm, boolean dontUseOrientation) {

        Map<Integer, float[]> upStreamRowSums = new HashMap<>();
        Map<Integer, float[]> upStreamColSums = new HashMap<>();
        Map<Integer, float[]> downStreamRowSums;
        Map<Integer, float[]> downStreamColSums;
        if (dontUseOrientation) {
            downStreamAPAMatrices = upStreamAPAMatrices;
            downStreamRowSums = upStreamRowSums;
            downStreamColSums = upStreamColSums;
        } else {
            downStreamAPAMatrices = new HashMap<>();
            downStreamRowSums = new HashMap<>();
            downStreamColSums = new HashMap<>();
        }

        for (Feature2D loop : loops) {
            int upStreamBin = (int) (loop.getMidPt1() / resolution);
            int downStreamBin = (int) (loop.getMidPt2() / resolution);

            initialize(upStreamBin, upStreamAPAMatrices, upStreamRowSums, upStreamColSums, matrixWidth);
            initialize(downStreamBin, downStreamAPAMatrices, downStreamRowSums, downStreamColSums, matrixWidth);

            addToMatrix(upStreamAPAMatrices.get(upStreamBin), scm, loop, window, resolution, matrixWidth);
            addToMatrix(downStreamAPAMatrices.get(downStreamBin), scm, loop, window, resolution, matrixWidth);

            if (isAggNorm) {
                doAggregateNormalization(loop,
                        upStreamRowSums.get(upStreamBin),
                        upStreamColSums.get(upStreamBin), vector, resolution, window);
                doAggregateNormalization(loop,
                        downStreamRowSums.get(downStreamBin),
                        downStreamColSums.get(downStreamBin), vector, resolution, window);
            }
        }

        if (isAggNorm) {
            finalAggregateNormalization(upStreamAPAMatrices, upStreamRowSums, upStreamColSums);
            if (!dontUseOrientation) {
                finalAggregateNormalization(downStreamAPAMatrices, downStreamRowSums, downStreamColSums);
            }
        }
    }

    private void finalAggregateNormalization(Map<Integer, float[][]> matrices, Map<Integer, float[]> rowSums,
                                             Map<Integer, float[]> colSums) {
        for (int key : matrices.keySet()) {
            APADataExporter.normalizeBySum(rowSums.get(key));
            APADataExporter.normalizeBySum(colSums.get(key));
            float[][] matrix = APADataExporter.normedCopyFloats(matrices.get(key),
                    rowSums.get(key), colSums.get(key));
            matrices.put(key, matrix);
        }
    }

    private void initialize(int bin, Map<Integer, float[][]> matrices,
                            Map<Integer, float[]> rowSums, Map<Integer, float[]> colSums, int matrixWidth) {
        if (!matrices.containsKey(bin)) {
            matrices.put(bin, new float[matrixWidth][matrixWidth]);
            rowSums.put(bin, new float[matrixWidth]);
            colSums.put(bin, new float[matrixWidth]);
        }
    }

    public static void addToMatrix(float[][] output, SparseContactMatrixWithMasking scm,
                                   Feature2D loop, int window, int resolution, int matrixWidth) {
        int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
        int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
        scm.addLocalBoundedRegion(output, binXStart, binYStart, matrixWidth);
    }

    private void doAggregateNormalization(Feature2D loop, float[] rowSums, float[] colSums,
                                          double[] nv, int resolution, int window) {
        int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
        int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
        APAUtils.addLocalSums(rowSums, nv, binXStart);
        APAUtils.addLocalSums(colSums, nv, binYStart);
    }

    public List<AnchorAPAScore> getAnchorAPAScores(Chromosome chromosome, int resolution, boolean dontUseOrientation) {
        List<AnchorAPAScore> scores = new ArrayList<>();
        int width = 100;
        for (int bin : upStreamAPAMatrices.keySet()) {
            scores.add(new AnchorAPAScore(chromosome, resolution,
                    bin, width, "Forward_" + bin,
                    getAPAScore(upStreamAPAMatrices.get(bin), resolution),
                    true));
        }
        if (!dontUseOrientation) {
            for (int bin : downStreamAPAMatrices.keySet()) {
                scores.add(new AnchorAPAScore(chromosome, resolution,
                        bin, width, "Reverse_" + bin,
                        getAPAScore(downStreamAPAMatrices.get(bin), resolution),
                        false));
            }
        }
        return scores;
    }

    private float getAPAScore(float[][] data, int resolution) {
        float numerator;
        int n = data.length;
        if (resolution < 50) {
            int width = 50 / resolution;
            numerator = mean(data, n / 2 - width, n / 2 + width + 1, n / 2 - width, n / 2 + width + 1);
        } else {
            numerator = data[n / 2][n / 2];
        }

        float denom = mean(data, 3 * n / 4, n, 0, n / 4);

        return (numerator + 1) / (denom + 1);
    }

    private float mean(float[][] data, int r0, int r1, int c0, int c1) {
        float sum = 0;
        int count = 0;
        for (int r = r0; r < r1; r++) {
            for (int c = c0; c < c1; c++) {
                if (data[r][c] > 0) {
                    sum += data[r][c];
                    count++;
                }
            }
        }
        if (count > 0) {
            return sum / count;
        } else {
            return 0;
        }
    }
}
