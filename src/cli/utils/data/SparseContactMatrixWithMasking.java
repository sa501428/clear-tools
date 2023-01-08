package cli.utils.data;

import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class SparseContactMatrixWithMasking {

    private final Map<Integer, Map<Integer, ContactRecord>> data = new HashMap<>();

    public SparseContactMatrixWithMasking(MatrixZoomData zd, Collection<Feature2D> loops, int resolution,
                                          int window, int matrixWidth, NormalizationType norm) {
        Set<Integer> allIndices = getIndices(loops, resolution, window, matrixWidth);

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (allIndices.contains(cr.getBinX()) && allIndices.contains(cr.getBinY())) {
                    if (!data.containsKey(cr.getBinX())) {
                        data.put(cr.getBinX(), new HashMap<>());
                    }
                    data.get(cr.getBinX()).put(cr.getBinY(), cr);
                }
            }
        }
    }

    private static Set<Integer> getIndices(Collection<Feature2D> loops, int resolution, int window, int matrixWidth) {
        Set<Integer> indices = new HashSet<>();
        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);

            for (int i = 0; i < matrixWidth; i++) {
                indices.add(i + binXStart);
                indices.add(i + binYStart);
            }
        }
        return indices;
    }

    public float getContact(int x, int y) {
        if (data.containsKey(x)) {
            if (data.get(x).containsKey(y)) {
                return data.get(x).get(y).getCounts();
            }
        }
        return 0;
    }

    public void eraseAll() {
        for (Map<Integer, ContactRecord> map : data.values()) {
            map.clear();
        }
        data.clear();
    }

    public float[][] getRegion(int binXStart, int binYStart, int binXEnd, int binYEnd) {
        int numRows = binXEnd - binXStart;
        int numCols = binYEnd - binYStart;
        float[][] matrix = new float[numRows][numCols];
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                matrix[i][j] = getContact(i + binXStart, j + binYStart);
            }
        }
        return matrix;
    }

    public void addLocalBoundedRegion(float[][] output, int binXStart, int binYStart, int matrixWidth) {
        for (int r = 0; r < matrixWidth; r++) {
            for (int c = 0; c < matrixWidth; c++) {
                output[r][c] += getContact(r + binXStart, c + binYStart);
            }
        }
    }
}
