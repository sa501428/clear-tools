package cli.utils.data;

import javastraw.reader.block.ContactRecord;

import java.util.HashMap;
import java.util.Map;

public abstract class SparseContactRecordStorage {

    protected final Map<Integer, Map<Integer, ContactRecord>> data = new HashMap<>();

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
