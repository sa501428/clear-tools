package cli.enhance;

import javastraw.feature2D.Feature2D;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.List;

public class EnhanceUtils {

    public static void addLocalBoundedRegion(int[][] matrix, MatrixZoomData zd, long binXStart, long binYStart,
                                             int window, int matrixWidth, final Object key) {

        long binXEnd = binXStart + (window + 1);
        long binYEnd = binYStart + (window + 1);

        List<Block> blocks;
        synchronized (key) {
            blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd,
                    NormalizationHandler.NONE, false);
        }
        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) {
                            // only called for small regions - should not exceed int
                            int relativeX = (int) (rec.getBinX() - binXStart);
                            int relativeY = (int) (rec.getBinY() - binYStart);
                            if (relativeX >= 0 && relativeX < matrixWidth) {
                                if (relativeY >= 0 && relativeY < matrixWidth) {
                                    matrix[relativeX][relativeY] += rec.getCounts();
                                }
                            }
                        }
                    }
                }
            }
        }
        // force cleanup
        blocks.clear();
        blocks = null;
        //System.gc();
    }
}
