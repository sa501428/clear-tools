package cli.utils.general;

import javastraw.feature2D.Feature2D;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.List;

public class HiCValue {

    public static float getExactPixel(MatrixZoomData zd, Feature2D loop, int resolution, NormalizationType norm) {
        long binXStart = loop.getMidPt1() / resolution;
        long binYStart = loop.getMidPt2() / resolution;
        return getExactPixel(zd, binXStart, binXStart+2, binYStart, binYStart+2, norm);
    }

    public static float getExactPixel(MatrixZoomData zd, long binXStart, long binXEnd,
                                             long binYStart, long binYEnd, NormalizationType normalizationType) {
        List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd,
                normalizationType, false);
        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) {
                            int relativeX = (int) (rec.getBinX() - binXStart);
                            int relativeY = (int) (rec.getBinY() - binYStart);
                            if (relativeX == 0 && relativeY == 0) {
                                return rec.getCounts();
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
        System.err.println("Error retrieving pixel value");
        return Float.NaN;
    }


}
