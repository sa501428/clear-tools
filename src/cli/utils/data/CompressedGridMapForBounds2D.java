package cli.utils.data;

import javastraw.feature2D.Feature2D;

import java.util.*;

public class CompressedGridMapForBounds2D {

    private final int resolution;
    private final int window;
    private final int matrixWidthL;
    private final int compression = 1000;
    private final Map<Point2D, List<Bounds2DInfo>> loopListAsMap = new HashMap<>();

    public CompressedGridMapForBounds2D(List<Feature2D> loops, int resolution, int window, int matrixWidthL,
                                        Set<Integer> allRowIndices, Set<Integer> allColIndices) {

        this.resolution = resolution;
        this.window = window;
        this.matrixWidthL = matrixWidthL;


        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);
            int binXEnd = binXStart + matrixWidthL;
            int binYEnd = binYStart + matrixWidthL;
            Bounds2DInfo bounds2DInfo = new Bounds2DInfo(binXStart, binXEnd, binYStart, binYEnd);

            int compressedBinXStart = binXStart / compression;
            int compressedBinYStart = binYStart / compression;
            int compressedBinXEnd = binXEnd / compression;
            int compressedBinYEnd = binYEnd / compression;

            for (int r = compressedBinXStart; r <= compressedBinXEnd; r++) {
                for (int c = compressedBinYStart; c <= compressedBinYEnd; c++) {
                    if (!loopListAsMap.containsKey(new Point2D(r, c))) {
                        loopListAsMap.put(new Point2D(r, c), new LinkedList<>());
                    }
                    loopListAsMap.get(new Point2D(r, c)).add(bounds2DInfo);
                }
            }

            for (int r = binXStart; r < binXEnd; r++) {
                allRowIndices.add(r);
            }

            for (int c = binYStart; c < binYEnd; c++) {
                allColIndices.add(c);
            }
        }
    }

    public boolean contains(int binX, int binY) {
        int cX = binX / compression;
        int cY = binY / compression;
        for (Bounds2DInfo bounds2DInfo : loopListAsMap.get(new Point2D(cX, cY))) {
            if (bounds2DInfo.contains(binX, binY)) {
                return true;
            }
        }
        return false;
    }

    public List<Bounds2DInfo> get(int binX, int binY) {
        int cX = binX / compression;
        int cY = binY / compression;
        List<Bounds2DInfo> result = new LinkedList<>();
        for (Bounds2DInfo bounds2DInfo : loopListAsMap.get(new Point2D(cX, cY))) {
            if (bounds2DInfo.contains(binX, binY)) {
                result.add(bounds2DInfo);
            }
        }
        return result;
    }


    // Map<Integer, List<BoundsInfo>>


}
