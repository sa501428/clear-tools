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
            int binXStart = Math.max(0, (int) ((loop.getMidPt1() / resolution) - window));
            int binYStart = Math.max(0, (int) ((loop.getMidPt2() / resolution) - window));
            int binXEnd = Math.min(binXStart + matrixWidthL, Integer.MAX_VALUE);
            int binYEnd = Math.min(binYStart + matrixWidthL, Integer.MAX_VALUE);
            Bounds2DInfo bounds2DInfo = new Bounds2DInfo(binXStart, binXEnd, binYStart, binYEnd);

            int compressedBinXStart = binXStart / compression;
            int compressedBinYStart = binYStart / compression;
            int compressedBinXEnd = binXEnd / compression;
            int compressedBinYEnd = binYEnd / compression;

            for (int r = compressedBinXStart; r <= compressedBinXEnd; r++) {
                for (int c = compressedBinYStart; c <= compressedBinYEnd; c++) {
                    Point2D point = new Point2D(r, c);
                    loopListAsMap.computeIfAbsent(point, k -> new LinkedList<>()).add(bounds2DInfo);
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
        Point2D point = new Point2D(cX, cY);
        List<Bounds2DInfo> boundsList = loopListAsMap.get(point);
        if (boundsList != null) {
            for (Bounds2DInfo bounds2DInfo : boundsList) {
                if (bounds2DInfo.contains(binX, binY)) {
                    return true;
                }
            }
        }
        return false;
    }

    public List<Bounds2DInfo> get(int binX, int binY) {
        int cX = binX / compression;
        int cY = binY / compression;
        Point2D point = new Point2D(cX, cY);
        List<Bounds2DInfo> result = new LinkedList<>();
        List<Bounds2DInfo> boundsList = loopListAsMap.get(point);
        if (boundsList != null) {
            for (Bounds2DInfo bounds2DInfo : boundsList) {
                if (bounds2DInfo.contains(binX, binY)) {
                    result.add(bounds2DInfo);
                }
            }
        }
        return result;
    }
}
