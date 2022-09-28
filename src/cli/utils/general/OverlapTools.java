package cli.utils.general;

import javastraw.feature2D.Feature2D;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class OverlapTools {

    public static Set<Feature2D> getMatchesWithOverlap(Feature2D pixel, List<Feature2D> features, int buffer) {
        Set<Feature2D> pixelList = new HashSet<>();
        for (Feature2D px : features) {
            if (hasOverlap(px, pixel, buffer)) {
                pixelList.add(px);
            }
        }
        return pixelList;
    }

    public static Set<Feature2D> getExactMatches(Feature2D pixel, List<Feature2D> features) {
        Set<Feature2D> pixelList = new HashSet<>();
        for (Feature2D px : features) {
            if (isExact(px, pixel)) {
                pixelList.add(px);
            }
        }
        return pixelList;
    }

    private static boolean isExact(Feature2D px, Feature2D px2) {
        return px.getStart1() == px2.getStart1()
                && px.getStart2() == px2.getStart2()
                && px.getEnd1() == px2.getEnd1()
                && px.getEnd2() == px2.getEnd2();
    }

    private static boolean hasOverlap(Feature2D px1, Feature2D original, int buffer) {
        return getWidth(px1.getStart1(), px1.getEnd1(), original.getStart1() - buffer, original.getEnd1() + buffer) *
                getWidth(px1.getStart2(), px1.getEnd2(), original.getStart2() - buffer, original.getEnd2() + buffer) > 0;
    }

    private static int getWidth(long p1, long p2, long g1, long g2) {
        long start = Math.max(p1, g1);
        long end = Math.min(p2, g2);
        return (int) Math.max(0, end - start);
    }
}
