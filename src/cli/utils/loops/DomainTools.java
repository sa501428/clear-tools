package cli.utils.loops;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class DomainTools {
    public static Feature2DList[] createLoops(Feature2DList domains) {
        Feature2DList corners = new Feature2DList();
        Feature2DList centers = new Feature2DList();
        domains.processLists((s, list) -> {
            List<Feature2D> outCorners = new ArrayList<>();
            List<Feature2D> outCenters = new ArrayList<>();
            for (Feature2D feature : list) {
                outCorners.add(getCorner(feature));
                outCenters.add(getCenter(feature));
            }
            corners.addByKey(s, outCorners);
            centers.addByKey(s, outCenters);
        });
        return new Feature2DList[]{corners, centers};
    }

    private static Feature2D getCenter(Feature2D feature) {
        long width = feature.getEnd1() - feature.getStart1();
        long anchor1 = (long) (feature.getStart1() + 0.25 * width);
        long anchor2 = (long) (feature.getStart1() + 0.75 * width);

        long start1 = anchor1 - 100;
        long end1 = anchor1 + 100;
        long start2 = anchor2 - 100;
        long end2 = anchor2 + 100;

        return new Feature2D(Feature2D.FeatureType.PEAK, feature.getChr1(), start1, end1,
                feature.getChr2(), start2, end2, Color.BLACK, new HashMap<>());
    }

    private static Feature2D getCorner(Feature2D feature) {
        long start1 = feature.getStart1() - 100;
        long end1 = feature.getStart1() + 100;

        long start2 = feature.getEnd2() - 100;
        long end2 = feature.getEnd2() + 100;

        return new Feature2D(Feature2D.FeatureType.PEAK, feature.getChr1(), start1, end1,
                feature.getChr2(), start2, end2, Color.BLACK, new HashMap<>());
    }
}
