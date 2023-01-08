package cli.utils.loops;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class OffsetTools {
    public static Feature2DList createLoops(Feature2DList loops, int offset,
                                            boolean addOffset, boolean subtractOffset) {
        Feature2DList output = new Feature2DList();
        loops.processLists((s, list) -> {
            List<Feature2D> out = new ArrayList<>();
            for (Feature2D feature : list) {
                if (subtractOffset) {
                    Feature2D upStream = getOffsetFeature(feature, -offset);
                    if (upStream != null) {
                        out.add(upStream);
                    }
                }
                if (addOffset) {
                    Feature2D downStream = getOffsetFeature(feature, offset);
                    if (downStream != null) {
                        out.add(downStream);
                    }
                }
            }
            output.addByKey(s, out);
        });
        return output;
    }

    private static Feature2D getOffsetFeature(Feature2D feature, int offset) {
        long start1 = feature.getStart1() + offset;
        long start2 = feature.getStart2() + offset;
        long end1 = feature.getEnd1() + offset;
        long end2 = feature.getEnd2() + offset;
        if (start1 > 0 && start2 > 0) {
            return new Feature2D(Feature2D.FeatureType.PEAK, feature.getChr1(), start1, end1,
                    feature.getChr2(), start2, end2, Color.BLACK, new HashMap<>());
        }
        return null;
    }
}
