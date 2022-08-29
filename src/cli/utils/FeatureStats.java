package cli.utils;

import javastraw.feature2D.Feature2D;

import java.util.LinkedList;

public class FeatureStats {

    public static long minStart1(LinkedList<Feature2D> group) {
        long minStart1 = group.peek().getStart1();
        for (Feature2D feature : group) {
            minStart1 = Math.min(minStart1, feature.getStart1());
        }
        return minStart1;
    }

    public static long minStart2(LinkedList<Feature2D> group) {
        long minStart2 = group.peek().getStart2();
        for (Feature2D feature : group) {
            minStart2 = Math.min(minStart2, feature.getStart2());
        }
        return minStart2;
    }


    public static long maxEnd1(LinkedList<Feature2D> group) {
        long maxEnd1 = group.peek().getEnd1();
        for (Feature2D feature : group) {
            maxEnd1 = Math.max(maxEnd1, feature.getEnd1());
        }
        return maxEnd1;
    }

    public static long maxEnd2(LinkedList<Feature2D> group) {
        long maxEnd2 = group.peek().getEnd2();
        for (Feature2D feature : group) {
            maxEnd2 = Math.max(maxEnd2, feature.getEnd2());
        }
        return maxEnd2;
    }
}
