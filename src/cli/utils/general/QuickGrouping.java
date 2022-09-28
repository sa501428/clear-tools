package cli.utils.general;

import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;

import java.util.*;

public class QuickGrouping {

    public static Map<SimpleLocation, List<Feature2D>> groupNearbyRecords(List<Feature2D> points, int scalar) {
        return groupNearbyRecords(new HashSet<>(points), scalar);
    }

    public static Map<SimpleLocation, List<Feature2D>> groupNearbyRecords(Set<Feature2D> initialPoints, int scalar) {
        Map<SimpleLocation, List<Feature2D>> locationMap = new HashMap<>();
        for (Feature2D feature : initialPoints) {
            SimpleLocation region = new SimpleLocation((int) (feature.getMidPt1() / scalar), (int) (feature.getMidPt2() / scalar));
            if (locationMap.containsKey(region)) {
                locationMap.get(region).add(feature);
            } else {
                List<Feature2D> values = new ArrayList<>();
                values.add(feature);
                locationMap.put(region, values);
            }
        }
        return locationMap;
    }

    public static Map<SimpleLocation, List<Feature2D>> groupNearbyRecordsWithOverlap(List<Feature2D> points, int scalar) {
        return groupNearbyRecordsWithOverlap(new HashSet<>(points), scalar);
    }

    public static Map<SimpleLocation, List<Feature2D>> groupNearbyRecordsWithOverlap(Set<Feature2D> initialPoints, int scalar) {
        Map<SimpleLocation, List<Feature2D>> locationMap = new HashMap<>();
        for (Feature2D feature : initialPoints) {
            int setR = (int) (feature.getMidPt1() / scalar);
            int setC = (int) (feature.getMidPt2() / scalar);

            // add to neighbors as well;
            for (int dr = -1; dr < 2; dr++) {
                for (int dc = -1; dc < 2; dc++) {
                    SimpleLocation region = new SimpleLocation(setR + dr, setC + dc);
                    if (locationMap.containsKey(region)) {
                        locationMap.get(region).add(feature);
                    } else {
                        List<Feature2D> values = new ArrayList<>();
                        values.add(feature);
                        locationMap.put(region, values);
                    }
                }
            }
        }
        return locationMap;
    }
}
