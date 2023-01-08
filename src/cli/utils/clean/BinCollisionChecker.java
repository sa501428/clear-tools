package cli.utils.clean;

import javastraw.feature2D.Feature2D;

import java.util.*;

public class BinCollisionChecker {
    public static Set<Feature2D> ensureOnlyOneLoopPerBin(Set<Feature2D> loops, int resolution) {
        Map<String, List<Feature2D>> collisions = new HashMap<>();
        for (Feature2D feature : loops) {
            String key = getKey(feature, resolution);
            if (!collisions.containsKey(key)) {
                collisions.put(key, new LinkedList<>());
            }
            collisions.get(key).add(feature);
        }
        Set<Feature2D> goodLoops = new HashSet<>();
        for (List<Feature2D> features : collisions.values()) {
            if (features.size() == 1) {
                goodLoops.add(features.get(0));
            }
        }
        collisions.clear();
        return goodLoops;
    }

    private static String getKey(Feature2D feature, int resolution) {
        return (int) (feature.getMidPt1() / resolution) + "_" + (int) (feature.getMidPt2() / resolution);
    }
}
