package cli.utils.sift;

import javastraw.reader.block.ContactRecord;

import java.util.*;

public class NMSUtils {


    public static void filterOutByOverlap(Set<ContactRecord> initialPoints, Set<SimpleLocation> regions, int scalar) {
        Set<ContactRecord> toRemove = new HashSet<>();
        Map<SimpleLocation, List<ContactRecord>> locationMap = new HashMap<>();
        for (ContactRecord cr : initialPoints) {
            SimpleLocation region = new SimpleLocation(cr.getBinX() / scalar, cr.getBinY() / scalar);
            if (regions.contains(region)) {
                if (locationMap.containsKey(region)) {
                    locationMap.get(region).add(cr);
                } else {
                    List<ContactRecord> values = new ArrayList<>();
                    values.add(cr);
                    locationMap.put(region, values);
                }
            } else {
                toRemove.add(cr);
            }
        }

        nonMaxSuppressionInGroup(initialPoints, toRemove, locationMap);
    }

    public static void filterOutByOverlap(Set<ContactRecord> initialPoints, int scalar) {
        Set<ContactRecord> toRemove = new HashSet<>();
        Map<SimpleLocation, List<ContactRecord>> locationMap = groupNearbyRecords(initialPoints, scalar);
        nonMaxSuppressionInGroup(initialPoints, toRemove, locationMap);
    }

    static Map<SimpleLocation, List<ContactRecord>> groupNearbyRecords(Set<ContactRecord> initialPoints, int scalar) {
        Map<SimpleLocation, List<ContactRecord>> locationMap = new HashMap<>();
        for (ContactRecord cr : initialPoints) {
            SimpleLocation region = new SimpleLocation(cr.getBinX() / scalar, cr.getBinY() / scalar);
            if (locationMap.containsKey(region)) {
                locationMap.get(region).add(cr);
            } else {
                List<ContactRecord> values = new ArrayList<>();
                values.add(cr);
                locationMap.put(region, values);
            }
        }
        return locationMap;
    }

    public static void nonMaxSuppressionInGroup(Set<ContactRecord> initialPoints, Set<ContactRecord> toRemove, Map<SimpleLocation, List<ContactRecord>> locationMap) {
        for (List<ContactRecord> overlaps : locationMap.values()) {
            if (overlaps.size() > 1) {
                float max = getMax(overlaps);
                for (ContactRecord cr2 : overlaps) {
                    if (cr2.getCounts() < max) {
                        toRemove.add(cr2);
                    }
                }
            }
        }
        initialPoints.removeAll(toRemove);
        locationMap.clear();
    }


    public static float getMax(List<ContactRecord> records) {
        float maxVal = records.get(0).getCounts();
        for (ContactRecord cr : records) {
            if (cr.getCounts() > maxVal) {
                maxVal = cr.getCounts();
            }
        }
        return maxVal;
    }

    public static boolean inRegions(ContactRecord cr, Set<SimpleLocation> regions, int scalar) {
        SimpleLocation region = new SimpleLocation(cr.getBinX() / scalar, cr.getBinY() / scalar);
        return regions.contains(region);
    }

}
