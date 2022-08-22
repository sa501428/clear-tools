package cli.utils.sift;

import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;

import java.util.*;

public class FeatureUtils {

    public static Map<Region, Integer> addPointsToCountMap(Map<Integer, Set<SimpleLocation>> resToLocations,
                                                           int[] resolutions) {

        Map<Region, Integer> countMap = new HashMap<>();
        int base = 100;
        for (int res : resolutions) {
            if (res == base) {
                for (SimpleLocation location : resToLocations.get(res)) {
                    countMap.put(location.toRegion(res), 1);
                }
            } else {
                for (SimpleLocation location : resToLocations.get(res)) {
                    boolean updatesWereMade = updateMapForOverlap(location, countMap, res);
                    if (!updatesWereMade) {
                        countMap.put(location.toRegion(res), 1);
                    }
                }
            }
        }

        return countMap;
    }

    public static boolean updateMapForOverlap(SimpleLocation location, Map<Region, Integer> countMap, int res) {
        boolean updatesWereMade = false;
        for (Region region : countMap.keySet()) {
            if (region.containedBy(location, res)) {
                countMap.put(region, countMap.get(region) + 1);
                updatesWereMade = true;
            }
        }
        return updatesWereMade;
    }

    public static Set<Region> getPointsWithMoreThan(Map<Region, Integer> countsForRecord, int cutoff) {
        Set<Region> finalSet = new HashSet<>();
        for (Region record : countsForRecord.keySet()) {
            if (countsForRecord.get(record) >= cutoff) {
                finalSet.add(record);
            }
        }
        return finalSet;
    }

    public static List<Feature2D> convertToFeature2Ds(Set<ContactRecordBox> records, Chromosome c1) {
        List<Feature2D> features = new ArrayList<>();
        for (ContactRecordBox record : records) {
            features.add(record.toFeature2D(c1));
        }
        return features;
    }

    public static Set<SimpleLocation> toLocationsAndClear(Set<ContactRecord> enrichedRegions) {
        Set<SimpleLocation> locations = new HashSet<>();
        for (ContactRecord record : enrichedRegions) {
            locations.add(new SimpleLocation(record));
        }
        enrichedRegions.clear();
        return locations;
    }
}
