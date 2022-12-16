package cli.utils.sift.collapse;

import cli.clt.loops.Sift;
import cli.utils.sift.ContactRecordBox;
import cli.utils.sift.SimpleLocation;
import javastraw.reader.block.ContactRecord;

import java.util.*;

public class MultiResCentroidCollapser {

    public static Set<ContactRecordBox> coalesce(Set<ContactRecordBox> regions, int minResolutions) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        Map<SimpleLocation, LinkedList<ContactRecordBox>> map = groupNearbyRecords(regions, 5000000);
        Set<ContactRecordBox> coalesced = new HashSet<>();

        for (LinkedList<ContactRecordBox> records : map.values()) {
            sortInPlace(records);

            while (!records.isEmpty()) {
                ContactRecordBox pixel = records.pollFirst();
                if (pixel != null) {
                    records.remove(pixel);

                    int gX1 = pixel.getGenomeX1();
                    int gY1 = pixel.getGenomeY1();
                    int gX2 = pixel.getGenomeX2();
                    int gY2 = pixel.getGenomeY2();

                    int prevSize = 0;
                    Set<ContactRecordBox> pixelList = new HashSet<>();
                    pixelList.add(pixel);

                    while (prevSize != pixelList.size()) {
                        prevSize = pixelList.size();
                        for (ContactRecordBox px : records) {
                            if (hasOverlapWith(px, gX1, gY1, gX2, gY2)) {
                                pixelList.add(px);
                                gX1 = Math.min(gX1, px.getGenomeX1());
                                gY1 = Math.min(gY1, px.getGenomeY1());
                                gX2 = Math.max(gX2, px.getGenomeX2());
                                gY2 = Math.max(gY2, px.getGenomeY2());
                            }
                        }
                        records.removeAll(pixelList);
                    }

                    assessNeighbors(coalesced, pixel, pixelList, minResolutions);
                }
            }
        }
        return coalesced;
    }

    private static void sortInPlace(LinkedList<ContactRecordBox> records) {
        records.sort((o1, o2) -> {
            if (o1.getResolution() == o2.getResolution()) {
                return -Float.compare(o1.getCounts(), o2.getCounts());
            } else {
                return Integer.compare(o1.getResolution(), o2.getResolution());
            }
        });
    }

    private static Map<SimpleLocation, LinkedList<ContactRecordBox>> groupNearbyRecords(Set<ContactRecordBox> initialPoints, int scalar) {
        Map<SimpleLocation, LinkedList<ContactRecordBox>> locationMap = new HashMap<>();
        for (ContactRecordBox cr : initialPoints) {
            SimpleLocation region = new SimpleLocation(cr.getGenomeX1() / scalar, cr.getGenomeY1() / scalar);
            if (locationMap.containsKey(region)) {
                locationMap.get(region).add(cr);
            } else {
                LinkedList<ContactRecordBox> values = new LinkedList<>();
                values.add(cr);
                locationMap.put(region, values);
            }
        }
        return locationMap;
    }

    private static boolean hasOverlapWith(ContactRecordBox px, int gX1, int gY1, int gX2, int gY2) {
        return getWidth(px.getGenomeX1(), px.getGenomeX2(), gX1, gX2) *
                getWidth(px.getGenomeY1(), px.getGenomeY2(), gY1, gY2) > 0;
    }

    private static int getWidth(int boxG1, int boxG2, int g1, int g2) {
        int start = Math.max(boxG1, g1);
        int end = Math.min(boxG2, g2);
        return Math.max(0, end - start);
    }

    private static void assessNeighbors(Set<ContactRecordBox> coalesced, ContactRecordBox pixel, Set<ContactRecordBox> pixelList, int minResolutions) {
        Set<Integer> resolutions = new HashSet<>();
        for (ContactRecordBox box : pixelList) {
            resolutions.add(box.getResolution());
        }
        if (resolutions.size() >= minResolutions) {
            coalesced.add(pixel);
        }
    }

    private static boolean pixelMoreEnrichedThanNeighbors(ContactRecord pixel, Set<ContactRecord> pixelList) {
        float sumTotal = getTotalCountsSum(pixelList) - pixel.getCounts();
        float average = sumTotal / (pixelList.size() - 1);
        return pixel.getCounts() / average > Sift.ENRICHMENT_VS_NEIGHBORS;
    }

    private static float getTotalCountsSum(Set<ContactRecord> records) {
        float total = 0;
        for (ContactRecord record : records) {
            total += record.getCounts();
        }
        return total;
    }
}
