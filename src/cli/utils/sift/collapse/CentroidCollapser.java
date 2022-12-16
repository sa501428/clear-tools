package cli.utils.sift.collapse;

import cli.clt.loops.Sift;
import cli.utils.sift.NMSUtils;
import cli.utils.sift.SimpleLocation;
import javastraw.reader.block.ContactRecord;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

public class CentroidCollapser {

    public static Set<ContactRecord> coalesce(Set<ContactRecord> regions, int buffer, int minNeighbors) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        Map<SimpleLocation, LinkedList<ContactRecord>> map = NMSUtils.groupNearbyRecords(regions, 250);
        Set<ContactRecord> coalesced = new HashSet<>();

        for (LinkedList<ContactRecord> records : map.values()) {
            sortInPlace(records);

            while (!records.isEmpty()) {
                ContactRecord pixel = records.pollFirst();
                if (pixel != null) {
                    records.remove(pixel);

                    int binX0 = pixel.getBinX() - buffer;
                    int binY0 = pixel.getBinY() - buffer;
                    int binX1 = pixel.getBinX() + buffer + 1;
                    int binY1 = pixel.getBinY() + buffer + 1;

                    int prevSize = 0;
                    Set<ContactRecord> pixelList = new HashSet<>();
                    pixelList.add(pixel);

                    while (prevSize != pixelList.size()) {
                        prevSize = pixelList.size();
                        for (ContactRecord px : records) {
                            if (contains(px, binX0, binY0, binX1, binY1)) {
                                pixelList.add(px);
                                binX0 = Math.min(binX0, px.getBinX() - buffer);
                                binY0 = Math.min(binY0, px.getBinY() - buffer);
                                binX1 = Math.max(binX1, px.getBinX() + buffer + 1);
                                binY1 = Math.max(binY1, px.getBinY() + buffer + 1);
                            }
                        }
                        records.removeAll(pixelList);
                    }

                    assessNeighbors(coalesced, pixel, pixelList, minNeighbors);
                }
            }
        }
        return coalesced;
    }

    private static void sortInPlace(LinkedList<ContactRecord> records) {
        records.sort((o1, o2) -> Float.compare(-o1.getCounts(), -o2.getCounts()));
    }

    private static boolean contains(ContactRecord px, int binX0, int binY0, int binX1, int binY1) {
        return binX0 <= px.getBinX() && binX1 > px.getBinX() && binY0 <= px.getBinY() && binY1 > px.getBinY();
    }

    private static void assessNeighbors(Set<ContactRecord> coalesced, ContactRecord pixel, Set<ContactRecord> pixelList, int minNeighbors) {
        if (pixelList.size() > minNeighbors) {
            if (pixelMoreEnrichedThanNeighbors(pixel, pixelList)) {
                coalesced.add(pixel);
            }
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
