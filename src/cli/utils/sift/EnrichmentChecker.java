package cli.utils.sift;

import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class EnrichmentChecker {


    public static void filterOutIfNotLocalMax(MatrixZoomData zdLow, Set<ContactRecord> initialPoints, int scalar) {

        List<BoundingBoxWithContacts> boxes = getBoundingBoxes(NMSUtils.getLocationMap(initialPoints, scalar * 200));


    }

    private static List<BoundingBoxWithContacts> getBoundingBoxes(Map<SimpleLocation, List<ContactRecord>> locationMap) {
        List<BoundingBoxWithContacts> boxes = new ArrayList<>();
        for (List<ContactRecord> contacts : locationMap.values()) {
            boxes.add(new BoundingBoxWithContacts(contacts));
        }
        return boxes;
    }
}
