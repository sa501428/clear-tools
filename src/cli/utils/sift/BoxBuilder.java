package cli.utils.sift;

import javastraw.reader.block.ContactRecord;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class BoxBuilder {


    /*
    public static void filterOutIfNotLocalMax(MatrixZoomData zdLow, Set<ContactRecord> initialPoints, NormalizationType norm) {

        List<BoundingBoxWithContacts> boxes = getBoundingBoxes(
                NMSUtils.groupNearbyRecords(initialPoints, 250));

        for (BoundingBoxWithContacts box : boxes) {
            Set<ContactRecord> toRemove = box.findPointsNotEnriched(zdLow, norm);
            initialPoints.removeAll(toRemove);
        }
    }
    */

    static List<BoundingBoxWithContacts> getBoundingBoxes(Map<SimpleLocation, List<ContactRecord>> locationMap) {
        List<BoundingBoxWithContacts> boxes = new ArrayList<>();
        for (List<ContactRecord> contacts : locationMap.values()) {
            boxes.add(new BoundingBoxWithContacts(contacts));
        }
        return boxes;
    }
}
