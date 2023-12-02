package cli.utils.data;

import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class SparseContactMatrixOfSpecificRegionsOnly extends SparseContactRecordStorage {

    public SparseContactMatrixOfSpecificRegionsOnly(MatrixZoomData zd, List<Feature2D> features, int resolution,
                                                    int genomeBuffer, NormalizationType norm) {
        Map<Integer, Set<Integer>> exactMask = new HashMap<>();
        populateAllBinsToTrack(features, resolution, genomeBuffer / resolution, exactMask);

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (exactMask.containsKey(cr.getBinX()) && exactMask.get(cr.getBinX()).contains(cr.getBinY())) {
                    if (!data.containsKey(cr.getBinX())) {
                        data.put(cr.getBinX(), new HashMap<>());
                    }
                    data.get(cr.getBinX()).put(cr.getBinY(), cr);
                }
            }
        }
        exactMask.clear();
    }


    private static void populateAllBinsToTrack(List<Feature2D> stripes,
                                               int resolution, int buffer, Map<Integer, Set<Integer>> exactMask) {
        for (Feature2D loop : stripes) {
            for (int gx = (int) (loop.getStart1() / resolution) - buffer;
                 gx <= (int) (loop.getEnd1() / resolution) + buffer; gx++) {

                if (!exactMask.containsKey(gx)) {
                    exactMask.put(gx, new HashSet<>());
                }

                for (int gy = (int) (loop.getStart2() / resolution) - buffer;
                     gy <= (int) (loop.getEnd2() / resolution) + buffer; gy++) {

                    exactMask.get(gx).add(gy);
                }
            }
        }
    }
}
