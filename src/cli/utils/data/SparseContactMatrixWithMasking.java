package cli.utils.data;

import cli.utils.apa.DistanceBoundCalculator;
import javastraw.expected.ExpectedUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class SparseContactMatrixWithMasking extends SparseContactRecordStorage {

    public SparseContactMatrixWithMasking(MatrixZoomData zd, Collection<Feature2D> loops, int resolution,
                                          int window, int matrixWidth, NormalizationType norm, boolean isIntra) {
        Set<Integer> allIndices = getIndices(loops, resolution, window, matrixWidth);
        DistanceBoundCalculator distanceBoundCalculator = DistanceBoundCalculator.fromFlatList(loops, window,
                resolution, isIntra);

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (distanceBoundCalculator.inDistanceRange(cr)) {
                    if (allIndices.contains(cr.getBinX()) && allIndices.contains(cr.getBinY())) {
                        if (!data.containsKey(cr.getBinX())) {
                            data.put(cr.getBinX(), new HashMap<>());
                        }
                        data.get(cr.getBinX()).put(cr.getBinY(), cr);
                    }
                }
            }
        }

        allIndices.clear();
    }

    private static Set<Integer> getIndices(Collection<Feature2D> loops, int resolution, int window, int matrixWidth) {
        Set<Integer> indices = new HashSet<>();
        for (Feature2D loop : loops) {
            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);

            for (int i = 0; i < matrixWidth; i++) {
                indices.add(i + binXStart);
                indices.add(i + binYStart);
            }
        }
        return indices;
    }
}
