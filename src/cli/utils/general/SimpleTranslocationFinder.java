package cli.utils.general;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;

public class SimpleTranslocationFinder {

    private static final int capacity = 10;
    private static final int zCutOff = 5;

    public static TranslocationSet find(Dataset ds, NormalizationType norm, int res) {
        TranslocationSet tSet = new TranslocationSet();
        Chromosome[] chroms = ds.getChromosomeHandler().getAutosomalChromosomesArray();

        Welford globalVals = new Welford();

        Map<Integer, Map<Integer, PriorityQueue<Float>>> interRegionToValues = new HashMap<>();

        for (int i = 0; i < chroms.length; i++) {
            interRegionToValues.put(i, new HashMap<>());
            for (int j = i + 1; j < chroms.length; j++) {
                Matrix matrix = ds.getMatrix(chroms[i], chroms[j]);
                if (matrix != null) {
                    MatrixZoomData zd = matrix.getZoomData(new HiCZoom(res));
                    if (zd != null) {
                        PriorityQueue<Float> minHeap = new PriorityQueue<>(capacity);
                        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
                        while (iterator.hasNext()) {
                            ContactRecord record = iterator.next();
                            if (record.getCounts() > 0) {
                                globalVals.addValue(record.getCounts());
                                if (minHeap.size() < capacity) {
                                    minHeap.offer(record.getCounts());
                                } else if (record.getCounts() > minHeap.peek()) {
                                    minHeap.poll();       // Remove the smallest number
                                    minHeap.offer(record.getCounts()); // Add the new number
                                }
                            }
                        }
                        if (!minHeap.isEmpty()) {
                            interRegionToValues.get(i).put(j, minHeap);
                        }
                    }
                }
                System.out.print(".");
            }
            System.out.println(".");
        }

        Zscore globalZscore = globalVals.getZscore();
        for (int i = 0; i < chroms.length; i++) {
            for (int j = i + 1; j < chroms.length; j++) {
                if (interRegionToValues.get(i).containsKey(j)) {
                    PriorityQueue<Float> minHeap = interRegionToValues.get(i).get(j);
                    if (!minHeap.isEmpty()) {
                        float smallest = minHeap.peek();
                        if (globalZscore.getZscore(smallest) > zCutOff) {
                            tSet.add(chroms[i], chroms[j]);
                        }
                    }
                } else {
                    tSet.add(chroms[i], chroms[j]);
                }
            }
        }

        return tSet;
    }
}
