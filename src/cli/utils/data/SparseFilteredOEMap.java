package cli.utils.data;

import cli.utils.general.ArrayTools;
import cli.utils.general.VectorCleaner;
import javastraw.expected.ExpectedUtils;
import javastraw.expected.LogExpectedZscoreSpline;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class SparseFilteredOEMap {

    private static final float Z_TOP_TEN = 1.28f;
    private static final float Z_TOP_FIFTEEN = 1.04f;
    private static final float Z_1 = 1f;

    private final Map<Integer, Map<Integer, FloatPair>> contactMap = new HashMap<>();
    private final float[] upStreamLogSignal;
    private final float[] downStreamLogSignal;


    public SparseFilteredOEMap(Dataset ds, MatrixZoomData zd, NormalizationType norm, Chromosome chrom, int resolution,
                               int minPeakDist, int maxPeakDist, HiCZoom zoom) {

        int numEntries = (int) ((chrom.getLength() / resolution) + 1);

        upStreamLogSignal = new float[numEntries];
        downStreamLogSignal = new float[numEntries];

        LogExpectedZscoreSpline poly = new LogExpectedZscoreSpline(zd, norm, chrom, resolution);

        double[] vector = ArrayTools.copy(ds.getNormalizationVector(chrom.getIndex(), zoom,
                NormalizationHandler.VC).getData().getValues().get(0));
        VectorCleaner.inPlaceZscore(vector);

        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 0) {
                if (vector[cr.getBinX()] > -1 && vector[cr.getBinY()] > -1 &&
                        vector[cr.getBinX() - 1] > -1 && vector[cr.getBinY() - 1] > -1 &&
                        vector[cr.getBinX() + 1] > -1 && vector[cr.getBinY() + 1] > -1) {

                    int dist = ExpectedUtils.getDist(cr);
                    if (dist > minPeakDist && dist < maxPeakDist) {

                        float oe = (float) ((cr.getCounts() + 1) / (poly.getExpectedFromUncompressedBin(dist) + 1));
                        float zscore = (float) poly.getZscoreForObservedUncompressedBin(dist, cr.getCounts());
                        if (oe > 1.5 || zscore > Z_1) { // oe > 2 || zscore > Z_TOP_TEN
                            populateMap(cr, cr.getCounts());
                            //upStreamSignal[cr.getBinX()] *= oe;
                            //downStreamSignal[cr.getBinY()] *= oe;
                            upStreamLogSignal[cr.getBinX()] += Math.log(oe);
                            downStreamLogSignal[cr.getBinY()] += Math.log(oe);
                        }
                    }
                }
            }
        }
    }


    public void populateMap(ContactRecord cr, float oe) {
        if (!contactMap.containsKey(cr.getBinX())) {
            contactMap.put(cr.getBinX(), new HashMap<>());
        }
        contactMap.get(cr.getBinX()).put(cr.getBinY(), new FloatPair(cr.getCounts(), oe));
    }

    public float getOEValue(int binX, int binY) {
        if (contactMap.containsKey(binX) && contactMap.get(binX).containsKey(binY)) {
            return contactMap.get(binX).get(binY).oe;
        }
        return 1;
    }

    public float getCountValue(int binX, int binY) {
        if (contactMap.containsKey(binX) && contactMap.get(binX).containsKey(binY)) {
            return contactMap.get(binX).get(binY).count;
        }
        return 0;
    }

    public void clear() {
        for (Map<Integer, FloatPair> map : contactMap.values()) {
            map.clear();
        }
        contactMap.clear();
    }

    public float[] getHorizontalSignal() {
        return upStreamLogSignal;
    }

    public float[] getVerticalSignal() {
        return downStreamLogSignal;
    }
}
