package cli.utils.data;

import cli.utils.general.ArrayTools;
import cli.utils.general.VectorCleaner;
import javastraw.expected.ExpectedUtils;
import javastraw.expected.LogExpectedZscoreSpline;
import javastraw.feature2D.Feature2D;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.awt.*;
import java.util.List;
import java.util.*;

public class SparseFilteredOEMap {

    private static final float Z_TOP_TEN = 1.28f;
    private final Map<Integer, Map<Integer, Float>> contactMap = new HashMap<>();
    private final int minPeakDist;
    private final int maxPeakDist;
    private final float[] upStreamLogSignal;
    private final float[] downStreamLogSignal;
    private final Chromosome chrom;
    private final int minLengthStripe;
    private final int resolution;

    public SparseFilteredOEMap(Dataset ds, MatrixZoomData zd, NormalizationType norm, Chromosome chrom, int resolution,
                               int minPeakDist, int maxPeakDist, HiCZoom zoom, int minLengthStripe) {
        this.chrom = chrom;
        this.minPeakDist = minPeakDist;
        this.maxPeakDist = maxPeakDist;
        this.minLengthStripe = minLengthStripe;
        this.resolution = resolution;

        int numEntries = (int) ((chrom.getLength() / resolution) + 1);

        upStreamLogSignal = new float[numEntries];
        downStreamLogSignal = new float[numEntries];
        //Arrays.fill(upStreamSignal, 1);
        //Arrays.fill(downStreamSignal, 1);

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
                        if (oe > 2 || zscore > Z_TOP_TEN) {
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

    private static boolean elevated5(float[] data, int i) {
        return data[i] > 1.1 * data[i - 1]
                && data[i] > 1.1 * data[i + 1]
                && data[i - 1] > 1.1 * data[i - 2]
                && data[i + 1] > 1.1 * data[i + 2]
                && data[i - 2] > 0
                && data[i + 2] > 0;
        //&& data[i + 1] > 1
        //&& data[i - 1] > 1);
    }

    private void populateMap(ContactRecord cr, float oe) {
        if (!contactMap.containsKey(cr.getBinX())) {
            contactMap.put(cr.getBinX(), new HashMap<>());
        }
        contactMap.get(cr.getBinX()).put(cr.getBinY(), oe);
    }

    private float getValue(int binX, int binY) {
        if (contactMap.containsKey(binX) && contactMap.get(binX).containsKey(binY)) {
            return contactMap.get(binX).get(binY);
        }
        return 1;
    }

    public void clear() {
        for (Map<Integer, Float> map : contactMap.values()) {
            map.clear();
        }
        contactMap.clear();
    }

    public List<Feature2D> getHorizontalStripes() {
        List<Integer> localPeaks = getLocalPeaks(upStreamLogSignal);
        List<Feature2D> stripes = new LinkedList<>();

        for (int i : localPeaks) {
            int stripeStart = -1;
            int stripeEnd = -1;
            int stripeLen = 0;
            for (int j = i + minPeakDist; j < i + maxPeakDist; j++) {
                if (enrichedOverHorizontalWindow(i, j)) {
                    if (stripeStart == -1) {
                        stripeStart = j;
                    }
                    stripeEnd = j;
                    stripeLen++;
                } else {
                    if (stripeLen >= minLengthStripe) {
                        stripes.add(makeHorizontalStripe(i, stripeStart, stripeEnd, resolution));
                    }
                    stripeStart = -1;
                    stripeEnd = -1;
                    stripeLen = 0;
                }
            }
            if (stripeLen >= minLengthStripe) {
                stripes.add(makeHorizontalStripe(i, stripeStart, stripeEnd, resolution));
            }
        }
        return stripes;
    }

    private boolean enrichedOverHorizontalWindow(int i, int j) {
        float currRowOE = getValue(i, j - 1) + getValue(i, j) + getValue(i, j + 1);
        float prevRowOE = getValue(i - 1, j - 1) + getValue(i - 1, j) + getValue(i - 1, j + 1);
        float nextRowOE = getValue(i + 1, j - 1) + getValue(i + 1, j) + getValue(i + 1, j + 1);
        return currRowOE > 1.25 * prevRowOE && currRowOE > 1.25 * nextRowOE;
    }

    private boolean enrichedOverVerticalWindow(int i, int j) {
        float currColOE = getValue(i - 1, j) + getValue(i, j) + getValue(i + 1, j);
        float prevColOE = getValue(i - 1, j - 1) + getValue(i, j - 1) + getValue(i + 1, j - 1);
        float nextColOE = getValue(i - 1, j + 1) + getValue(i, j + 1) + getValue(i + 1, j + 1);
        return currColOE > 1.25 * prevColOE && currColOE > 1.25 * nextColOE;
    }

    private Feature2D makeHorizontalStripe(int i, int stripeStart, int stripeEnd, int resolution) {
        return new Feature2D(Feature2D.FeatureType.NONE,
                chrom.getName(), (long) i * resolution, (long) (i + 1) * resolution,
                chrom.getName(), (long) stripeStart * resolution, (long) stripeEnd * resolution,
                Color.MAGENTA, new HashMap<>());
    }

    public List<Feature2D> getVerticalStripes() {
        List<Integer> localPeaks = getLocalPeaks(downStreamLogSignal);
        List<Feature2D> stripes = new LinkedList<>();

        for (int j : localPeaks) {
            int stripeStart = -1;
            int stripeEnd = -1;
            int stripeLen = 0;
            for (int i = j - maxPeakDist; i < j - minPeakDist; i++) {
                if (enrichedOverVerticalWindow(i, j)) {
                    if (stripeStart == -1) {
                        stripeStart = i;
                    }
                    stripeEnd = i;
                    stripeLen++;
                } else {
                    if (stripeLen >= minLengthStripe) {
                        stripes.add(makeVerticalStripe(j, stripeStart, stripeEnd, resolution));
                    }
                    stripeStart = -1;
                    stripeEnd = -1;
                    stripeLen = 0;
                }
            }
            if (stripeLen >= minLengthStripe) {
                stripes.add(makeVerticalStripe(j, stripeStart, stripeEnd, resolution));
            }
        }
        return stripes;
    }

    private Feature2D makeVerticalStripe(int j, int stripeStart, int stripeEnd, int resolution) {
        return new Feature2D(Feature2D.FeatureType.NONE,
                chrom.getName(), (long) stripeStart * resolution, (long) stripeEnd * resolution,
                chrom.getName(), (long) j * resolution, (long) (j + 1) * resolution,
                Color.MAGENTA, new HashMap<>());
    }

    private List<Integer> getLocalPeaks(float[] data) {
        List<Integer> peaks = new LinkedList<>();
        for (int i = 1; i < data.length - 1; i++) {
            if (elevated5(data, i)) {
                peaks.add(i);
            }
        }
        return peaks;
    }

    public float[] getHorizontalSignal() {
        return upStreamLogSignal;
    }

    public float[] getVerticalSignal() {
        return downStreamLogSignal;
    }
}
