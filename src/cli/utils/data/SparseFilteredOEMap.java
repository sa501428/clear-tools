package cli.utils.data;

import cli.utils.general.ArrayTools;
import cli.utils.general.VectorCleaner;
import cli.utils.stripes.StripeUtils;
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
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.util.List;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class SparseFilteredOEMap {

    private static final float Z_TOP_TEN = 1.28f;
    private static final float Z_TOP_FIFTEEN = 1.04f;
    private static final float Z_1 = 1f;

    private final Map<Integer, Map<Integer, FloatPair>> contactMap = new HashMap<>();
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

    private static boolean elevated5original(float[] data, int i) {
        return data[i] > 1.1 * data[i - 1]
                && data[i] > 1.1 * data[i + 1]
                && data[i - 1] > 1.05 * data[i - 2]
                && data[i + 1] > 1.05 * data[i + 2]
                && data[i - 2] > 0
                && data[i + 2] > 0;
        //&& data[i + 1] > 1
        //&& data[i - 1] > 1);
    }

    private static boolean elevated5(float[] data, int i) {
        return data[i] > data[i - 1]
                && data[i] > data[i + 1]
                && data[i - 1] > data[i - 2]
                && data[i + 1] > data[i + 2];
    }

    private void populateMap(ContactRecord cr, float oe) {
        if (!contactMap.containsKey(cr.getBinX())) {
            contactMap.put(cr.getBinX(), new HashMap<>());
        }
        contactMap.get(cr.getBinX()).put(cr.getBinY(), new FloatPair(cr.getCounts(), oe));
    }

    private float getOEValue(int binX, int binY) {
        if (contactMap.containsKey(binX) && contactMap.get(binX).containsKey(binY)) {
            return contactMap.get(binX).get(binY).oe;
        }
        return 1;
    }

    private float getCountValue(int binX, int binY) {
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





    private Feature2D makeHorizontalStripe(int i, int stripeStart, int stripeEnd, int resolution) {
        return new Feature2D(Feature2D.FeatureType.NONE,
                chrom.getName(), (long) i * resolution, (long) (i + 1) * resolution,
                chrom.getName(), (long) stripeStart * resolution, (long) stripeEnd * resolution,
                Color.MAGENTA, new HashMap<>());
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

    public List<Feature2D> getHorizontalStripes() {
        List<Integer> localPeaks = getLocalPeaks(upStreamLogSignal);

        List<Feature2D> globalStripes = new LinkedList<>();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            List<Feature2D> localStripes = new LinkedList<>();
            int i = index.getAndIncrement();
            while (i < localPeaks.size()) {
                //simpleHorizontalCall(localPeaks.get(i), localStripes);
                complexHorizontalCall(localPeaks.get(i), localStripes);
                i = index.getAndIncrement();
            }
            synchronized (globalStripes) {
                globalStripes.addAll(localStripes);
            }
        });
        return globalStripes;
    }

    public List<Feature2D> getVerticalStripes() {
        List<Integer> localPeaks = getLocalPeaks(downStreamLogSignal);

        List<Feature2D> globalStripes = new LinkedList<>();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            List<Feature2D> localStripes = new LinkedList<>();
            int i = index.getAndIncrement();
            while (i < localPeaks.size()) {
                //simpleVerticalCall(localPeaks.get(i), localStripes);
                complexVerticalCall(localPeaks.get(i), localStripes);
                i = index.getAndIncrement();
            }
            synchronized (globalStripes) {
                globalStripes.addAll(localStripes);
            }
        });
        return globalStripes;
    }

    private void complexHorizontalCall(int i, List<Feature2D> stripes) {
        float[][] dataSlice = getHorizontalCountSlice(i);
        float[][] dataOESlice = getHorizontalOESlice(i);
        List<int[]> stretches = StripeUtils.findContiguousStretches(dataSlice, dataOESlice, minLengthStripe);
        for (int[] stretch : stretches) {
            stripes.add(makeHorizontalStripe(i,
                    i + minPeakDist + stretch[0],
                    i + minPeakDist + stretch[1], resolution));
        }
    }

    private void complexVerticalCall(int j, List<Feature2D> stripes) {
        float[][] dataCountSlice = getVerticalCountSlice(j);
        float[][] dataOESlice = getVerticalOESlice(j);
        List<int[]> stretches = StripeUtils.findContiguousStretches(dataCountSlice, dataOESlice, minLengthStripe);
        for (int[] stretch : stretches) {
            stripes.add(makeVerticalStripe(j,
                    j - maxPeakDist + stretch[0],
                    j - maxPeakDist + stretch[1], resolution));
        }
    }

    private float[][] getHorizontalCountSlice(int i0) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int ii = 0; ii < 5; ii++) {
            for (int j = 0; j < maxPeakDist - minPeakDist; j++) {
                dataSlice[ii][j] = getCountValue(i0 - 2 + ii, j + i0 + minPeakDist);
            }
        }
        return dataSlice;
    }

    private float[][] getVerticalCountSlice(int j) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int k = 0; k < 5; k++) {
            for (int i = 0; i < maxPeakDist - minPeakDist; i++) {
                dataSlice[k][i] = getCountValue(j - maxPeakDist + i, j - 2 + k);
            }
        }
        return dataSlice;
    }

    private float[][] getHorizontalOESlice(int i0) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int ii = 0; ii < 5; ii++) {
            for (int j = 0; j < maxPeakDist - minPeakDist; j++) {
                dataSlice[ii][j] = getOEValue(i0 - 2 + ii, j + i0 + minPeakDist);
            }
        }
        return dataSlice;
    }

    private float[][] getVerticalOESlice(int j) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int k = 0; k < 5; k++) {
            for (int i = 0; i < maxPeakDist - minPeakDist; i++) {
                dataSlice[k][i] = getOEValue(j - maxPeakDist + i, j - 2 + k);
            }
        }
        return dataSlice;
    }

}
