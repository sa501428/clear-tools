package cli.utils.data;

import cli.utils.stripes.StripeUtils;
import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class StripeFinder {

    private final SparseFilteredOEMap map;
    private final int minLengthStripe;
    private final int minPeakDist;
    private final int maxPeakDist;
    private final Chromosome chrom;
    private final int resolution;

    public StripeFinder(SparseFilteredOEMap map, Chromosome chrom, int resolution, int minPeakDist, int maxPeakDist, int minLengthStripe) {
        this.map = map;
        this.chrom = chrom;
        this.resolution = resolution;
        this.minPeakDist = minPeakDist;
        this.maxPeakDist = maxPeakDist;
        this.minLengthStripe = minLengthStripe;
    }

    private static boolean elevated5(float[] data, int i) {
        return data[i] > data[i - 1]
                && data[i] > data[i + 1]
                && data[i - 1] > data[i - 2]
                && data[i + 1] > data[i + 2];
    }

    public List<Feature2D> getHorizontalStripes() {
        List<Integer> localPeaks = getLocalPeaks(map.getHorizontalSignal());

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
        List<Integer> localPeaks = getLocalPeaks(map.getVerticalSignal());

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

        for (int[] stretch : StripeUtils.onlySignificantStretches(stretches, dataSlice, dataOESlice)) {
            stripes.add(makeHorizontalStripe(i,
                    i + minPeakDist + stretch[0],
                    i + minPeakDist + stretch[1], resolution));
        }
    }

    private void complexVerticalCall(int j, List<Feature2D> stripes) {
        float[][] dataSlice = getVerticalCountSlice(j);
        float[][] dataOESlice = getVerticalOESlice(j);
        List<int[]> stretches = StripeUtils.findContiguousStretches(dataSlice, dataOESlice, minLengthStripe);

        for (int[] stretch : StripeUtils.onlySignificantStretches(stretches, dataSlice, dataOESlice)) {
            stripes.add(makeVerticalStripe(j,
                    j - maxPeakDist + stretch[0],
                    j - maxPeakDist + stretch[1], resolution));
        }
    }

    private float[][] getHorizontalCountSlice(int i0) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int ii = 0; ii < 5; ii++) {
            for (int j = 0; j < maxPeakDist - minPeakDist; j++) {
                dataSlice[ii][j] = map.getCountValue(i0 - 2 + ii, j + i0 + minPeakDist);
            }
        }
        return dataSlice;
    }

    private float[][] getVerticalCountSlice(int j) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int k = 0; k < 5; k++) {
            for (int i = 0; i < maxPeakDist - minPeakDist; i++) {
                dataSlice[k][i] = map.getCountValue(j - maxPeakDist + i, j - 2 + k);
            }
        }
        return dataSlice;
    }

    public float[][] getHorizontalOESlice(int i0) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int ii = 0; ii < 5; ii++) {
            for (int j = 0; j < maxPeakDist - minPeakDist; j++) {
                dataSlice[ii][j] = map.getOEValue(i0 - 2 + ii, j + i0 + minPeakDist);
            }
        }
        return dataSlice;
    }

    public float[][] getVerticalOESlice(int j) {
        float[][] dataSlice = new float[5][maxPeakDist - minPeakDist];
        for (int k = 0; k < 5; k++) {
            for (int i = 0; i < maxPeakDist - minPeakDist; i++) {
                dataSlice[k][i] = map.getOEValue(j - maxPeakDist + i, j - 2 + k);
            }
        }
        return dataSlice;
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
}
