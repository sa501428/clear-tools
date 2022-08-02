package cli.clt;

import cli.utils.ExpectedUtils;
import cli.utils.Welford;
import cli.utils.expected.LogExpectedModel;
import cli.utils.sift.SimpleLocation;
import cli.utils.sift.Zscore;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.io.File;
import java.util.List;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class HotSpot {
    private static final int MAX_DIST = 10000000;
    private static final int MIN_DIST = 100000;

    // The Z-Score Cutoff is currently hardcoded at 2, which has a confidence interval of 97.72%
    private static final float ZSCORE_CUTOFF = 2f;

    private static void printUsageAndExit() {
        System.out.println("hotspot [--res resolution] [--norm normalization] <file1.hic,file2.hic,...> <out_file>");
        System.exit(19);
    }

    public static void run(String[] args, CommandLineParser parser) {

        // hotspot [--res int] [--norm string] <file1.hic,file2.hic,...> <out_file>
        if (args.length != 3) {
            printUsageAndExit();
        }

        String[] files = args[1].split(",");
        String outputFileName = args[2];
        int resolution = parser.getResolutionOption(5000);
        String normString = parser.getNormalizationStringOption();

        final Feature2DList result = new Feature2DList();

        final int countThreshold;
        if (files.length > 3) {
            countThreshold = 3;
        } else {
            countThreshold = 2;
        }

        Dataset[] datasets = new Dataset[files.length];
        for (int k = 0; k < files.length; k++) {
            datasets[k] = HiCFileTools.extractDatasetForCLT(files[k],
                    false, false, false);
        }

        NormalizationType norm = datasets[0].getNormalizationHandler().getNormTypeFromString(normString);
        Chromosome[] chromosomes = datasets[0].getChromosomeHandler().getAutosomalChromosomesArray();
        AtomicInteger cIndex = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chrom = chromosomes[currIndex];
                // for iterating on chr17 alone
                // if (chrom.getIndex() == 17) {...
                List<Feature2D> hotspots = findTheHotspots(chrom, datasets, resolution, norm, countThreshold);
                if (hotspots.size() > 0) {
                    synchronized (result) {
                        result.addByKey(Feature2DList.getKey(chrom, chrom), hotspots);
                    }
                }
                // ... }
                currIndex = cIndex.getAndIncrement();
            }
        });

        result.exportFeatureList(new File(outputFileName + ".hotspot.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("hotspot complete");
    }

    private static List<Feature2D> findTheHotspots(Chromosome chrom, Dataset[] datasets, int resolution,
                                                   NormalizationType norm, int countThreshold) {
        Map<SimpleLocation, Welford> results = new HashMap<>();
        int minBin = MIN_DIST / resolution;
        int maxBin = MAX_DIST / resolution;

        for (Dataset ds : datasets) {
            Matrix matrix;
            synchronized (ds) {
                matrix = ds.getMatrix(chrom, chrom, resolution);
            }
            if (matrix == null) {
                System.err.println("matrix for " + chrom.getName() + " == null -> continuing");
                continue;
            }
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) {
                System.err.println("zd for " + chrom.getName() + " == null -> continuing");
                continue;
            }

            // iterating through chrom using type 1 iteration
            NormalizationType scaleNorm = ds.getNormalizationHandler().getNormTypeFromString("SCALE");
            NormalizationType vcNorm = ds.getNormalizationHandler().getNormTypeFromString("VC");

            double[] vector1 = ds.getNormalizationVector(chrom.getIndex(), new HiCZoom(resolution), scaleNorm).getData().getValues().get(0);
            double[] vector2 = ds.getNormalizationVector(chrom.getIndex(), new HiCZoom(resolution), vcNorm).getData().getValues().get(0);

            iterateThruAndGrabPercentContact(zd, maxBin, minBin, norm, results, vector1, vector2);
            matrix.clearCache();
        }

        List<SimpleLocation> removeList = new ArrayList<>();
        for (SimpleLocation key : results.keySet()) {
            double range = results.get(key).getRange();
            int counts = (int) results.get(key).getCounts();
            double min = results.get(key).getMin();
            double max = results.get(key).getMax();
            //value.addZeroIfBelow(files.length);

            if (range < .05 || range > 0.35 || counts < countThreshold || min > .05 || max > 0.35) {
                removeList.add(key);
            }
            //if (entry.getValue().getCounts() < NUM_NONZERO_VALUES_THRESHOLD)
            //    removeList.add(entry.getKey());
        }
        int n1 = results.size();
        int n2 = removeList.size();
        // test print


        for (SimpleLocation key : removeList) {
            results.remove(key);
        }
        // test print

        int n3 = results.size();

        removeList.clear();


        List<Feature2D> hotspots = new ArrayList<>();
        if (results.values().size() > 1) {

            // comment out below in order to use Z-Score Cutoff
            for (Map.Entry<SimpleLocation, Welford> entry : results.entrySet()) {

                Welford welford = entry.getValue();
                Map<String, String> attributes = new HashMap<>();
                attributes.put("sigma", "" + welford.getStdDev());
                attributes.put("range", "" + welford.getRange());
                attributes.put("min", "" + welford.getMin());
                attributes.put("max", "" + welford.getMax());
                long startX = (long) entry.getKey().getBinX() * resolution;
                long endX = startX + resolution;
                long startY = (long) entry.getKey().getBinY() * resolution;
                long endY = startY + resolution;
                Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, chrom.getName(), startX, endX, chrom.getName(), startY, endY, Color.BLACK, attributes);
                hotspots.add(feature);

            }
            // comment out ^

//            // uncomment below in order to use Z-Score Cutoff
//            Zscore zscore = getOverallZscoreMetric(results.values());
//            for (Map.Entry<SimpleLocation, Welford> entry : results.entrySet()) {
//                Welford welford = entry.getValue();
//                if (zscore.getZscore(welford.getStdDev()) >= ZSCORE_CUTOFF) {
//                    //if (zscore.getZscore(welford.getRange()) >= ZSCORE_CUTOFF) {
//                    Map<String, String> attributes = new HashMap<>();
//                    attributes.put("sigma", "" + welford.getStdDev());
//                    attributes.put("range", "" + welford.getRange());
//                    attributes.put("min", "" + welford.getMin());
//                    attributes.put("max", "" + welford.getMax());
//                    long startX = (long) entry.getKey().getBinX() * resolution;
//                    long endX = startX + resolution;
//                    long startY = (long) entry.getKey().getBinY() * resolution;
//                    long endY = startY + resolution;
//                    Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, chrom.getName(), startX, endX, chrom.getName(), startY, endY, Color.BLACK, attributes);
//                    hotspots.add(feature);
//                }
//            }
//            // uncomment ^

        }

        int n4 = hotspots.size();
        System.out.println("Hotspots " + chrom.getName() + " pre-filter: " + n1 +
                " removed: " + n2 + " post-filter: " + n3 + " final: " + n4);
        //return hotspots;
        return coalesceAndRetainCentroids(new HashSet<>(hotspots), 3 * resolution);
    }

    private static Zscore getOverallZscoreMetric(Collection<Welford> values) {
        Welford overallWelford = new Welford();
        for (Welford welford : values) {
            overallWelford.addValue(welford.getStdDev());
            //overallWelford.addValue(welford.getRange());
        }
        System.out.println("Overall welford: " + overallWelford.getSummary());
        return overallWelford.getZscore();
    }

    private static void iterateThruAndGrabPercentContact(MatrixZoomData zd, int maxBin, int minBin,
                                                         NormalizationType norm,
                                                         Map<SimpleLocation, Welford> results, double[] vector1, double[] vector2) {

        LogExpectedModel expected = new LogExpectedModel(zd, norm, maxBin, false, 0);

        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            if (cr.getCounts() > 0) {
                int dist = ExpectedUtils.getDist(cr);
                if (vector1[cr.getBinX()] > 1 && vector1[cr.getBinY()] > 1 && vector2[cr.getBinX()] > 1 && vector2[cr.getBinY()] > 1) {
                    if (dist > minBin && dist < maxBin) {
                        float percentContact = expected.getPercentContact(dist, cr.getCounts());
                        percentContact = Math.min(1, Math.max(0, percentContact));
                        //percentContact2 = Math.exp(percentContact0 - 1);

                        SimpleLocation location = new SimpleLocation(cr);
                        results.putIfAbsent(location, new Welford());
                        results.get(location).addValue(percentContact);
                    }
                }
            }
        }
        System.out.print(".");
        //System.out.println("Finished recording locations for this file");
    }

    public static List<Feature2D> coalesceAndRetainCentroids(Set<Feature2D> features, int gRadius) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        LinkedList<Feature2D> featureLL = new LinkedList<>(features);
        List<Feature2D> coalesced = new ArrayList<>();

        //possible alternative: compare by max
        featureLL.sort((o1, o2) -> Double.compare(Double.parseDouble(o1.getAttribute("max")), Double.parseDouble(o2.getAttribute("max"))));
        Collections.reverse(featureLL);


        while (!featureLL.isEmpty()) {
            Feature2D pixel = featureLL.pollFirst();
            if (pixel != null) {
                coalesced.add(pixel);
                featureLL.remove(pixel);
                long x = pixel.getMidPt1();
                long y = pixel.getMidPt2();

                Set<Feature2D> toRemove = new HashSet<>();
                for (Feature2D px : featureLL) {
                    if (distance(x - px.getMidPt1(),
                            y - px.getMidPt2()) <= gRadius) {
                        toRemove.add(px);
                    }
                }
                featureLL.removeAll(toRemove);
            }
        }
        features.clear();
        return coalesced;
    }

    public static long distance(long x, long y) {
        return Math.max(Math.abs(x), Math.abs(y));
    }
}
