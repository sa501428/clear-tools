package cli.clt.loops;

import cli.clt.CommandLineParser;
import cli.utils.hotspot.HotSpotUtils;
import cli.utils.sift.SimpleLocation;
import javastraw.expected.ExpectedModel;
import javastraw.expected.ExpectedUtils;
import javastraw.expected.LogExpectedSpline;
import javastraw.expected.Welford;
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
    public static String usage = "hotspot [--res int] [--norm string] <file1.hic,file2.hic,...> <outfile>";

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
                List<Feature2D> hotspots = findTheHotspots(chrom, datasets, resolution, norm, countThreshold);
                if (hotspots.size() > 0) {
                    synchronized (result) {
                        result.addByKey(Feature2DList.getKey(chrom, chrom), hotspots);
                    }
                }
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
            Matrix matrix = ds.getMatrix(chrom, chrom, resolution);
            if (matrix != null) {
                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                if (zd != null) {
                    // iterating through chrom using type 1 iteration
                    NormalizationType scaleNorm = ds.getNormalizationHandler().getNormTypeFromString("SCALE");
                    NormalizationType vcNorm = ds.getNormalizationHandler().getNormTypeFromString("VC");

                    double[] vector1 = ds.getNormalizationVector(chrom.getIndex(), new HiCZoom(resolution), scaleNorm).getData().getValues().get(0);
                    double[] vector2 = ds.getNormalizationVector(chrom.getIndex(), new HiCZoom(resolution), vcNorm).getData().getValues().get(0);

                    iterateThruAndGrabPercentContact(zd, maxBin, minBin, norm, results, vector1, vector2, chrom, resolution);
                }
                matrix.clearCache();
            }
        }

        initialHighIntensityFilter(results, countThreshold);

        Set<ContactRecord> candidateHotSpotsSet = secondPassLowRangeFilter(results);
        Set<ContactRecord> ubiquitousPeaksSet = findUbiquitousLoops(results);

        Set<ContactRecord> records = new HashSet<>(candidateHotSpotsSet);
        records.addAll(ubiquitousPeaksSet);
        // records is currently the union of the candidateHotSpotsSet and ubiquitousPeaksSet
        HotSpotUtils.coalesceAndRetainCentroids(records, 30000 / resolution);
        records.removeAll(ubiquitousPeaksSet);
        // now records is the final set of hotspots (as ContactRecords)

        List<Feature2D> hotspots = new ArrayList<>();
        if (records.size() > 1) {
            for (ContactRecord record : records) {
                SimpleLocation location = new SimpleLocation(record.getBinX(), record.getBinY());
                Welford welford = results.get(location);
                Map<String, String> attributes = getStats(welford);
                long startX = (long) location.getBinX() * resolution;
                long endX = startX + resolution;
                long startY = (long) location.getBinY() * resolution;
                long endY = startY + resolution;
                Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, chrom.getName(), startX, endX, chrom.getName(), startY, endY, Color.BLACK, attributes);
                hotspots.add(feature);
            }
        }
        return hotspots;
    }

    private static Set<ContactRecord> findUbiquitousLoops(Map<SimpleLocation, Welford> results) {
        Set<ContactRecord> ubiquitousLoopsSet = new HashSet<>();
        for (SimpleLocation key : results.keySet()) {
            double range = results.get(key).getRange();
            double min = results.get(key).getMin();
            if (range < .05 && min > .05) {
                ubiquitousLoopsSet.add(new ContactRecord(key.getBinX(), key.getBinY(),
                        (float) results.get(key).getMax()));
            }
        }
        return ubiquitousLoopsSet;
    }

    private static Set<ContactRecord> secondPassLowRangeFilter(Map<SimpleLocation, Welford> results) {
        Set<ContactRecord> candidateRecordsSet = new HashSet<>();
        for (SimpleLocation key : results.keySet()) {
            double range = results.get(key).getRange();
            double min = results.get(key).getMin();
            if (range > .05 && min < .05) {
                candidateRecordsSet.add(new ContactRecord(key.getBinX(), key.getBinY(),
                        (float) results.get(key).getMax()));
            }

            /*
            if (range < .05 || min > .05) {
                results.remove(key);
                removeList.add(key);
            }
            */
        }
        return candidateRecordsSet;
    }

    private static void initialHighIntensityFilter(Map<SimpleLocation, Welford> results, int countThreshold) {
        List<SimpleLocation> removeList = new ArrayList<>();
        for (SimpleLocation key : results.keySet()) {
            double range = results.get(key).getRange();
            int counts = (int) results.get(key).getCounts();
            double max = results.get(key).getMax();
            if (range > 0.35 || counts < countThreshold || max > 0.35) {
                removeList.add(key);
            }
        }
        removeAllAndClear(results, removeList);
    }

    private static Set<ContactRecord> convert(Set<SimpleLocation> locations, Map<SimpleLocation, Welford> results) {
        Set<ContactRecord> records = new HashSet<>();
        for (SimpleLocation location : locations) {
            records.add(new ContactRecord(location.getBinX(), location.getBinY(),
                    (float) results.get(location).getMax()));
        }
        return records;
    }

    private static void removeAllAndClear(Map<SimpleLocation, Welford> results, List<SimpleLocation> removeList) {
        for (SimpleLocation key : removeList) {
            results.remove(key);
        }
        removeList.clear();
    }

    private static Map<String, String> getStats(Welford welford) {
        Map<String, String> attributes = new HashMap<>();
        attributes.put("sigma", "" + welford.getStdDev());
        attributes.put("range", "" + welford.getRange());
        attributes.put("min", "" + welford.getMin());
        attributes.put("max", "" + welford.getMax());
        return attributes;
    }

    private static void iterateThruAndGrabPercentContact(MatrixZoomData zd, int maxBin, int minBin,
                                                         NormalizationType norm,
                                                         Map<SimpleLocation, Welford> results,
                                                         double[] vector1, double[] vector2,
                                                         Chromosome chrom, int resolution) {

        ExpectedModel poly = new LogExpectedSpline(zd, norm, chrom, resolution);

        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            if (cr.getCounts() > 0) {
                int dist = ExpectedUtils.getDist(cr);
                if (vector1[cr.getBinX()] > 1 && vector1[cr.getBinY()] > 1 && vector2[cr.getBinX()] > 1 && vector2[cr.getBinY()] > 1) {
                    if (dist > minBin && dist < maxBin) {
                        float percentContact = poly.getPercentContact(cr);
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
}
