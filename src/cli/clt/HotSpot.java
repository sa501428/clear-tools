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
    private static final int MIN_DIST = 25000;

    // The Z-Score Cutoff is currently hardcoded at 1.65, which has a confidence interval of __%
    private static final float ZSCORE_CUTOFF = 2f;

    private static void printUsageAndExit() {
        /* example print: ("apa [--min-dist minval] [--max-dist max_val] [--window window] [-r resolution]" +
                " [-k NONE/VC/VC_SQRT/KR] [--corner-width corner_width] [--include-inter include_inter_chr] [--ag-norm]" +
                " <input.hic> <loops.bedpe> <outfile>"); */
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

        Dataset[] datasets = new Dataset[files.length];
        for (int k = 0; k < files.length; k++) {
            datasets[k] = HiCFileTools.extractDatasetForCLT(files[k],
                    false, false, false);
        }

        NormalizationType norm = datasets[0].getNormalizationHandler().getNormTypeFromString(normString);
        Chromosome[] chromosomes = datasets[0].getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        AtomicInteger cIndex = new AtomicInteger(0);
        // this code can be commented out when running small-scale tests on local machine

        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chrom = chromosomes[currIndex];
                List<Feature2D> hotspots = findTheHotspots(chrom, datasets, resolution, norm);
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
                                                   NormalizationType norm) {
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
            iterateThruAllTheValues(zd, maxBin, minBin, norm, results);
            matrix.clearCache();
        }

        List<SimpleLocation> removeList = new ArrayList<>();
        for (SimpleLocation key : results.keySet()) {
            double range = results.get(key).getRange();
            int counts = (int) results.get(key).getCounts();
            //value.addZeroIfBelow(files.length);
            if (range < .005 || range > 0.5 || counts < 2) {
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
            Zscore zscore = getOverallZscoreMetric(results.values());
            for (Map.Entry<SimpleLocation, Welford> entry : results.entrySet()) {
                Welford welford = entry.getValue();
                if (zscore.getZscore(welford.getStdDev()) >= ZSCORE_CUTOFF) {
                    Map<String, String> attributes = new HashMap<>();
                    attributes.put("sigma", "" + welford.getStdDev());
                    attributes.put("range", "" + welford.getRange());
                    long startX = (long) entry.getKey().getBinX() * resolution;
                    long endX = startX + resolution;
                    long startY = (long) entry.getKey().getBinY() * resolution;
                    long endY = startY + resolution;
                    Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, chrom.getName(), startX, endX, chrom.getName(), startY, endY, Color.BLACK, attributes);
                    hotspots.add(feature);
                }
            }
        }

        int n4 = hotspots.size();
        System.out.println("Hotspots " + chrom.getName() + " pre-filter: " + n1 +
                " removed: " + n2 + " post-filter: " + n3 + " final: " + n4);
        return hotspots;
    }

    private static Zscore getOverallZscoreMetric(Collection<Welford> values) {
        Welford overallWelford = new Welford();
        for (Welford welford : values) {
            //overallWelford.addValue(welford.getStdDev());
            overallWelford.addValue(welford.getRange());
        }
        System.out.println("Overall welford: " + overallWelford.getSummary());
        return overallWelford.getZscore();
    }

    private static void iterateThruAllTheValues(MatrixZoomData zd, int maxBin, int minBin,
                                                NormalizationType norm,
                                                Map<SimpleLocation, Welford> results) {

        LogExpectedModel expected = new LogExpectedModel(zd, norm, maxBin);

        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            if (cr.getCounts() > 0) {
                int dist = ExpectedUtils.getDist(cr);
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
        System.out.print(".");
        //System.out.println("Finished recording locations for this file");
    }
}
