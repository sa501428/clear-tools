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
import javastraw.reader.norm.NormalizationPicker;
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
        Dataset ds = HiCFileTools.extractDatasetForCLT(files[0], false, false, true);
        int resolution = parser.getResolutionOption(5000);
        String norm = parser.getNormalizationStringOption();

        final Feature2DList result = new Feature2DList();

        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        AtomicInteger cIndex = new AtomicInteger(0);

        // this code can be commented out when running small-scale tests on local machine

        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chrom = chromosomes[currIndex];
                List<Feature2D> hotspots = findTheHotspots(chrom, files, resolution, norm);
                synchronized (result) {
                    result.addByKey(Feature2DList.getKey(chrom, chrom), hotspots);
                }
                currIndex = cIndex.getAndIncrement();
            }
        });

        result.exportFeatureList(new File(outputFileName + ".hotspot.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("hotspot complete");
    }

    private static List<Feature2D> findTheHotspots(Chromosome chrom, String[] files, int resolution,
                                                   String normStringOption) {
        // 2 ints (positions) row and column
        Map<SimpleLocation, Welford> results = new HashMap<>();
        List<SimpleLocation> removeList = new ArrayList<>();
        List<Feature2D> hotspots = new ArrayList<>();
        Map<String, String> attributes = new HashMap<>();

        //int NUM_NONZERO_VALUES_THRESHOLD = Math.max(files.length / 2, 2);
        //System.out.println("validCountThreshold: " + NUM_NONZERO_VALUES_THRESHOLD);

        int minBin = MIN_DIST / resolution;
        int maxBin = MAX_DIST / resolution;

        for (String file : files) {
            Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, false, false);

            NormalizationType norm = normSelector(ds, normStringOption);

            Matrix matrix = ds.getMatrix(chrom, chrom);
            if (matrix == null) {
                System.out.println("matrix for " + chrom.getName() + " of file " + file + " == null -> continuing");
                continue;
            }
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) {
                System.out.println("zd for " + chrom.getName() + " of file " + file + " == null -> continuing");
                continue;
            }

            // iterating through chrom using type 1 iteration
            iterateThruAllTheValues(zd, maxBin, minBin, norm, results);
        }

        System.out.println("Total recorded locations before removal: " + results.size());

        for (SimpleLocation key : results.keySet()) {
            Welford value = results.get(key);
            //value.addZeroIfBelow(files.length);
            if (value.getRange() < .005) {
                removeList.add(key);
            }
            //if (entry.getValue().getCounts() < NUM_NONZERO_VALUES_THRESHOLD)
            //    removeList.add(entry.getKey());
        }
        // test print
        System.out.println("removeList size: " + removeList.size());

        for (SimpleLocation key : removeList) {
            results.remove(key);
        }
        // test print
        System.out.println("num. remaining entries after removal: " + results.size());
        removeList.clear();

        Zscore zscore = getOverallZscoreMetric(results.values());

        for (Map.Entry<SimpleLocation, Welford> entry : results.entrySet()) {
            Welford welford = entry.getValue();
            if (zscore.getZscore(welford.getStdDev()) >= ZSCORE_CUTOFF) {
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

        System.out.println("final number of hotspots: " + hotspots.size());
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


    public static NormalizationType normSelector(Dataset ds, String normStringOption) {
        NormalizationType norm;
        if (normStringOption != null && normStringOption.length() > 0) {
            try {
                norm = ds.getNormalizationHandler().getNormTypeFromString(normStringOption);
            } catch (Exception e) {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normStringOption, "SCALE", "KR", "NONE"});
            }
        } else {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "NONE"});
        }
        System.out.println("Norm being used: " + norm.getLabel());
        return norm;
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
        System.out.println("finished recording locations for this file");
    }
}
