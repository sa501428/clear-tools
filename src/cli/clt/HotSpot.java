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

//    private static void getMatrices(int resolution, int window, String strNorm, String file,
//                                    Map<Integer, Map<String, float[][]>> results) {
//        /*
//        Accepts parameters of the HOTSPOT command line tool and returns a list of 2D float arrays that represent the Hi-C maps of a (temporarily) pre-defined region
//        corresponding to the files. The 2D float arrays will have pixel sizes equal window argument
//         */
//        boolean useCache = false;
//        // create a hic dataset object
//        Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, useCache, true);
//        // choose norm: we know the datasets we're using will have SCALE available
//        // NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{strNorm, "NONE"});
//        NormalizationType norm = NormalizationHandler.SCALE;
//
//        // Instantiates the 2D float array that will be used to represent the Hi-C map of individual tissue genomes
//
//        // todo remove this hardcode. This is just used for displaying a singular region of interest as a matrix.
//        //  Proof of concept that stdDeviationFinder works and matrix prints.
//        //  Eventually want to slide and record non-zero (or higher std dev?) values to bedpe file
//
//        Chromosome c5 = ds.getChromosomeHandler().getChromosomeFromName("chr5");
//        Matrix matrix = ds.getMatrix(c5, c5);
//        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
//        int binXStart = 119000000 / resolution;
//        int binXEnd = binXStart + window;
//        int binYStart = 119000000 / resolution;
//        int binYEnd = binYStart + window;
//        try {
//            ExpectedValueFunction df = ds.getExpectedValuesOrExit(zd.getZoom(), norm, c5, true, false);
//            float[][] currentWindowMatrix = ExtractingOEDataUtils.extractObsOverExpBoundedRegion(zd,
//                    binXStart, binXEnd, binYStart, binYEnd, window, window, norm, df, c5.getIndex(), 50, true,
//                    true, ExtractingOEDataUtils.ThresholdType.TRUE_OE, 1, 0);
//
//            results.get(5).put(file, currentWindowMatrix);
//        } catch (IOException e) {
//            System.err.println("error extracting local bounded region float matrix");
//        }

    //        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
//
//        for (int i = 0; i < chromosomes.length; i++) {
//            Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[i]);
//            if (matrix == null) continue;
//            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
//            if (zd == null) continue;
//
//            long binXStart = 0; // todo make this slide across the map
//            long binXEnd = binXStart + window;
//            long binYStart = 0;
//            long binYEnd = binXStart + window;
//
//            float[][] currentMatrix = new float[window][window];
//            results.get(chromosomes[i].getIndex()).put(file, currentMatrix);
//
//            // Iterates through the blocks of the chromosome pair, eventually grabbing the bin coordinates and count to record them in currentMatrix
//            Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
//            while (iterator.hasNext()) {
//                ContactRecord record = iterator.next();
//                int binX = record.getBinX();
//                int binY = record.getBinY();
//                currentMatrix[binX][binY] = record.getCounts();
//            }
//        }
//    }
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

        final int DEFAULT_RES = 2000;
        String[] files = args[1].split(",");
        String outputFileName = args[2];
        Dataset ds = HiCFileTools.extractDatasetForCLT(files[0], false, false, true);
        int resolution = parser.getResolutionOption(DEFAULT_RES);
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

        // this code is used for running small-scale tests on local machine
//        Chromosome c21 = ds.getChromosomeHandler().getChromosomeFromName("chr21");
//        List<Feature2D> hotspots = findTheHotspots(c21, files, parser.getResolutionOption(DEFAULT_RES),
//                parser.getNormalizationStringOption());
//        result.addByKey(Feature2DList.getKey(c21, c21), hotspots);

        result.exportFeatureList(new File(outputFileName + ".hotspot.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("hotspot complete");
    }

    private static List<Feature2D> findTheHotspots(Chromosome chrom, String[] files, int resolutionOption,
                                                   String normStringOption) {
        // 2 ints (positions) row and column
        Map<SimpleLocation, Welford> results = new HashMap<>();
        List<SimpleLocation> removeList = new ArrayList();
        List<Feature2D> hotspots = new ArrayList();
        Map<String, String> attributes = new HashMap();

        int NUM_NONZERO_VALUES_THRESHOLD = Math.max(files.length / 2, 2);
        System.out.println("validCountThreshold: " + NUM_NONZERO_VALUES_THRESHOLD);

        int minBin = MIN_DIST / resolutionOption;
        int maxBin = MAX_DIST / resolutionOption;

        ////// for every dataset, extract the chromosome
        // iterate on it type 1

        // load the dataset // iterating on them 1 at a time
        // iteration type Option #1

        // any contact more than 10MB from the diagonal can be skipped / ignored
        // also skip any contact < 25kb
        // only do this for contacts where the counts > 0

        // inside where the contacts are,


        ////// END for every dataset, extract the chromosome


        // iterate and delete any Welford less than (3 for now, but let this be a variable that we can test) values

        for (String file : files) {
            Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, false, false);

            NormalizationType norm = normSelector(ds, normStringOption);

            Matrix matrix = ds.getMatrix(chrom, chrom);
            if (matrix == null) {
                System.out.println("matrix for " + chrom.getName() + " of file " + file + " == null -> continuing");
                continue;
            }
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolutionOption));
            if (zd == null) {
                System.out.println("zd for " + chrom.getName() + " of file " + file + " == null -> continuing");
                continue;
            }

            // iterating through chrom using type 1 iteration
            iterateThruAllTheValues(zd, maxBin, minBin, norm, results);
            // test print
            System.out.println("Total recorded locations before removal: " + results.size());

            for (Map.Entry<SimpleLocation, Welford> entry : results.entrySet()) {
                if (entry.getValue().getCounts() < NUM_NONZERO_VALUES_THRESHOLD)
                    removeList.add(entry.getKey());
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
                    attributes.put("std", "" + welford.getStdDev());
                    long startX = (long) entry.getKey().getBinX() * resolutionOption;
                    long endX = startX + resolutionOption;
                    long startY = (long) entry.getKey().getBinY() * resolutionOption;
                    long endY = startY + resolutionOption;
                    Feature2D feature = new Feature2D(Feature2D.FeatureType.PEAK, chrom.getName(), startX, endX, chrom.getName(), startY, endY, Color.BLACK, attributes);
                    hotspots.add(feature);
                }
            }
        }

        System.out.println("final number of hotspots: " + hotspots.size());
        return hotspots;
    }

    private static Zscore getOverallZscoreMetric(Collection<Welford> values) {
        Welford overallWelford = new Welford();
        for (Welford welford : values) {
            overallWelford.addValue(welford.getStdDev());
        }
        System.out.println("overall welford standard deviation: " + overallWelford.getStdDev());
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

                    float percentContact0 = expected.getPercentContact(dist, cr.getCounts());

                    float percentContact1 = Math.min(1, Math.max(0, percentContact0));
                    //percentContact2 = Math.exp(percentContact0 - 1);

                    SimpleLocation location = new SimpleLocation(cr);
                    results.putIfAbsent(location, new Welford());
                    results.get(location).addValue(percentContact1);
                }
            }
        }
        System.out.println("finished recording locations for this file");
    }
}
