package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.HiCUtils;
import cli.utils.RecapTools;
import cli.utils.flags.RegionConfiguration;
import cli.utils.flags.Utils;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.expected.QuickMedian;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Recap {

    private NormalizationType norm;

    public Recap(String[] args, CommandLineParser parser) {


        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5);
        }

        // recap <loops.bedpe> <outfolder> <file1.hic,file2.hic,...> <name1,name2,...>
        String loopListPath = args[1];
        File outFolder = UNIXTools.makeDir(new File(args[2]));

        String[] filepaths = args[3].split(",");
        String[] names = args[4].split(",");

        Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[0], false, false,
                true);

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                false, null, false);

        String possibleNorm = parser.getNormalizationStringOption();
        try {
            if (possibleNorm != null && possibleNorm.length() > 0) {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            } else {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR"});
            }
        } catch (Exception e) {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "VC_SQRT", "VC"});
        }
        System.out.println("Using normalization: " + norm.getLabel());

        ds.clearCache(false);
        ds = null;

        int numRows = loopList.getNumTotalFeatures();
        System.out.println("Number of loops: " + numRows);

        int resolution = parser.getResolutionOption(1000);
        if (resolution < 1) {
            resolution = 1000;
        }

        int window = parser.getWindowSizeOption(0);
        if (window < 1) {
            window = 5;
        }

        System.out.println("Using resolution: " + resolution);

        Feature2DList refinedLoops = recapStats(filepaths, names, loopList, handler, resolution, window, norm);
        refinedLoops.exportFeatureList(new File(outFolder, "recap.bedpe"), false, Feature2DList.ListFormat.NA);
        RecapTools.exportAllMatrices(refinedLoops, names, outFolder);
        System.out.println("pinpoint complete");
    }

    private static Feature2DList recapStats(String[] filepaths, String[] names, Feature2DList loopList,
                                            ChromosomeHandler handler, int resolution, int window,
                                            NormalizationType norm) {

        if (Main.printVerboseComments) {
            System.out.println("Start Recap/Compile process");
        }

        HiCZoom zoom = new HiCZoom(resolution);
        int matrixWidth = 1 + 2 * window;

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), false);

        int numTotalLoops = loopList.getNumTotalFeatures();
        final Object key = new Object();

        for (int di = 0; di < filepaths.length; di++) {

            final Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[di], false, false,
                    true);
            String prefix = names[di] + "_";

            final AtomicInteger currChromPair = new AtomicInteger(0);
            final AtomicInteger currNumLoops = new AtomicInteger(0);

            ParallelizationTools.launchParallelizedCode(() -> {

                int threadPair = currChromPair.getAndIncrement();
                while (threadPair < chromosomePairCounter) {
                    RegionConfiguration config = chromosomePairs.get(threadPair);
                    Chromosome chrom1 = config.getChr1();
                    Chromosome chrom2 = config.getChr2();

                    List<Feature2D> loops = loopList.get(chrom1.getIndex(), chrom2.getIndex());
                    if (loops != null && loops.size() > 0) {
                        MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom1, chrom2, zoom);

                        if (zd == null) {
                            System.err.println("ZD is null " + chrom1.getName() + "_" + chrom2.getName());
                            System.exit(9);
                        }

                        int maxBinDist = Math.max(getMaxDistance(loops, resolution, window), 9000000 / resolution);
                        double[] expected;
                        try {
                            expected = ds.getExpectedValues(
                                    zoom, norm, false).getExpectedValuesWithNormalization(
                                    chrom1.getIndex()).getValues().get(0);
                        } catch (Exception e) {
                            expected = null;
                        }
                        if (expected == null) {
                            System.out.println("Calculating log-normal expected vector");
                            expected = ExpectedUtils.calculateExpected(zd, norm, maxBinDist, true);
                        }

                        float pseudoCount = getMedianExpectedAt(maxBinDist - 2 * window, window, expected);
                        double superDiagonal = expected[1];

                        try {
                            for (Feature2D loop : loops) {
                                float[][] obsMatrix = new float[matrixWidth][matrixWidth];
                                Utils.addLocalizedData(obsMatrix, zd, loop, matrixWidth, resolution, window, norm, key);
                                float[][] eMatrix = new float[matrixWidth][matrixWidth];
                                Utils.fillInExpectedMatrix(eMatrix, loop, matrixWidth, expected, chrom1.getIndex(),
                                        resolution, window);
                                // MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_raw.npy")).getAbsolutePath(), output);

                                Map<String, String> attributes = RecapTools.getStats(obsMatrix, eMatrix,
                                        window, superDiagonal, pseudoCount);
                                for (String akey : attributes.keySet()) {
                                    loop.addStringAttribute(prefix + akey, attributes.get(akey));
                                }

                                if (currNumLoops.incrementAndGet() % 100 == 0) {
                                    System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                }
                            }
                            zd.clearCache();

                            System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                        } catch (Exception e) {
                            e.printStackTrace();
                            System.exit(76);
                        }
                    }
                    threadPair = currChromPair.getAndIncrement();
                }
            });
            ds.clearCache(false);
        }

        return loopList;
    }

    private static int getMaxDistance(List<Feature2D> loops, int resolution, int window) {
        long maxDist = 0;
        for (Feature2D loop : loops) {
            long dist = Math.abs((loop.getStart1() / resolution - window) - (loop.getEnd2() / resolution + window));
            if (dist > maxDist) {
                maxDist = dist;
            }
        }
        return (int) (maxDist + 4 * window);
    }

    private static float getMedianExpectedAt(int d0, int dx, double[] expectedVector) {
        List<Double> values = new ArrayList<>();
        for (int dist = d0 - dx; dist < d0 + dx; dist++) {
            double expected = expectedVector[dist];
            if (expected > 0) {
                values.add(expected);
            }
        }
        float median = QuickMedian.fastMedian(values);
        if (median > 0) return median;
        return 0f;
    }
}
