package cli.clt;

import cli.Main;
import cli.utils.FeatureStats;
import cli.utils.clean.BinCollisionChecker;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import cli.utils.general.QuickGrouping;
import cli.utils.general.Utils;
import cli.utils.general.ZscoreTools;
import javastraw.expected.ExpectedModel;
import javastraw.expected.LogExpectedSpline;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Sieve {

    // [-strict][-peek]
    // ; peek just saves values\n\t\tstrict requires each resolution to meet the criteria
    public static String usage = "sieve [-k NORM] <loops.bedpe> <out.stem> <file.hic> [res1,...]\n" +
            "\t\tretain loop if at a loop-y location";

    public Sieve(String[] args, CommandLineParser parser, String command) {
        // sieve <loops.bedpe> <output.bedpe> <file1.hic> <res1,res2,...>
        if (args.length != 5 && args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        String loopListPath = args[1];
        String outStem = args[2];
        String filepath = args[3];
        int[] resolutions = new int[]{1000, 2000, 5000};
        if (args.length > 4) {
            resolutions = parseInts(args[4]);
        }

        Dataset ds = HiCFileTools.extractDatasetForCLT(filepath, false, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler, true, null, false);

        String possibleNorm = parser.getNormalizationStringOption();
        NormalizationType norm = NormalizationHandler.VC;
        if (possibleNorm != null && possibleNorm.length() > 0) {
            if (possibleNorm.equalsIgnoreCase("none")) {
                norm = NormalizationHandler.NONE;
            } else {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            }
        }
        System.out.println("Using normalization: " + norm.getLabel());

        int window = parser.getWindowSizeOption(0);
        if (window < 2) {
            window = 5;
        }

        Feature2DList result = sieveFilter(ds, loopList, handler, resolutions, window, norm);
        result.exportFeatureList(new File(outStem + ".attributes.bedpe"), false, Feature2DList.ListFormat.NA);

        Feature2DList[] goodAndBad = filterByScore(result);
        goodAndBad[0].exportFeatureList(new File(outStem + ".good.loops.bedpe"), false, Feature2DList.ListFormat.NA);
        goodAndBad[1].exportFeatureList(new File(outStem + ".not.loops.bedpe"), false, Feature2DList.ListFormat.NA);

        System.out.println("sieve complete");
    }

    private static Feature2DList sieveFilter(Dataset ds, Feature2DList loopList,
                                             ChromosomeHandler handler, int[] resolutions, int window,
                                             NormalizationType norm) {

        if (Main.printVerboseComments) {
            System.out.println("Start Sieve process");
        }

        final Feature2DList newLoopList = new Feature2DList();
        int buffer = 2 * window;

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), false);

        final AtomicInteger currChromPair = new AtomicInteger(0);
        final AtomicInteger numLoopsDone = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(() -> {
            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chrom1 = config.getChr1();
                Chromosome chrom2 = config.getChr2();

                Set<Feature2D> loopsToAssessGlobal = new HashSet<>(loopList.get(chrom1.getIndex(), chrom2.getIndex()));
                Matrix matrix = ds.getMatrix(chrom1, chrom2);
                int numLoopsForThisChromosome = loopsToAssessGlobal.size();

                if (matrix != null) {
                    for (int resolution : resolutions) {
                        HiCZoom zoom = new HiCZoom(resolution);
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            Set<Feature2D> loopsToAssessThisRound = filterForAppropriateResolution(loopsToAssessGlobal, resolution);

                            if (loopsToAssessThisRound.size() > 0) {
                                Collection<List<Feature2D>> loopGroups = QuickGrouping.groupNearbyRecords(
                                        loopsToAssessThisRound, 500 * resolution).values();

                                ExpectedModel poly = new LogExpectedSpline(zd, norm, chrom1, resolution);

                                for (List<Feature2D> group : loopGroups) {
                                    int minR = (int) ((FeatureStats.minStart1(group) / resolution) - buffer);
                                    int minC = (int) ((FeatureStats.minStart2(group) / resolution) - buffer);
                                    int maxR = (int) ((FeatureStats.maxEnd1(group) / resolution) + buffer);
                                    int maxC = (int) ((FeatureStats.maxEnd2(group) / resolution) + buffer);
                                    float[][] regionMatrix = Utils.getRegion(zd, minR, minC, maxR, maxC, norm);
                                    for (Feature2D loop : group) {
                                        int absCoordBinX = (int) (loop.getMidPt1() / resolution);
                                        int absCoordBinY = (int) (loop.getMidPt2() / resolution);
                                        int dist = Math.abs(absCoordBinX - absCoordBinY);
                                        int midX = absCoordBinX - minR;
                                        int midY = absCoordBinY - minC;

                                        float zScore = (float) ZscoreTools.getLocalZscore(regionMatrix, midX, midY, window);
                                        float observed = regionMatrix[midX][midY];
                                        float oe = (float) (observed / poly.getExpectedFromUncompressedBin(dist));

                                        loop.addStringAttribute(resolution + "_sieve_obs_over_expected", "" + oe);
                                        loop.addStringAttribute(resolution + "_sieve_local_zscore", "" + zScore);
                                    }
                                    regionMatrix = null;
                                }
                                System.out.print(".");
                            }
                        }
                        matrix.clearCacheForZoom(zoom);
                    }
                    matrix.clearCache();
                }

                synchronized (newLoopList) {
                    newLoopList.addByKey(Feature2DList.getKey(chrom1, chrom2),
                            new ArrayList<>(loopsToAssessGlobal));
                }

                if (numLoopsForThisChromosome > 0) {
                    int num = numLoopsDone.addAndGet(numLoopsForThisChromosome);
                    System.out.println("\n" + chrom1.getName() + " done\nNumber of loops processed overall: " + num);
                }

                threadPair = currChromPair.getAndIncrement();
            }
        });

        return newLoopList;
    }

    private static Set<Feature2D> filterForAppropriateResolution(Set<Feature2D> loops, int resolution) {
        Set<Feature2D> goodLoops = new HashSet<>();
        for (Feature2D loop : loops) {
            if (Math.max(loop.getWidth1(), loop.getWidth2()) <= resolution) {
                goodLoops.add(loop);
            }
        }
        return BinCollisionChecker.ensureOnlyOneLoopPerBin(goodLoops, resolution);
    }

    private Feature2DList[] filterByScore(Feature2DList result) {
        Feature2DList good = new Feature2DList();
        Feature2DList bad = new Feature2DList();
        result.processLists((s, list) -> {
            List<Feature2D> goodLoops = new ArrayList<>();
            List<Feature2D> badLoops = new ArrayList<>();
            for (Feature2D feature : list) {
                if (isEnrichedLocalAndGlobalMaxima(feature)) {
                    goodLoops.add(feature);
                } else {
                    badLoops.add(feature);
                }
            }
            good.addByKey(s, goodLoops);
            bad.addByKey(s, badLoops);
        });
        return new Feature2DList[]{good, bad};
    }

    private boolean isEnrichedLocalAndGlobalMaxima(Feature2D feature) {

        float zScore1kb = getAttribute(feature, "1000_sieve_local_zscore", 0);
        float zScore2kb = getAttribute(feature, "2000_sieve_local_zscore", 0);
        float zScore5kb = getAttribute(feature, "5000_sieve_local_zscore", 0);
        float oe1kb = getAttribute(feature, "1000_sieve_obs_over_expected", 0);
        float oe2kb = getAttribute(feature, "2000_sieve_obs_over_expected", 0);
        float oe5kb = getAttribute(feature, "5000_sieve_obs_over_expected", 0);

        boolean singleResStrongSignalPrediction = (zScore1kb > 2 && oe1kb > 2)
                || (zScore2kb > 2 && oe2kb > 2)
                || (zScore5kb > 2 && oe5kb > 2);

        boolean pz1 = zScore1kb > 1.5 && oe1kb > 1.5;
        boolean pz2 = zScore2kb > 1.5 && oe2kb > 1.5;
        boolean pz5 = zScore5kb > 1.5 && oe5kb > 1.5;
        boolean dualResWeakSignalPrediction = (pz1 && pz2) || (pz2 && pz5) || (pz1 && pz5);

        return singleResStrongSignalPrediction || dualResWeakSignalPrediction;
    }

    private float getAttribute(Feature2D feature, String key, int defaultValue) {
        if (feature.hasAttributeKey(key)) {
            try {
                return Float.parseFloat(feature.getAttribute(key));
            } catch (Exception ignored) {
            }
        }
        return defaultValue;
    }

    private int[] parseInts(String input) {
        String[] inputs = input.split(",");
        int[] values = new int[inputs.length];
        for (int i = 0; i < values.length; i++) {
            values[i] = Integer.parseInt(inputs[i]);
        }
        Arrays.sort(values);
        return values;
    }
}
