package cli.clt.loops;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.FeatureStats;
import cli.utils.data.SparseContactMatrixWithMasking;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import cli.utils.general.QuickGrouping;
import cli.utils.general.ZscoreTools;
import javastraw.expected.LogExpectedZscoreSpline;
import javastraw.expected.Welford;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
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

    private static final int zHighCutoff = 2;
    private static final float oeHighCutoff = (float) Math.log(2);
    // [-strict][-peek]
    // ; peek just saves values\n\t\tstrict requires each resolution to meet the criteria
    public static String usage = "sieve[-easy][-skip-global] [--threads num_threads] [-k NORM] <loops.bedpe> <out.stem> <file.hic> [res1,...]\n" +
            "\t\tretain loop if at a loop-y location\n" +
            "\t\tsieve-post-filter <loops.bedpe> <out.stem> <genomeID>";
    private static int zLowCutoff = 1;
    public static String GLOBAL_OE = "_sieve_obs_over_global_expected";
    public static String GLOBAL_PERCENT = "_sieve_global_percent_contact";
    public static String LOCAL_OE = "_sieve_obs_over_local_expected";
    public static String GLOBAL_Z = "_sieve_global_zscore";
    public static String LOCAL_Z = "_sieve_local_zscore";
    private static float oeLowCutoff = (float) Math.log(1.5);
    private static int numThreads = 8;

    public int[] resolutions = new int[]{1000, 2000, 5000, 10000};
    private static boolean skipGlobal = false;

    public Sieve(String[] args, CommandLineParser parser, String command) {
        // sieve <loops.bedpe> <output.bedpe> <file1.hic> <res1,res2,...>
        if (args.length != 5 && args.length != 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        String loopListPath = args[1];
        String outStem = args[2];
        String hicPath = args[3];

        if (command.contains("easy")) {
            zLowCutoff = 0;
            oeLowCutoff = 0;
        }
        if (command.contains("skip") && command.contains("global")) {
            skipGlobal = true;
        }

        numThreads = parser.getNumThreads(numThreads);

        if (command.contains("post")) {
            // just filter using values in list
            ChromosomeHandler handler = ChromosomeTools.loadChromosomes(hicPath);
            Feature2DList result = Feature2DParser.loadFeatures(loopListPath, handler, true, null, false);
            Feature2DList[] goodAndBad = filterByScore(result);
            goodAndBad[0].exportFeatureList(new File(outStem + ".good.loops.bedpe"), false, Feature2DList.ListFormat.NA);
            goodAndBad[1].exportFeatureList(new File(outStem + ".weak.loops.bedpe"), false, Feature2DList.ListFormat.NA);
            goodAndBad[2].exportFeatureList(new File(outStem + ".not.loops.bedpe"), false, Feature2DList.ListFormat.NA);
        } else {

            if (args.length > 4) {
                resolutions = parseInts(args[4]);
            }

            Dataset ds = HiCFileTools.extractDatasetForCLT(hicPath, false, false, true);
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

            int window = parser.getWindowSizeOption(10);
            if (window < 5) {
                window = 10;
            }

            Feature2DList result = sieveFilter(ds, loopList, handler, resolutions, window, norm);
            result.exportFeatureList(new File(outStem + ".attributes.bedpe"), false, Feature2DList.ListFormat.NA);

            if (command.contains("easy")) {
                Feature2DList goodCalls = filterBySimpleScore(result);
                goodCalls.exportFeatureList(new File(outStem + ".minimal.loops.bedpe"), false, Feature2DList.ListFormat.NA);
            }
        }
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

        ParallelizationTools.launchParallelizedCode(numThreads, () -> {
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

                            NormalizationVector nv1 = ds.getNormalizationVector(chrom1.getIndex(), zoom, norm);
                            NormalizationVector nv2 = ds.getNormalizationVector(chrom2.getIndex(), zoom, norm);

                            if (nv1 == null) {
                                System.err.println("Error getting normalization " + norm.getLabel() + " for " + chrom1.getName());
                            } else if (nv2 == null) {
                                System.err.println("Error getting normalization " + norm.getLabel() + " for " + chrom2.getName());
                            } else if (loopsToAssessGlobal.size() > 0) {
                                setDefaultAttributes(loopsToAssessGlobal, resolution);
                                Collection<List<Feature2D>> loopGroups = QuickGrouping.groupNearbyRecords(
                                        loopsToAssessGlobal, 500 * resolution).values();

                                LogExpectedZscoreSpline poly = null;
                                if (!skipGlobal) {
                                    poly = new LogExpectedZscoreSpline(zd, norm, chrom1, resolution);
                                }
                                SparseContactMatrixWithMasking sparseMatrix = new SparseContactMatrixWithMasking(zd,
                                        loopsToAssessGlobal, resolution, buffer, 2 * buffer + 1, norm);

                                for (List<Feature2D> group : loopGroups) {
                                    int minR = (int) ((FeatureStats.minStart1(group) / resolution) - buffer);
                                    int minC = (int) ((FeatureStats.minStart2(group) / resolution) - buffer);
                                    int maxR = (int) ((FeatureStats.maxEnd1(group) / resolution) + buffer);
                                    int maxC = (int) ((FeatureStats.maxEnd2(group) / resolution) + buffer);
                                    float[][] regionMatrix = sparseMatrix.getRegion(minR, minC, maxR, maxC);
                                    for (Feature2D loop : group) {
                                        int absCoordBinX = (int) (loop.getMidPt1() / resolution);
                                        int absCoordBinY = (int) (loop.getMidPt2() / resolution);
                                        int dist = Math.abs(absCoordBinX - absCoordBinY);
                                        int midX = absCoordBinX - minR;
                                        int midY = absCoordBinY - minC;

                                        float observed = regionMatrix[midX][midY];

                                        Welford localWelford = ZscoreTools.getLocalWelford(regionMatrix, midX, midY, window);
                                        float localOE = (float) (observed / localWelford.getMean());
                                        float localZScore = (float) localWelford.getZscore().getZscore(observed);
                                        loop.addStringAttribute(resolution + LOCAL_Z, "" + localZScore);
                                        loop.addStringAttribute(resolution + LOCAL_OE, "" + localOE);

                                        if (poly != null) {
                                            float globalOE = (float) (observed / poly.getExpectedFromUncompressedBin(dist));
                                            float globalZScore = (float) poly.getZscoreForObservedUncompressedBin(dist, observed);
                                            loop.addStringAttribute(resolution + GLOBAL_Z, "" + globalZScore);
                                            loop.addStringAttribute(resolution + GLOBAL_OE, "" + globalOE);

                                            loop.addStringAttribute(resolution + GLOBAL_PERCENT, "" + poly.getPercentContact(dist, observed));
                                        }
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

    private static void setDefaultAttributes(Set<Feature2D> loops, int resolution) {
        for (Feature2D loop : loops) {
            loop.addStringAttribute(resolution + LOCAL_OE, "NaN");
            loop.addStringAttribute(resolution + LOCAL_Z, "NaN");
            if (!skipGlobal) {
                loop.addStringAttribute(resolution + GLOBAL_OE, "NaN");
                loop.addStringAttribute(resolution + GLOBAL_Z, "NaN");
            }
        }
    }

    private Feature2DList filterBySimpleScore(Feature2DList result) {
        Feature2DList good = new Feature2DList();
        result.processLists((s, list) -> {
            List<Feature2D> goodLoops = new ArrayList<>();
            for (Feature2D feature : list) {
                if (isSimpleEnrichedLikelyLoop(feature, resolutions)) {
                    goodLoops.add(feature);
                }
            }
            good.addByKey(s, goodLoops);
        });
        return good;
    }

    private Feature2DList[] filterByScore(Feature2DList result) {
        Feature2DList good = new Feature2DList();
        Feature2DList weak = new Feature2DList();
        Feature2DList bad = new Feature2DList();
        result.processLists((s, list) -> {
            List<Feature2D> goodLoops = new ArrayList<>();
            List<Feature2D> weakLoops = new ArrayList<>();
            List<Feature2D> badLoops = new ArrayList<>();
            for (Feature2D feature : list) {
                if (isEnrichedLikelyLoop(feature)) {
                    goodLoops.add(feature);
                } else if (isWeaklyEnrichedMaybeLoopish(feature)) {
                    weakLoops.add(feature);
                } else {
                    badLoops.add(feature);
                }
            }
            good.addByKey(s, goodLoops);
            weak.addByKey(s, weakLoops);
            bad.addByKey(s, badLoops);
        });
        return new Feature2DList[]{good, weak, bad};
    }

    private boolean isEnrichedLikelyLoop(Feature2D feature) {

        double localZ1K = getAttribute(feature, 1000 + LOCAL_Z, -10);
        double localZ2K = getAttribute(feature, 2000 + LOCAL_Z, -10);
        double localZ5K = getAttribute(feature, 5000 + LOCAL_Z, -10);

        double localOE1K = Math.log(getAttribute(feature, 1000 + LOCAL_OE, .1f));
        double localOE2K = Math.log(getAttribute(feature, 2000 + LOCAL_OE, .1f));
        double localOE5K = Math.log(getAttribute(feature, 5000 + LOCAL_OE, .1f));

        double globalOE1K = Math.log(getAttribute(feature, 1000 + GLOBAL_OE, .1f));
        double globalOE2K = Math.log(getAttribute(feature, 2000 + GLOBAL_OE, .1f));
        double globalOE5K = Math.log(getAttribute(feature, 5000 + GLOBAL_OE, .1f));

        boolean weak1 = (localZ1K > zLowCutoff) && (localOE1K > oeLowCutoff);
        boolean weak2 = (localZ2K > zLowCutoff) && (localOE2K > oeLowCutoff);
        boolean weak5 = (localZ5K > zLowCutoff) && (localOE5K > oeLowCutoff);

        boolean medium1 = weak1 && (globalOE1K > oeLowCutoff);
        boolean medium2 = weak2 && (globalOE2K > oeLowCutoff);
        boolean medium5 = weak5 && (globalOE5K > oeLowCutoff);

        boolean strong1 = (localZ1K > zHighCutoff) && (localOE1K > oeHighCutoff) && (globalOE1K > oeHighCutoff);
        boolean strong2 = (localZ2K > zHighCutoff) && (localOE2K > oeHighCutoff) && (globalOE2K > oeHighCutoff);
        boolean strong5 = (localZ5K > zHighCutoff) && (localOE5K > oeHighCutoff) && (globalOE5K > oeHighCutoff);

        boolean localZAvg = (localZ1K + localZ2K + localZ5K) / 3 > 1;
        boolean localOEAvg = (localOE1K + localOE2K + localOE5K) / 3 > oeHighCutoff;
        boolean hasAverageEnrichment = localZAvg && localOEAvg;

        return (strong1 && (medium2 || medium5 || hasAverageEnrichment))
                || (strong2 && (medium1 || medium5 || hasAverageEnrichment))
                || (strong5 && (medium1 || medium2 || hasAverageEnrichment))
                || (hasAverageEnrichment && ((medium1 && medium2) || (weak2 && weak5) || (weak1 && weak5)));
    }

    private boolean isSimpleEnrichedLikelyLoop(Feature2D feature, int[] resolutions) {
        boolean isEnriched = true;
        for (int res : resolutions) {
            double localZ = getAttribute(feature, res + LOCAL_Z, -10);
            double localOE = Math.log(getAttribute(feature, res + LOCAL_OE, .1f));
            boolean weak = (localZ > zLowCutoff) && (localOE > oeLowCutoff);
            isEnriched = isEnriched && weak;
        }
        return isEnriched;
    }

    private boolean isWeaklyEnrichedMaybeLoopish(Feature2D feature) {
        for (int res : resolutions) {
            //float zScore = getAttribute(feature, res + LOCAL_Z, 0);
            double oe1 = Math.log(getAttribute(feature, res + GLOBAL_OE, 0.1f));
            double oe2 = Math.log(getAttribute(feature, res + LOCAL_OE, 0.1f));
            double zVal = getAttribute(feature, res + LOCAL_Z, -10);
            if (oe1 > 1.25 || oe2 > 1.1 || zVal > 0.5) {
                return true;
            }
        }
        return false;
    }

    private float getAttribute(Feature2D feature, String key, float defaultValue) {
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
