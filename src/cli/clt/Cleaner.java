package cli.clt;

import cli.Main;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import cli.utils.general.VectorCleaner;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Cleaner {

    private static final int MIN_LOOP_SIZE = 25000;

    public static void run(String[] args) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[1], false, true, true);
        String bedpeFile = args[2];
        String outFile = args[3];

        ChromosomeHandler handler = dataset.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(bedpeFile, handler,
                true, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        Feature2DList cleanList = cleanupLoops(dataset, loopList, handler);
        cleanList.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.FINAL);
    }

    private static Feature2DList cleanupLoops(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler) {

        int resolution = 1000;
        HiCZoom zoom = new HiCZoom(resolution);

        NormalizationType vcNorm = dataset.getNormalizationHandler().getNormTypeFromString("VC");

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), true);

        final AtomicInteger currChromPair = new AtomicInteger(0);
        final Object key = new Object();
        final Feature2DList goodLoopsList = new Feature2DList();

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                if (loops != null && loops.size() > 0) {
                    Set<Feature2D> goodLoops = new HashSet<>();

                    double[] vector1b = dataset.getNormalizationVector(chr1.getIndex(), zoom, vcNorm).getData().getValues().get(0);
                    VectorCleaner.inPlaceClean(vector1b);

                    double[] vector2b = vector1b;
                    if (chr1.getIndex() != chr2.getIndex()) {
                        vector2b = dataset.getNormalizationVector(chr2.getIndex(), zoom, vcNorm).getData().getValues().get(0);
                        VectorCleaner.inPlaceClean(vector2b);
                    }

                    try {
                        for (Feature2D loop : loops) {
                            if (passesMinLoopSize(loop)) {
                                if (normHasHighValPixel(loop.getStart1(), loop.getEnd1(), resolution, vector1b)
                                        && normHasHighValPixel(loop.getStart2(), loop.getEnd2(), resolution, vector2b)) {
                                    goodLoops.add(loop);
                                }
                            }
                        }
                    } catch (Exception e) {
                        System.err.println(e.getMessage());
                    }

                    synchronized (key) {
                        goodLoopsList.addByKey(Feature2DList.getKey(chr1, chr2), new ArrayList<>(goodLoops));
                    }
                }
                //System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus.get())) + "% ");
                threadPair = currChromPair.getAndIncrement();
            }
        });

        return goodLoopsList;
    }

    public static boolean passesMinLoopSize(Feature2D loop) {
        return dist(loop) > MIN_LOOP_SIZE;
    }

    public static long dist(Feature2D loop) {
        return (Math.min(Math.abs(loop.getEnd1() - loop.getStart2()),
                Math.abs(loop.getMidPt1() - loop.getMidPt2())) / MIN_LOOP_SIZE) * MIN_LOOP_SIZE;
    }

    private static boolean normHasHighValPixel(long start, long end, int resolution, double[] vector) {

        int x0 = (int) (start / resolution);
        int xF = (int) (end / resolution) + 1;

        for (int k = x0; k < xF + 1; k++) {
            if (vector[k] > 1 && vector[k] >= vector[k - 1] && vector[k] >= vector[k + 1]) {
                return true;
            }
        }

        return false;
    }
}
