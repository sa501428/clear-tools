package cli.clt;

import cli.Main;
import cli.utils.HiCUtils;
import cli.utils.VectorCleanerUtils;
import cli.utils.flags.RegionConfiguration;
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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Cleaner {

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

        NormalizationType scaleNorm = dataset.getNormalizationHandler().getNormTypeFromString("SCALE");
        NormalizationType vcNorm = dataset.getNormalizationHandler().getNormTypeFromString("VC");


        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getAutosomalChromosomesArray());

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
                    List<Feature2D> goodLoops = new ArrayList<>();

                    double[] vector1 = dataset.getNormalizationVector(chr1.getIndex(), zoom, scaleNorm).getData().getValues().get(0);
                    double[] vector1b = dataset.getNormalizationVector(chr1.getIndex(), zoom, vcNorm).getData().getValues().get(0);
                    VectorCleanerUtils.inPlaceClean(vector1);
                    VectorCleanerUtils.inPlaceClean(vector1b);

                    double[] vector2 = vector1;
                    double[] vector2b = vector1b;
                    if (chr1.getIndex() != chr2.getIndex()) {
                        vector2 = dataset.getNormalizationVector(chr2.getIndex(), zoom, scaleNorm).getData().getValues().get(0);
                        vector2b = dataset.getNormalizationVector(chr2.getIndex(), zoom, vcNorm).getData().getValues().get(0);
                        VectorCleanerUtils.inPlaceClean(vector2);
                        VectorCleanerUtils.inPlaceClean(vector2b);
                    }

                    try {
                        for (Feature2D loop : loops) {
                            if (normsAreGood(loop.getStart1(), loop.getWidth1(), resolution, vector1, vector1b)
                                    && normsAreGood(loop.getStart2(), loop.getWidth2(), resolution, vector2, vector2b)) {
                                goodLoops.add(loop);
                            }
                        }
                    } catch (Exception e) {
                        System.err.println(e.getMessage());
                    }

                    synchronized (key) {
                        goodLoopsList.addByKey(Feature2DList.getKey(chr1, chr2), goodLoops);
                    }
                }
                //System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus.get())) + "% ");
                threadPair = currChromPair.getAndIncrement();
            }
        });

        return goodLoopsList;
    }

    private static boolean normsAreGood(long start, long width, int resolution, double[] vector1, double[] vector2) {

        int x0 = (int) (start / resolution) - 1;
        int window = (int) (width / resolution) + 3;

        boolean normValuesAreGood = true;

        for (int k = x0; k < x0 + window + 1; k++) {
            if (valueIsBad(vector1[k]) || valueIsBad(vector2[k])) {
                normValuesAreGood = false;
            }
        }

        return normValuesAreGood;
    }

    private static boolean valueIsBad(double v) {
        return Double.isNaN(v) && v < 2;
    }
}
