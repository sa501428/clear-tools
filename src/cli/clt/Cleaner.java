package cli.clt;

import cli.Main;
import cli.utils.clean.LoopTools;
import cli.utils.clean.OracleScorer;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import cli.utils.general.VectorCleaner;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Cleaner {

    public static String usage = "clean <input.hic> <loops.bedpe> <output.bedpe>\n" +
            "clean [--threshold float] <genomeID> <loops.bedpe> <output.bedpe>";

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        double threshold = parser.getThresholdOption(0.5);
        Dataset dataset = null;
        ChromosomeHandler handler;
        if (args[1].endsWith(".hic")) {
            dataset = HiCFileTools.extractDatasetForCLT(args[1], false, true, true);
            handler = dataset.getChromosomeHandler();
        } else {
            handler = ChromosomeTools.loadChromosomes(args[1]);
        }

        String[] bedpeFiles = args[2].split(",");
        String[] outFiles = args[3].split(",");

        if (bedpeFiles.length != outFiles.length) {
            System.err.println("Number of input and output entries don't match");
            System.exit(92);
        }

        for (int z = 0; z < bedpeFiles.length; z++) {
            Feature2DList loopList = LoopTools.loadFilteredBedpe(bedpeFiles[z], handler, true);
            Feature2DList cleanList;
            if (dataset != null) {
                cleanList = cleanupLoops(dataset, loopList, handler);
            } else {
                cleanList = OracleScorer.filter(loopList, threshold);
            }
            cleanList.exportFeatureList(new File(outFiles[z]), false, Feature2DList.ListFormat.NA);
        }
    }

    private static Feature2DList cleanupLoops(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler) {

        Set<HiCZoom> resolutions = getResolutions(loopList);
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

                    Map<Integer, double[]> vectorMap = loadVectors(dataset, chr1, vcNorm, resolutions);
                    try {
                        for (Feature2D loop : loops) {
                            if (normIsOk(loop.getMidPt1(), (int) loop.getWidth1(), vectorMap)
                                    && normIsOk(loop.getMidPt2(), (int) loop.getWidth2(), vectorMap)) {
                                goodLoops.add(loop);
                            }
                        }
                    } catch (Exception e) {
                        System.err.println(e.getMessage());
                    }

                    synchronized (key) {
                        goodLoopsList.addByKey(Feature2DList.getKey(chr1, chr2), new ArrayList<>(goodLoops));
                    }
                }
                threadPair = currChromPair.getAndIncrement();
            }
        });

        return goodLoopsList;
    }

    private static Map<Integer, double[]> loadVectors(Dataset dataset, Chromosome chrom, NormalizationType vcNorm,
                                                      Set<HiCZoom> zooms) {
        Map<Integer, double[]> vectorMap = new HashMap<>();
        for (HiCZoom zoom : zooms) {
            double[] vector = dataset.getNormalizationVector(chrom.getIndex(), zoom,
                    vcNorm).getData().getValues().get(0);
            VectorCleaner.inPlaceClean(vector);
            vectorMap.put(zoom.getBinSize(), vector);
        }
        return vectorMap;
    }

    private static Set<HiCZoom> getResolutions(Feature2DList loopList) {
        Set<HiCZoom> zooms = new HashSet<>();
        Set<Integer> resolutions = getResolutionSet(loopList);
        for (Integer res : resolutions) {
            zooms.add(new HiCZoom(res));
        }
        return zooms;
    }

    private static Set<Integer> getResolutionSet(Feature2DList loopList) {
        Set<Integer> resolutions = new HashSet<>();
        loopList.processLists((s, list) -> {
            for (Feature2D feature2D : list) {
                resolutions.add((int) feature2D.getWidth1());
                resolutions.add((int) feature2D.getWidth2());
            }
        });
        return resolutions;
    }

    private static boolean normIsOk(long pos, int resolution, Map<Integer, double[]> vMap) {
        double[] vector = vMap.get(resolution);
        return vector[(int) (pos / resolution)] > 0;
    }
}
