package cli.clt;

import cli.Main;
import cli.utils.FeatureStats;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import cli.utils.general.QuickGrouping;
import cli.utils.general.Utils;
import cli.utils.pinpoint.ConvolutionTools;
import cli.utils.pinpoint.LandScape;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Pinpoint {
    private static final NormalizationType NONE = NormalizationHandler.NONE;
    public static String usage = "pinpoint [--only-one] [--res int] <input.hic> <loops.bedpe> <output>";

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);
        String loopListPath = args[2];
        String outFile = args[3];

        ChromosomeHandler handler = dataset.getChromosomeHandler();

        boolean onlyGetOne = parser.getOnlyOneOption();

        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                true, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        int resolution = parser.getResolutionOption(-1);
        if (resolution < 1) {
            resolution = HiCUtils.getHighestResolution(dataset.getBpZooms()).getBinSize();
        }

        final Feature2DList pinpointLoopsNoNorm = new Feature2DList();
        final Feature2DList pinpointBounds = new Feature2DList();

        localize(dataset, loopList, handler, resolution, onlyGetOne, pinpointLoopsNoNorm, pinpointBounds);
        pinpointLoopsNoNorm.exportFeatureList(new File(outFile + "_raw.bedpe"), false, Feature2DList.ListFormat.NA);
        pinpointBounds.exportFeatureList(new File(outFile + "_bounds.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("pinpoint complete");
    }

    private static void localize(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler,
                                 int resolution, boolean onlyGetOne,
                                 Feature2DList finalLoops, Feature2DList finalBounds) {

        if (Main.printVerboseComments) {
            System.out.println("Pinpointing location for loops");
        }

        int compressedKernelSize = Math.max((400 / resolution) + 1, 3); // effectively multiplied by 20 later
        int kernelSize = Math.max((2000 / resolution) + 1, 3);

        final int[] globalMaxWidth = new int[1];
        loopList.processLists((s, list) -> {
            int maxWidth = 0;
            for (Feature2D feature : list) {
                maxWidth = (int) Math.max(maxWidth, feature.getWidth1() / resolution);
                maxWidth = (int) Math.max(maxWidth, feature.getWidth2() / resolution);
            }
            globalMaxWidth[0] = Math.max(globalMaxWidth[0], maxWidth);
        });

        int halfMatrixWidth = (int) (1.5 * Math.max(globalMaxWidth[0], 10)) + 1;
        int matrixWidth = 2 * halfMatrixWidth + 1;

        //GPUController gpuController = Circe.buildGPUController(kernelSize, matrixWidth, kernelSize / 2 + 1);
        final float[][] kernel = ConvolutionTools.getManhattanKernel(kernelSize);
        final float[][] compressedKernel = ConvolutionTools.getManhattanKernel(compressedKernelSize);

        HiCZoom zoom = new HiCZoom(resolution);

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), true);

        int numTotalLoops = loopList.getNumTotalFeatures();
        final AtomicInteger currChromPair = new AtomicInteger(0);
        final AtomicInteger currNumLoops = new AtomicInteger(0);

        final Object key = new Object();


        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                Matrix matrix = dataset.getMatrix(chr1, chr2);
                if (matrix != null) {
                    List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                    if (loops != null && loops.size() > 0) {
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) {
                            try {
                                Collection<List<Feature2D>> loopGroups = QuickGrouping.groupNearbyRecords(
                                        loops, matrixWidth * resolution).values();

                                List<Feature2D> pinpointedLoops = new ArrayList<>();
                                List<Feature2D> pinpointedBounds = new ArrayList<>();

                                for (List<Feature2D> group : loopGroups) {
                                    Feature2D loop = group.get(0);
                                    String saveString = generateLoopInfo(loop);
                                    int minR = (int) ((FeatureStats.minStart1(group) / resolution) - halfMatrixWidth);
                                    int minC = (int) ((FeatureStats.minStart2(group) / resolution) - halfMatrixWidth);
                                    int maxR = (int) ((FeatureStats.maxEnd1(group) / resolution) + halfMatrixWidth);
                                    int maxC = (int) ((FeatureStats.maxEnd2(group) / resolution) + halfMatrixWidth);

                                    List<ContactRecord> records = Utils.getRecords(zd, minR, minC, maxR, maxC, NONE);

                                    LandScape.extractMaxima(records, minR, minC, resolution,
                                            pinpointedLoops, pinpointedBounds, loop, saveString,
                                            onlyGetOne, matrixWidth, kernel, compressedKernel);
                                    if (currNumLoops.addAndGet(group.size()) % 100 == 0) {
                                        System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                    }
                                }

                                synchronized (key) {
                                    finalLoops.addByKey(Feature2DList.getKey(chr1, chr2), pinpointedLoops);
                                    finalBounds.addByKey(Feature2DList.getKey(chr1, chr2), pinpointedBounds);
                                }

                                System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                            } catch (Exception e) {
                                System.err.println(e.getMessage());
                            }
                        }
                    }
                    matrix.clearCache();
                }
                threadPair = currChromPair.getAndIncrement();
            }
        });
    }

    private static String generateLoopInfo(Feature2D loop) {
        String info = loop.simpleString();
        String[] saveStrings = info.split("\\s+");
        info = String.join("_", saveStrings);
        return info;
    }
}
