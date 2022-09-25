package cli.clt;

import cli.Main;
import cli.utils.cc.ConnectedComponents;
import cli.utils.flags.RegionConfiguration;
import cli.utils.flags.Utils;
import cli.utils.general.HiCUtils;
import cli.utils.pinpoint.ConvolutionTools;
import cli.utils.pinpoint.LocalNorms;
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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Pinpoint {
    private static final NormalizationType NONE = NormalizationHandler.NONE;
    public static String usage = "pinpoint [--res int] <input.hic> <loops.bedpe> <output.bedpe>";

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
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

        Feature2DList refinedLoops = localize(dataset, loopList, handler, resolution, onlyGetOne);

        String originalLoops = outFile.replace(".bedpe", "");
        originalLoops += "_with_original.bedpe";

        loopList.exportFeatureList(new File(originalLoops), false, Feature2DList.ListFormat.NA);

        refinedLoops.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
        System.out.println("pinpoint complete");
    }

    private static Feature2DList localize(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler,
                                          int resolution, boolean onlyGetOne) {

        if (Main.printVerboseComments) {
            System.out.println("Pinpointing location for loops");
        }

        int kernelSize = Math.max(10, 200 / resolution);

        final int[] globalMaxWidth = new int[1];
        loopList.processLists((s, list) -> {
            int maxWidth = 0;
            for (Feature2D feature : list) {
                maxWidth = (int) Math.max(maxWidth, feature.getWidth1() / resolution);
                maxWidth = (int) Math.max(maxWidth, feature.getWidth2() / resolution);
            }
            globalMaxWidth[0] = Math.max(globalMaxWidth[0], maxWidth);
        });

        int window = Math.max(globalMaxWidth[0], 10);
        int matrixWidth = 3 * window + 1;

        //GPUController gpuController = Circe.buildGPUController(kernelSize, matrixWidth, kernelSize / 2 + 1);

        final float[][] kernel = ConvolutionTools.getManhattanKernel(kernelSize);
        //float maxK = ArrayTools.getMax(kernel);

        HiCZoom zoom = new HiCZoom(resolution);

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), true);

        int numTotalLoops = loopList.getNumTotalFeatures();
        final AtomicInteger currChromPair = new AtomicInteger(0);
        final AtomicInteger currNumLoops = new AtomicInteger(0);

        final Object key = new Object();
        final Feature2DList refinedLoops = new Feature2DList();

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
                                List<Feature2D> pinpointedLoops = new ArrayList<>();
                                for (Feature2D loop : loops) {

                                    int binXStart = (int) ((loop.getStart1() / resolution) - window);
                                    int binYStart = (int) ((loop.getStart2() / resolution) - window);

                                    float[][] output = new float[matrixWidth][matrixWidth];
                                    Utils.addLocalBoundedRegion(output, zd, binXStart, binYStart, matrixWidth, NONE);

                                    String saveString = loop.simpleString();
                                    String[] saveStrings = saveString.split("\\s+");
                                    saveString = String.join("_", saveStrings);

                                    float[][] kde = ConvolutionTools.sparseConvolution(output, kernel);
                                    output = null; // clear output
                                    LocalNorms.normalizeLocally(kde);
                                    ConnectedComponents.extractMaxima(kde, binXStart, binYStart, resolution,
                                            pinpointedLoops, loop, saveString, onlyGetOne);
                                    kde = null;

                                    if (currNumLoops.incrementAndGet() % 100 == 0) {
                                        System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                    }
                                }

                                synchronized (key) {
                                    refinedLoops.addByKey(Feature2DList.getKey(chr1, chr2), pinpointedLoops);
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
        return refinedLoops;
    }
}
