package cli.clt;

import cli.Main;
import cli.apa.APAUtils;
import cli.apa.RegionConfiguration;
import cli.utils.ConvolutionTools;
import cli.utils.HiCUtils;
import cli.utils.cc.ConnectedComponents;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Pinpoint {
    public static void run(String[] args) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[1], false, true);
        String bedpeFile = args[2];
        String outFolder = args[3];

        ChromosomeHandler handler = dataset.getChromosomeHandler();

        Feature2DList loopList = Feature2DParser.loadFeatures(bedpeFile, handler,
                false, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        UNIXTools.makeDir(outFolder);

        amplifyLoops(dataset, loopList, handler, outFolder);
    }

    private static void amplifyLoops(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler,
                                     String outFolder) {

        if (Main.printVerboseComments) {
            System.out.println("Pinpointing location for loops");
        }

        HiCZoom zoom = HiCUtils.getHighestResolution(dataset.getBpZooms());
        int resolution = zoom.getBinSize();

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getAutosomalChromosomesArray());

        final AtomicInteger currChromPair = new AtomicInteger(0);
        final Object key = new Object();
        final Feature2DList refinedLoops = new Feature2DList();

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                if (loops.size() > 0) {
                    List<Feature2D> pinpointedLoops = new ArrayList<>();
                    MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chr1, chr2, zoom);

                    if (zd != null) {
                        try {
                            for (Feature2D loop : loops) {

                                int window = (int) (Math.max(loop.getWidth1(), loop.getWidth2()) / resolution + 1);

                                int binXStart = (int) ((loop.getStart1() / resolution) - window);
                                int binYStart = (int) ((loop.getStart2() / resolution) - window);

                                int matrixWidth = 3 * window + 1;
                                int[][] output = new int[matrixWidth][matrixWidth];
                                APAUtils.addRawLocalBoundedRegion(output, zd,
                                        binXStart, binYStart, window, matrixWidth, key);

                                String saveString = loop.simpleString();
                                String[] saveStrings = saveString.split("\\s+");
                                saveString = String.join("_", saveStrings);

                                //MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_raw.npy")).getAbsolutePath(), output);
                                float[][] kde = ConvolutionTools.sparseConvolution(output);
                                output = null; // clear output
                                //MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_kde.npy")).getAbsolutePath(), kde);

                                ConnectedComponents.extractMaxima(kde, binXStart, binYStart, resolution,
                                        pinpointedLoops, loop, outFolder, saveString);

                                kde = null;
                            }
                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                        }

                        synchronized (key) {
                            refinedLoops.addByKey(Feature2DList.getKey(chr1, chr2), pinpointedLoops);
                        }
                    }
                }
                //System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus.get())) + "% ");
                threadPair = currChromPair.getAndIncrement();
            }
        });
        System.out.println("pinpoint complete");
    }
}
