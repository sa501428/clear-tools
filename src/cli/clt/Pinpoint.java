package cli.clt;

import cli.Main;
import cli.apa.APAUtils;
import cli.apa.RegionConfiguration;
import cli.utils.ConvolutionTools;
import cli.utils.HiCUtils;
import cli.utils.cc.ConnectedComponents;
import cli.utils.cc.Location2D;
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

import java.io.File;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Pinpoint {
    public static void run(String[] args, int resolutionOption) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);
        String loopListPath = args[2];
        String outFile = args[3];

        ChromosomeHandler handler = dataset.getChromosomeHandler();

        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                false, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        int resolution = resolutionOption;
        if (resolution < 1) {
            resolution = HiCUtils.getHighestResolution(dataset.getBpZooms()).getBinSize();
        }

        Feature2DList refinedLoops = localize(dataset, loopList, handler, resolution);
        refinedLoops.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
        System.out.println("pinpoint complete");
    }

    private static Feature2DList localize(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler,
                                          int resolution) {

        if (Main.printVerboseComments) {
            System.out.println("Pinpointing location for loops");
        }

        HiCZoom zoom = new HiCZoom(resolution);

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getAutosomalChromosomesArray());

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

                List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                if (loops != null && loops.size() > 0) {
                    MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chr1, chr2, zoom);
                    if (zd != null) {
                        try {
                            Set<Location2D> locations = new HashSet<>();
                            for (Feature2D loop : loops) {

                                int window = (15000 / resolution + 1);

                                int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
                                int binYStart = (int) ((loop.getMidPt2() / resolution) - window);

                                int matrixWidth = 2 * window + 1;
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
                                        locations, loop, saveString);

                                kde = null;

                                if (currNumLoops.incrementAndGet() % 100 == 0) {
                                    System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                }
                            }
                            zd.clearCache();

                            List<Feature2D> pinpointedLoops = convertToFeature2D(locations);
                            locations.clear();
                            synchronized (key) {
                                refinedLoops.addByKey(Feature2DList.getKey(chr1, chr2), pinpointedLoops);
                            }
                            pinpointedLoops.clear();

                            System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                        }
                    }
                }
                threadPair = currChromPair.getAndIncrement();
            }
        });
        return refinedLoops;
    }

    private static List<Feature2D> convertToFeature2D(Set<Location2D> locations) {
        List<Feature2D> features = new ArrayList<>();
        for (Location2D location : locations) {
            features.add(location.toFeature2D());
        }
        return features;
    }
}
