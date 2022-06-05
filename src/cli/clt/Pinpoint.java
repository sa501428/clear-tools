package cli.clt;

import cli.Main;
import cli.apa.APAUtils;
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
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class Pinpoint {
    public static void run(String[] args) {
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(5);
        }

        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[1], false, true, true);
        String loopListPath = args[2];
        String outFolder = args[3];

        ChromosomeHandler handler = dataset.getChromosomeHandler();

        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                false, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        UNIXTools.makeDir(outFolder);

        Feature2DList refinedLoops = localize(dataset, loopList, handler, outFolder);
        refinedLoops.exportFeatureList(new File(outFolder, "pinpoint.bedpe"),
                false, Feature2DList.ListFormat.NA);
        System.out.println("pinpoint complete");
    }

    private static Feature2DList localize(final Dataset dataset, Feature2DList loopList, ChromosomeHandler handler,
                                          String outFolder) {

        if (Main.printVerboseComments) {
            System.out.println("Pinpointing location for loops");
        }

        HiCZoom zoom = HiCUtils.getHighestResolution(dataset.getBpZooms());
        int resolution = zoom.getBinSize();

        int numTotalLoops = loopList.getNumTotalFeatures();

        final AtomicInteger currNumLoops = new AtomicInteger(0);

        final Feature2DList refinedLoops = new Feature2DList();
        final Object key = new Object();

        long interval1 = 0L, interval2 = 0L, interval3 = 0L;

        for (Chromosome chromosome : handler.getChromosomeArrayWithoutAllByAll()) {
            List<Feature2D> loops = loopList.get(chromosome.getIndex(), chromosome.getIndex());
            if (loops != null && loops.size() > 0) {
                MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chromosome, chromosome, zoom);
                if (zd != null) {
                    try {
                        List<Feature2D> pinpointedLoops = new ArrayList<>();
                        for (Feature2D loop : loops) {

                            int window = (int) (Math.max(loop.getWidth1(), loop.getWidth2()) / resolution + 1);

                            int binXStart = (int) ((loop.getStart1() / resolution) - window);
                            int binYStart = (int) ((loop.getStart2() / resolution) - window);

                            long time0 = System.nanoTime();
                            int matrixWidth = 3 * window + 1;
                            int[][] output = new int[matrixWidth][matrixWidth];
                            APAUtils.addRawLocalBoundedRegion(output, zd,
                                    binXStart, binYStart, window, matrixWidth, key);
                            long time1 = System.nanoTime();

                            String saveString = loop.simpleString();
                            String[] saveStrings = saveString.split("\\s+");
                            saveString = String.join("_", saveStrings);

                            //MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_raw.npy")).getAbsolutePath(), output);

                            long time2 = System.nanoTime();
                            float[][] kde = ConvolutionTools.sparseConvolution(output);
                            long time3 = System.nanoTime();
                            output = null; // clear output
                            //MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_kde.npy")).getAbsolutePath(), kde);

                            long time4 = System.nanoTime();
                            ConnectedComponents.extractMaxima(kde, binXStart, binYStart, resolution,
                                    pinpointedLoops, loop, outFolder, saveString);
                            long time5 = System.nanoTime();

                            kde = null;

                            interval1 += (time1 - time0);
                            interval2 += (time3 - time2);
                            interval3 += (time5 - time4);

                            if (currNumLoops.incrementAndGet() % 5 == 0) {
                                System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                System.out.println("Times " + interval1 + "  " + interval2 + "  " + interval3);
                            }
                        }
                        zd.clearCache();

                        synchronized (key) {
                            refinedLoops.addByKey(Feature2DList.getKey(chromosome, chromosome), pinpointedLoops);
                        }

                        System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                    } catch (Exception e) {
                        System.err.println(e.getMessage());
                    }
                }
            }
        }
        return refinedLoops;
    }
}
