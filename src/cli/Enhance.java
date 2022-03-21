package cli;

import cli.apa.*;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Enhance {
    public static void run(String[] args, int resolution, String normString) {
        if(args.length < 4){
            Main.printGeneralUsageAndExit(5);
        }

        String outFolder = args[1];
        String bedpeFile = args[2];
        String[] hicFiles = new String[args.length - 3];
        System.arraycopy(args, 3, hicFiles, 0, hicFiles.length);

        final Object[] keys = new Object[hicFiles.length];
        Dataset[] datasets = new Dataset[hicFiles.length];
        for(int q = 0; q < hicFiles.length; q++){
            datasets[q] = HiCFileTools.extractDatasetForCLT(hicFiles[q] , false, true);
            keys[q] = new Object();
        }

        ChromosomeHandler handler = datasets[0].getChromosomeHandler();
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(datasets[0],
                new String[]{normString, "KR", "SCALE", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        Feature2DList loopList = Feature2DParser.loadFeatures(bedpeFile, handler, false, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        UNIXTools.makeDir(outFolder);

        amplifyLoops(datasets, loopList, handler, resolution, norm, outFolder, keys);
    }

    private static void amplifyLoops(final Dataset[] datasets, Feature2DList loopList, ChromosomeHandler handler,
                                     int resolution, NormalizationType norm, String outFolder, Object[] keys) {

        if(Main.printVerboseComments) {
            System.out.println("Processing AMPLIFI for resolution " + resolution);
        }
        HiCZoom zoom = new HiCZoom(HiCZoom.HiCUnit.BP, resolution);

        final AtomicInteger currentProgressStatus = new AtomicInteger(0);

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        int pairCounter = 0;
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                RegionConfiguration config = new RegionConfiguration(chromosomes[i], chromosomes[j]);
                chromosomePairs.put(pairCounter, config);
                pairCounter++;
            }
        }
        final int chromosomePairCounter = pairCounter;
        final AtomicInteger maxProgressStatus = new AtomicInteger(pairCounter);
        final AtomicInteger chromosomePair = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = chromosomePair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                if(loops.size() > 0) {
                    for (int li = 0; li < loops.size(); li++) {
                        Feature2D loop = loops.get(li);
                        int window = (int) (Math.max(loop.getWidth1(), loop.getWidth2()) / resolution + 1);
                        int matrixWidth = 2 * window + 1;
                        double[][] output = new double[matrixWidth][matrixWidth];

                        for (int di = 0; di < datasets.length; di++) {
                            final Dataset ds = datasets[di];
                            MatrixZoomData zd;
                            synchronized (ds) {
                                zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, zoom);
                            }

                            if (zd != null) {
                                try {
                                    APAUtils.addLocalizedData(output, zd, loop, matrixWidth, resolution, window,
                                            norm, keys[di]);
                                } catch (Exception e) {
                                    System.err.println(e.getMessage());
                                    System.err.println("Unable to find data for loop: " + loop);
                                }
                            }
                        }

                        String saveString = loop.simpleString();
                        String[] saveStrings = saveString.split("\\s+");
                        saveString = String.join("_", saveStrings);

                        MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + ".npy")).getAbsolutePath(),
                                output);
                        /*
                        MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_agg_norm.npy")).getAbsolutePath(),
                                AggNorm.getAggNormedMatrix(output));
                        */
                    }
                }

                for (int di = 0; di < datasets.length; di++) {
                    final Dataset ds = datasets[di];
                    MatrixZoomData zd;
                    synchronized (ds) {
                        zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, zoom);
                    }

                    if (zd != null) {
                        zd.clearCache();
                    }
                }

                //System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus.get())) + "% ");
                threadPair = chromosomePair.getAndIncrement();
            }
        });

        System.out.println("Amplifi complete");
        //if no data return null
    }
}
