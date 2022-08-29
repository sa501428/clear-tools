package cli.clt;

import cli.Main;
import cli.utils.flags.RegionConfiguration;
import cli.utils.flags.Utils;
import cli.utils.general.HiCUtils;
import cli.utils.general.WritingTools;
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
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Enhance {
    public static void run(String[] args, int resolution, boolean exportNPY) {
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
            datasets[q] = HiCFileTools.extractDatasetForCLT(hicFiles[q], false, true, true);
            keys[q] = new Object();
        }

        ChromosomeHandler handler = datasets[0].getChromosomeHandler();

        Feature2DList loopList = Feature2DParser.loadFeatures(bedpeFile, handler, false, null, false);

        System.out.println("Number of loops: " + loopList.getNumTotalFeatures());

        UNIXTools.makeDir(outFolder);

        amplifyLoops(datasets, loopList, handler, resolution, outFolder, keys, exportNPY);
    }

    private static void amplifyLoops(final Dataset[] datasets, Feature2DList loopList, ChromosomeHandler handler,
                                     int resolution, String outFolder, Object[] keys,
                                     boolean exportNPY) {

        if (Main.printVerboseComments) {
            System.out.println("Processing AMPLIFI for resolution " + resolution);
        }
        HiCZoom zoom = new HiCZoom(HiCZoom.HiCUnit.BP, resolution);

        final AtomicInteger currentProgressStatus = new AtomicInteger(0);

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        int pairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), true);

        final int chromosomePairCounter = pairCounter;
        final AtomicInteger maxProgressStatus = new AtomicInteger(pairCounter);
        final AtomicInteger chromosomePair = new AtomicInteger(0);

        final Object mndKey = new Object();
        List<String> filePaths = Collections.synchronizedList(new ArrayList<>());

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = chromosomePair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chr1 = config.getChr1();
                Chromosome chr2 = config.getChr2();

                List<Feature2D> loops = loopList.get(chr1.getIndex(), chr2.getIndex());
                if(loops.size() > 0) {
                    String newMND = (new File(outFolder, config.getPairKey() + ".mnd")).getAbsolutePath();
                    try {
                        BufferedWriter bwMND = new BufferedWriter(new OutputStreamWriter(Files.newOutputStream(Paths.get(newMND))));
                        for (Feature2D loop : loops) {
                            int window = (int) (Math.max(loop.getWidth1(), loop.getWidth2()) / resolution + 1);
                            window = Math.max(window, 10000 / resolution);
                            int binXStart = (int) ((loop.getMidPt1() / resolution) - window);
                            int binYStart = (int) ((loop.getMidPt2() / resolution) - window);

                            int matrixWidth = 2 * window + 1;
                            float[][] output = new float[matrixWidth][matrixWidth];

                            for (int di = 0; di < datasets.length; di++) {
                                final Dataset ds = datasets[di];
                                MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, zoom);
                                if (zd != null) {
                                    try {
                                        Utils.addLocalBoundedRegion(output, zd,
                                                binXStart, binYStart, matrixWidth,
                                                NormalizationHandler.NONE, keys[di]);
                                    } catch (Exception e) {
                                        System.err.println(e.getMessage());
                                        System.err.println("Unable to find data for loop: " + loop);
                                    }
                                }
                            }

                            String saveString = loop.simpleString();
                            String[] saveStrings = saveString.split("\\s+");
                            saveString = String.join("_", saveStrings);

                            if (exportNPY) {
                                MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + ".npy")).getAbsolutePath(),
                                        output);
                            }
                            WritingTools.writeToMND(output, resolution, chr1.getName(), chr2.getName(),
                                    binXStart, binYStart, bwMND);
                            output = null;
                        }
                        bwMND.close();
                    } catch (Exception e){
                        e.printStackTrace();
                        System.exit(9);
                    }
                    synchronized (mndKey) {
                        filePaths.add(newMND);
                    }
                }

                for (final Dataset ds : datasets) {
                    Matrix matrix = ds.getMatrix(chr1, chr2);
                    if (matrix != null) {
                        matrix.clearCacheForZoom(zoom);
                    }
                }

                //System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus.get())) + "% ");
                threadPair = chromosomePair.getAndIncrement();
            }
        });

        System.out.println("MND lists complete");

        String mndPath = WritingTools.buildCatScript(filePaths, outFolder);
        String genomeID = WritingTools.cleanGenome(datasets[0].getGenomeId());
        String newHiCFile = new File(outFolder, "enhance.hic").getAbsolutePath();
        String resolutionsToBuild = WritingTools.getResolutionsToBuild(resolution);
        try {
            System.out.println("Run HiCTools Pre: \n   pre -n -r" + resolutionsToBuild + " " +
                    mndPath + " " + newHiCFile + " " + genomeID);
        } catch (Exception e){
            e.printStackTrace();
            System.exit(31);
        }

        //if no data return null
    }
}
