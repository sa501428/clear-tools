package cli.clt;

import cli.Main;
import cli.utils.expected.ExpectedModel;
import cli.utils.expected.LogBinnedExpectedModel;
import cli.utils.expected.LogExpectedSpline;
import cli.utils.flags.RegionConfiguration;
import cli.utils.flags.Utils;
import cli.utils.general.HiCUtils;
import cli.utils.recap.RecapTools;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Recap {

    private NormalizationType norm;

    public Recap(String[] args, CommandLineParser parser) {

        // if input doesn't follow required Recap syntax/format, then exit
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5);
        }

        // recap <loops.bedpe> <outfolder> <file1.hic,file2.hic,...> <name1,name2,...>
        // makeDir creates a new directory if there isn't one already, stores directory in outFolder var

        // recap <loops.bedpe> <outfolder> <file1.hic,file2.hic,...> <name1,name2,...>
        String loopListPath = args[1];
        File outFolder = UNIXTools.makeDir(new File(args[2]));

        // filepaths = [file1.hic, file2.hic, ...]
        // names = [name1, name2, ...]
        String[] filepaths = args[3].split(",");
        String[] names = args[4].split(",");
        
        Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[0], false, false, true);


        // grabs chromosomes...
        // todo: what is this Feature2DList object? I think it is storing arg[0] (bedpe file) as an object
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                false, null, false);

        // get a norm and print norm
        String possibleNorm = parser.getNormalizationStringOption();
        try {
            if (possibleNorm != null && possibleNorm.length() > 0) {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            } else {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR"});
            }
        } catch (Exception e) {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "VC_SQRT", "VC"});
        }
        System.out.println("Using normalization: " + norm.getLabel());

        ds.clearCache(false);
        ds = null;

        // I'm guessing getNumTotalFeatures gets the number of loops in the looplist?
        int numRows = loopList.getNumTotalFeatures();
        System.out.println("Number of loops: " + numRows);

        // since parser stores the options and their values from when you input arguments in the command line
        // parser.getResolutionOption takes the option stored
        // stores it in int resolution
        // by default, 1kb res
        int resolution = parser.getResolutionOption(1000);
        if (resolution < 1) {
            resolution = 1000;
        }
        System.out.println("Using resolution: " + resolution);
        boolean isDeepLoopAnalysis = parser.getIsLoopAnalysis();
        
        int window = parser.getWindowSizeOption(0);
        if (window < 1) {
            if (isDeepLoopAnalysis) {
                window = 10;
            } else {
                window = 1;
            }
        }
        
        Feature2DList refinedLoops = recapStats(filepaths, names, loopList, handler, resolution, window, norm, isDeepLoopAnalysis);
        refinedLoops.exportFeatureList(new File(outFolder, "recap.bedpe"), false, Feature2DList.ListFormat.NA);
        RecapTools.exportAllMatrices(handler.getChromosomeArrayWithoutAllByAll(), refinedLoops, names, outFolder, isDeepLoopAnalysis, window);
        System.out.println("recap complete");
    }

    private static Feature2DList recapStats(String[] filepaths, String[] names, Feature2DList loopList,
                                            ChromosomeHandler handler, int resolution, int window,
                                            NormalizationType norm, boolean isDeepLoopAnalysis) {

        // VerboseComments is a boolean that, if true, will have extra print statements describing task progress
        // probably for debugging
        if (Main.printVerboseComments) {
            System.out.println("Start Recap/Compile process");
        }

        // todo: ask Muhammad about this... what exactly is the zoom object
        HiCZoom zoom = new HiCZoom(resolution);
        int matrixWidth = 1 + 2 * window;

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), false);

        int numTotalLoops = loopList.getNumTotalFeatures();
        final Object key = new Object();

        for (int di = 0; di < filepaths.length; di++) {

            System.out.println("File (" + (di + 1) + "/" + filepaths.length + "): " + filepaths[di]);

            final Dataset ds = HiCFileTools.extractDatasetForCLT(filepaths[di], false, false,
                    true);
            String prefix = names[di] + "_";

            final AtomicInteger currChromPair = new AtomicInteger(0);
            final AtomicInteger currNumLoops = new AtomicInteger(0);

            ParallelizationTools.launchParallelizedCode(() -> {

                int threadPair = currChromPair.getAndIncrement();
                while (threadPair < chromosomePairCounter) {
                    RegionConfiguration config = chromosomePairs.get(threadPair);
                    Chromosome chrom1 = config.getChr1();
                    Chromosome chrom2 = config.getChr2();
                    
                    List<Feature2D> loops = loopList.get(chrom1.getIndex(), chrom2.getIndex());
                    if (loops != null && loops.size() > 0) {
                        Matrix matrix = ds.getMatrix(chrom1, chrom2, resolution);
                        if (matrix == null) {
                            System.err.println("Matrix is null " + chrom1.getName() + "_" + chrom2.getName());
                            System.exit(9);
                        }

                        MatrixZoomData zd = matrix.getZoomData(zoom);

                        if (zd == null) {
                            System.err.println("ZD is null " + chrom1.getName() + "_" + chrom2.getName());
                            System.exit(9);
                        }

                        int maxBinDist = Math.max(getMaxDistance(loops, resolution, window), 9000000 / resolution);
                        LogBinnedExpectedModel expected = new LogBinnedExpectedModel(zd, norm, maxBinDist, 0);
                        LogExpectedSpline spline = expected.getSpline();

                        float pseudoCount = getMedianExpectedAt(maxBinDist - 2 * window, expected);
                        double superDiagonal = expected.getExpectedFromUncompressedBin(1);

                        try {
                            for (Feature2D loop : loops) {
                                float[][] obsMatrix = new float[matrixWidth][matrixWidth];
                                Utils.addLocalizedData(obsMatrix, zd, loop, matrixWidth, resolution, window, norm, key);
                                // MatrixTools.saveMatrixTextNumpy((new File(outFolder, saveString + "_raw.npy")).getAbsolutePath(), output);

                                Map<String, String> attributes = RecapTools.getStats(obsMatrix,
                                        window, superDiagonal, pseudoCount, isDeepLoopAnalysis, spline,
                                        expected, loop, resolution);
                                for (String akey : attributes.keySet()) {
                                    loop.addStringAttribute(prefix + akey, attributes.get(akey));
                                }

                                if (currNumLoops.incrementAndGet() % 1000 == 0) {
                                    System.out.print(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");
                                }
                            }

                            System.out.println(((int) Math.floor((100.0 * currNumLoops.get()) / numTotalLoops)) + "% ");

                        } catch (Exception e) {
                            e.printStackTrace();
                            System.exit(76);
                        }
                        matrix.clearCache();
                    }
                    threadPair = currChromPair.getAndIncrement();
                }
            });
            ds.clearCache(false);
        }

        return loopList;
    }

    private static int getMaxDistance(List<Feature2D> loops, int resolution, int window) {
        long maxDist = 0;
        for (Feature2D loop : loops) {
            long dist = Math.abs((loop.getStart1() / resolution - window) - (loop.getEnd2() / resolution + window));
            if (dist > maxDist) {
                maxDist = dist;
            }
        }
        return (int) (maxDist + 4 * window);
    }

    private static float getMedianExpectedAt(int d0, ExpectedModel expectedVector) {
        return (float) expectedVector.getExpectedFromUncompressedBin(d0);
    }
}
