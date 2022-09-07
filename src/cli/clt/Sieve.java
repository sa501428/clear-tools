package cli.clt;

import cli.Main;
import cli.utils.flags.RegionConfiguration;
import cli.utils.general.HiCUtils;
import javastraw.expected.ExpectedModel;
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
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Sieve {

    private NormalizationType norm = NormalizationHandler.VC;

    /**
     * retain hi-res loops if at a 'loop-y' location per low-res assessment
     */
    public Sieve(String[] args, CommandLineParser parser) {
        // sieve <loops.bedpe> <output.bedpe> <file1.hic> <res1,res2,...>
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5);
        }

        String loopListPath = args[1];
        String outfile = args[2];
        String filepath = args[3];
        int[] resolutions = parseInts(args[4]);

        Dataset ds = HiCFileTools.extractDatasetForCLT(filepath, false, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler,
                false, null, false);

        String possibleNorm = parser.getNormalizationStringOption();
        if (possibleNorm != null && possibleNorm.length() > 0) {
            norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
        }
        System.out.println("Using normalization: " + norm.getLabel());

        int window = parser.getWindowSizeOption(0);
        if (window < 2) {
            window = 2;
        }

        Feature2DList result = retainMaxPeaks(ds, loopList, handler, resolutions, window, norm);
        result.exportFeatureList(new File(outfile), false, Feature2DList.ListFormat.NA);

        System.out.println("sieve complete");
    }

    private static Feature2DList retainMaxPeaks(Dataset ds, Feature2DList loopList,
                                                ChromosomeHandler handler, int[] resolutions, int window,
                                                NormalizationType norm) {

        if (Main.printVerboseComments) {
            System.out.println("Start Recap/Compile process");
        }

        final Feature2DList newLoopList = new Feature2DList();

        Map<Integer, RegionConfiguration> chromosomePairs = new ConcurrentHashMap<>();
        final int chromosomePairCounter = HiCUtils.populateChromosomePairs(chromosomePairs,
                handler.getChromosomeArrayWithoutAllByAll(), false);

        final AtomicInteger currChromPair = new AtomicInteger(0);
        final AtomicInteger currNumLoops = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(() -> {
            int threadPair = currChromPair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                RegionConfiguration config = chromosomePairs.get(threadPair);
                Chromosome chrom1 = config.getChr1();
                Chromosome chrom2 = config.getChr2();

                List<Feature2D> loops = loopList.get(chrom1.getIndex(), chrom2.getIndex());
                loops = new ArrayList<>(new HashSet<>(loops));
                Set<Feature2D> loopsToKeep = new HashSet<>();
                if (loops.size() > 0) {
                    Matrix matrix = ds.getMatrix(chrom1, chrom2);
                    if (matrix != null) {
                        for (int resolution : resolutions) {
                            HiCZoom zoom = new HiCZoom(resolution);
                            MatrixZoomData zd = matrix.getZoomData(zoom);
                            if (zd != null) {
                                try {
                                    // figure out if this loop is present at the resolution given

                                    synchronized (newLoopList) {
                                        newLoopList.addByKey(Feature2DList.getKey(chrom1, chrom2),
                                                new ArrayList<>(loopsToKeep));
                                    }
                                } catch (Exception e) {
                                    e.printStackTrace();
                                    System.exit(76);
                                }
                            }
                            matrix.clearCacheForZoom(zoom);
                        }
                        matrix.clearCache();
                    }

                }
                threadPair = currChromPair.getAndIncrement();
            }
        });


        return newLoopList;
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

    private int[] parseInts(String input) {
        String[] inputs = input.split(",");
        int[] values = new int[inputs.length];
        for (int i = 0; i < values.length; i++) {
            values[i] = Integer.parseInt(inputs[i]);
        }
        Arrays.sort(values);
        return values;
    }
}
