package cli.clt;

import cli.Main;
import cli.utils.sift.ExtremePixels;
import cli.utils.sift.FeatureUtils;
import cli.utils.sift.Region;
import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
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
import java.util.concurrent.atomic.AtomicInteger;


public class Sift {

    private static final int[] resolutions = new int[]{100, 200, 500, 1000, 2000, 5000, 10000};
    private NormalizationType norm = NormalizationHandler.NONE;


    /*
     * Iteration notes
     * 1 - original pass
     * 2 - similar, fix stuff
     * 3 - z-score at 3, diff filtering
     * 4 - z-score t0 2.5
     * 5 - change z-score to 2
     * 6 - change hires from 100 to 200, and add 500 and 2000 to the resolution runs (already have 1000, 5000 )
     * 7 - remove 500 bp res
     * 8 - coalesce pixels radius 13kb
     * 9 - coalesce pixels radius 5kb
     * 10 - use multiple normalizations *** seems best so far, but pretty conservative
     * 11 - only validate at 5000, not 1000/2000
     * 12 - higher z-score cutoff for 5kb
     * 13 - just the hires calls, collapse at 2k
     * 14 - scale/vc vector filtering, collapse at 5k
     * 15 - restore global max filtering
     * 16 - linear distance for 5k low-res expected
     * 17 - incorporate local enrichment; only use 1 norm for global enrichment at low-res,
     *      remove vector filtering but keep denom filtering
     * 18 - no max neighbor
     * 19 - change to 20% enrichment, not 25%; also put pre-filtering earlier
     * 20 - change low res z-score to 1
     * 21 - change low res z-score to 1.5
     * 22 - same as 19 (low res z-
     * score 2)
     * 23 - 2x min local enrichment
     *
     * 30 - try with 1k and 5k, with pc filter
     * 31 - restrict 1k partially
     * 32 - remove 1k
     * 33 - simplify local filtering (medium status)
     * 34 - add back 1kb
     * 40 - adjust neighbor windows to 1 and 2 (good)
     * 50 - restructure (worse??)
     * 51 - include 2kb, 10kb
     * 52 - higher cutoff
     *
     * 55 - try stepwise filtering, full range
     * 56 - use 1 as min for z score / extreme values
     * 57 - norm vector 0.7
     * 58 - norm vector 0.7 -> 1; low z score 2 -> 1.5
     * 60 - parallelize
     */

    public Sift(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            Main.printGeneralUsageAndExit(5);
        }

        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);

        String normString = parser.getNormalizationStringOption();
        if (normString != null && normString.length() > 1) {
            norm = ds.getNormalizationHandler().getNormTypeFromString(normString);
        }
        System.out.println("Using norm: " + norm.getLabel());

        Feature2DList refinedLoops = siftThroughCalls(ds);
        refinedLoops.exportFeatureList(new File(args[2] + ".sift.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("sift complete");
    }

    private Feature2DList siftThroughCalls(Dataset ds) {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList output = new Feature2DList();

        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();

        AtomicInteger cIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chromosome = chromosomes[currIndex];
                List<Feature2D> sharpLoops = findSiftedFeatures(chromosome, ds);
                if (sharpLoops.size() > 0) {
                    synchronized (output) {
                        output.addByKey(Feature2DList.getKey(chromosome, chromosome), sharpLoops);
                    }
                }
                currIndex = cIndex.getAndIncrement();
            }
        });

        return output;
    }

    private List<Feature2D> findSiftedFeatures(Chromosome chromosome, Dataset ds) {
        Matrix matrix = ds.getMatrix(chromosome, chromosome);
        if (matrix == null) return new ArrayList<>();

        final Map<Integer, Set<SimpleLocation>> pixelsForResolutions = new HashMap<>();

        AtomicInteger rIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currResIndex = rIndex.getAndIncrement();
            while (currResIndex < resolutions.length) {
                int lowRes = resolutions[currResIndex];

                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(lowRes));
                if (zd != null) {
                    Set<SimpleLocation> points = ExtremePixels.getExtremePixelsForResolution(ds, zd,
                            chromosome, lowRes, norm);
                    matrix.clearCacheForZoom(new HiCZoom(lowRes));

                    synchronized (pixelsForResolutions) {
                        pixelsForResolutions.put(lowRes, points);
                        System.out.println(lowRes + " completed (" + points.size() + ")");
                    }
                }
                currResIndex = rIndex.getAndIncrement();
            }
        });

        matrix.clearCache();

        Map<Region, Integer> countsForRecord = FeatureUtils.addPointsToCountMap(pixelsForResolutions,
                resolutions);
        pixelsForResolutions.clear();
        Set<Region> finalPoints = FeatureUtils.getPointsWithMoreThan(countsForRecord, 2);
        countsForRecord.clear();
        return FeatureUtils.convertToFeature2Ds(finalPoints, chromosome);
    }
}
