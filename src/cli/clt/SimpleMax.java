package cli.clt;

import cli.utils.FeatureStats;
import cli.utils.flags.Utils;
import cli.utils.general.FusionTools;
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
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.io.File;
import java.util.List;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class SimpleMax {

    public static void printUsageAndExit() {
        System.out.println("simple-max [-r resolution] [-k norm] <file.hic> <loops.bedpe> <output.bedpe>");
        System.exit(7);
    }

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        String inputBedpeFile = args[2];
        String outputPath = args[3];
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);

        int resolution = parser.getResolutionOption(1000);
        if (resolution < 1) {
            resolution = 1000;
        }

        NormalizationType norm = NormalizationHandler.NONE;
        String possibleNorm = parser.getNormalizationStringOption();
        if (possibleNorm != null && possibleNorm.length() > 0) {
            try {
                norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
            } catch (Exception e) {
                norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{possibleNorm, "SCALE", "KR", "NONE"});
            }
        }
        System.out.println("Normalization being used: " + norm.getLabel());

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpeFile, handler,
                false, null, false);

        Feature2DList localized = findSimpleMax(ds, resolution, loopList, norm);
        localized.exportFeatureList(new File(outputPath), false, Feature2DList.ListFormat.NA);
        System.out.println("simple-max complete");
    }

    private static Feature2DList findSimpleMax(Dataset ds, int resolution, Feature2DList loopList,
                                               NormalizationType norm) {

        Feature2DList finalLoopList = new Feature2DList();

        AtomicInteger currChromIndex = new AtomicInteger(0);
        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        ParallelizationTools.launchParallelizedCode(() -> {

            int threadIndex = currChromIndex.getAndIncrement();
            while (threadIndex < chromosomes.length) {

                Chromosome chrom1 = chromosomes[threadIndex];

                List<Feature2D> loops = loopList.get(chrom1.getIndex(), chrom1.getIndex());
                if (loops != null && loops.size() > 0) {
                    Matrix matrix = ds.getMatrix(chrom1, chrom1, resolution);
                    if (matrix != null) {
                        // load matrix in "resolution" resolution of that chromosome in sparse format
                        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                        if (zd != null) { // if it's not empty
                            try {
                                Set<Feature2D> newLoops = new HashSet<>();
                                Collection<LinkedList<Feature2D>> loopGroups = FusionTools.groupNearbyRecords(
                                        new HashSet<>(loops), 500000).values();
                                for (LinkedList<Feature2D> group : loopGroups) {
                                    long minR = (FeatureStats.minStart1(group) / resolution) - 1;
                                    long minC = (FeatureStats.minStart2(group) / resolution) - 1;
                                    long maxR = (FeatureStats.maxEnd1(group) / resolution) + 1;
                                    long maxC = (FeatureStats.maxEnd2(group) / resolution) + 1;
                                    float[][] regionMatrix = Utils.getRegion(zd, minR, minC, maxR, maxC, norm);
                                    for (Feature2D loop : group) {
                                        getTheMaxPixel(regionMatrix, loop, resolution, newLoops, minR, minC);
                                    }
                                    regionMatrix = null;
                                }

                                synchronized (finalLoopList) {
                                    finalLoopList.addByKey(Feature2DList.getKey(chrom1, chrom1), new ArrayList<>(newLoops));
                                }

                            } catch (Exception e) {
                                e.printStackTrace();
                                System.exit(76);
                            }
                        }
                        matrix.clearCache();
                    }
                }

                threadIndex = currChromIndex.getAndIncrement();
            }
        });

        return finalLoopList;
    }

    private static void getTheMaxPixel(float[][] regionMatrix, Feature2D loop, int resolution,
                                       Set<Feature2D> newLoops, long minR, long minC) {

        int r0 = (int) ((loop.getStart1() / resolution) - minR);
        int c0 = (int) ((loop.getStart2() / resolution) - minC);
        int rF = (int) ((loop.getEnd1() / resolution) + 1 - minR);
        int cF = (int) ((loop.getEnd2() / resolution) + 1 - minC);

        int[] binCoords = getMaxPixelInBox(regionMatrix, r0, rF, c0, cF);

        // arbitrarily chose binCoordsNONE, since NONE and VC have same coordinates
        long startX = (binCoords[0] + minR) * resolution;
        long endX = startX + resolution;
        long startY = (binCoords[1] + minC) * resolution;
        long endY = startY + resolution;
        Feature2D feature = new Feature2D(loop.getFeatureType(), loop.getChr1(), startX, endX,
                loop.getChr2(), startY, endY, Color.BLACK, loop.getAttributes());
        newLoops.add(feature);
    }

    private static int[] getMaxPixelInBox(float[][] matrix, int r0, int rF, int c0, int cF) {
        float maxCounts = matrix[r0][c0];
        int[] coordinates = new int[]{r0, c0};
        for (int r = r0; r < rF; r++) {
            for (int c = c0; c < cF; c++) {
                if (matrix[r][c] > maxCounts) { // will skip NaNs
                    maxCounts = matrix[r][c];
                    coordinates[0] = r;
                    coordinates[1] = c;
                }
            }
        }
        return coordinates;
    }
}
