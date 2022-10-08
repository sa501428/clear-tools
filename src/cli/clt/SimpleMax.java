package cli.clt;

import cli.utils.FeatureStats;
import cli.utils.general.QuickGrouping;
import cli.utils.general.Utils;
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

    private static final NormalizationType VC = NormalizationHandler.VC;
    private static final NormalizationType NONE = NormalizationHandler.NONE;
    public static String usage = "simple-max [-r resolution] [-k norm] <file.hic> <loops.bedpe> <output.bedpe>\n" +
            "\t\tget the maximum pixel within the feature inside the bounds";
    private static final int buffer = 5;

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
        HiCZoom zoom = new HiCZoom(resolution);

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
                        MatrixZoomData zd = matrix.getZoomData(zoom);
                        if (zd != null) { // if it's not empty
                            try {
                                double[] nv = ds.getNormalizationVector(chrom1.getIndex(), zoom, VC).getData().getValues().get(0);
                                Set<Feature2D> newLoops = getMaximaForRegions(new HashSet<>(loops), resolution, buffer, zd, nv);

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

    public static Set<Feature2D> getMaximaForRegions(Set<Feature2D> loops, int resolution, int buffer, MatrixZoomData zd, double[] nv) {
        Set<Feature2D> newLoops = new HashSet<>();
        Collection<List<Feature2D>> loopGroups = QuickGrouping.groupNearbyRecords(
                loops, resolution * 200).values();
        for (List<Feature2D> group : loopGroups) {
            int minR = (int) ((FeatureStats.minStart1(group) / resolution) - buffer);
            int minC = (int) ((FeatureStats.minStart2(group) / resolution) - buffer);
            int maxR = (int) ((FeatureStats.maxEnd1(group) / resolution) + buffer + 1);
            int maxC = (int) ((FeatureStats.maxEnd2(group) / resolution) + buffer + 1);
            float[][] regionMatrix = Utils.getRegion(zd, minR, minC, maxR, maxC, NONE);
            for (Feature2D loop : group) {
                getTheMaxPixel(regionMatrix, loop, resolution, newLoops, minR, minC, nv, 2);
            }
            regionMatrix = null;
        }
        return newLoops;
    }

    public static void getTheMaxPixel(float[][] regionMatrix, Feature2D loop, int resolution,
                                      Set<Feature2D> newLoops, int minR, int minC, double[] nv,
                                      int window) {

        int r0 = (int) ((loop.getStart1() / resolution) - window - minR);
        int c0 = (int) ((loop.getStart2() / resolution) - window - minC);
        int rF = (int) ((loop.getEnd1() / resolution) + window + 1 - minR);
        int cF = (int) ((loop.getEnd2() / resolution) + window + 1 - minC);

        int[] binCoords;
        if (resolution > 10) { // or median 1
            binCoords = getMaxPixelInBox(regionMatrix, r0, rF, c0, cF, nv, minR, minC);
        } else {
            binCoords = null; // TODO getCentroidInBox(regionMatrix, r0, rF, c0, cF, nv, minR, minC);
        }

        if (binCoords != null) {
            long startX = (long) (binCoords[0] + minR) * resolution;
            long endX = startX + resolution;
            long startY = (long) (binCoords[1] + minC) * resolution;
            long endY = startY + resolution;
            Feature2D feature = new Feature2D(loop.getFeatureType(), loop.getChr1(), startX, endX,
                    loop.getChr2(), startY, endY, Color.BLACK, loop.getAttributes());
            newLoops.add(feature);
        }
    }

    public static int[] getMaxPixelInBox(float[][] matrix, int r0, int rF, int c0, int cF,
                                         double[] nv, int minR, int minC) {
        float maxCounts = 0;
        int[] coordinates = new int[]{-1, -1};
        for (int r = r0; r < rF; r++) {
            for (int c = c0; c < cF; c++) {
                if (matrix[r][c] > 0) {
                    float norm = (float) Math.sqrt(nv[minR + r] * nv[minC + c]);
                    if (norm > 0) {
                        float realVal = matrix[r][c] / norm;
                        if (realVal > maxCounts) {
                            maxCounts = realVal;
                            coordinates[0] = r;
                            coordinates[1] = c;
                        }
                    }
                }
            }
        }
        if (maxCounts > 0) {
            return coordinates;
        }
        return null;
    }
}
