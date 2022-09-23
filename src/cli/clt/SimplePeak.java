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
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class SimplePeak {

    private static final NormalizationType VC = NormalizationHandler.VC;
    private static final NormalizationType NONE = NormalizationHandler.NONE;
    public static String usage = "simple [-r resolution] [-k norm] <file.hic> <loops.bedpe> <output.bedpe>";

    public static void printUsageAndExit() {
        System.out.println("simple <file.hic> <loops.bedpe> <output.bedpe>");
        System.exit(7);
    }

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        String inputBedpeFile = args[2];
        String outputPath = args[3];
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpeFile, handler,
                false, null, false);

        Feature2DList localized = findSimpleMax(ds, loopList);
        localized.exportFeatureList(new File(outputPath), false, Feature2DList.ListFormat.NA);
        System.out.println("simple-max complete");
    }

    private static Feature2DList findSimpleMax(Dataset ds, Feature2DList loopList) {

        Feature2DList finalLoopList = new Feature2DList();
        final int buffer = 5;

        AtomicInteger currChromIndex = new AtomicInteger(0);
        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        ParallelizationTools.launchParallelizedCode(() -> {

            int threadIndex = currChromIndex.getAndIncrement();
            while (threadIndex < chromosomes.length) {

                Chromosome chrom1 = chromosomes[threadIndex];
                Set<Feature2D> retainedLoops = new HashSet<>(loopList.get(chrom1.getIndex(), chrom1.getIndex()));
                Set<Feature2D> finalLoops = new HashSet<>();

                if (retainedLoops.size() > 0) {
                    for (int resolution : new int[]{1000, 500, 200, 10}) {
                        HiCZoom zoom = new HiCZoom(resolution);
                        Matrix matrix = ds.getMatrix(chrom1, chrom1);
                        if (matrix != null) {
                            // load matrix in "resolution" resolution of that chromosome in sparse format
                            MatrixZoomData zd = matrix.getZoomData(zoom);
                            if (zd != null) { // if it's not empty
                                try {
                                    double[] nv = ds.getNormalizationVector(chrom1.getIndex(), zoom, VC).getData().getValues().get(0);
                                    Set<Feature2D> currentLoops = new HashSet<>();
                                    Collection<LinkedList<Feature2D>> loopGroups = FusionTools.groupNearbyRecords(
                                            retainedLoops, resolution * 200).values();
                                    for (LinkedList<Feature2D> group : loopGroups) {
                                        int minR = (int) ((FeatureStats.minStart1(group) / resolution) - buffer);
                                        int minC = (int) ((FeatureStats.minStart2(group) / resolution) - buffer);
                                        int maxR = (int) ((FeatureStats.maxEnd1(group) / resolution) + buffer + 1);
                                        int maxC = (int) ((FeatureStats.maxEnd2(group) / resolution) + buffer + 1);
                                        float[][] regionMatrix = Utils.getRegion(zd, minR, minC, maxR, maxC, NONE);
                                        for (Feature2D loop : group) {
                                            SimpleMax.getTheMaxPixel(regionMatrix, loop, resolution, currentLoops, minR, minC, nv, 2);
                                        }
                                        regionMatrix = null;
                                    }

                                    retainedLoops.clear();
                                    if (resolution == 10) {
                                        finalLoops.addAll(currentLoops);
                                    } else {
                                        retainedLoops.addAll(currentLoops);
                                    }
                                    currentLoops.clear();


                                } catch (Exception e) {
                                    e.printStackTrace();
                                    System.exit(76);
                                }
                                matrix.clearCacheForZoom(zoom);
                            }
                            matrix.clearCache();

                            synchronized (finalLoopList) {
                                finalLoopList.addByKey(Feature2DList.getKey(chrom1, chrom1), new ArrayList<>(finalLoops));
                            }
                        }
                    }
                }

                threadIndex = currChromIndex.getAndIncrement();
            }
        });

        return finalLoopList;
    }
}
