package cli.clt;

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
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class SimpleMax {

    private static final NormalizationType VC = NormalizationHandler.VC;

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

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpeFile, handler,
                false, null, false);


        Feature2DList localized = findSimpleMax(ds, resolution, loopList);
        localized.exportFeatureList(new File(outputPath), false, Feature2DList.ListFormat.NA);
        System.out.println("simple-max complete");
    }

    private static Feature2DList findSimpleMax(Dataset ds, int resolution, Feature2DList loopList) {
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
                        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                        if (zd != null) {
                            try {
                                List<Feature2D> newLoops = new ArrayList<>();
                                for (Feature2D loop : loops) {
                                    Feature2D newLoop = getTheMaxPixel(loop, zd, resolution);
                                    newLoops.add(newLoop);
                                }
                                synchronized (finalLoopList) {
                                    finalLoopList.addByKey(Feature2DList.getKey(chrom1, chrom1), newLoops);
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

    private static Feature2D getTheMaxPixel(Feature2D loop, MatrixZoomData zd, int resolution) {
        // todo Justin
        // type 2 access of the hic data
        // loop has genome coordinates
        // grab the matrix corresponding the loop
        // find the max pixel when there is no NONE normalization
        // find the max pixel when there is VC normalization
        // ideally these should be the same
        // if they are not, print them out and let's debug
        // if they are, make a new loop based on this new location

        /*
        float[][] obsMatrix = new float[matrixWidth][matrixWidth];
        Utils.addLocalizedData(obsMatrix, zd, loop, matrixWidth, resolution, window, norm, key);


         */
        return null;
    }

    private List<Feature2D> locationsOfMaxVals(int resolution, Chromosome chromosome,
                                               Dataset ds, List<Feature2D> feature2DList,
                                               NormalizationType norm) {
        List<Feature2D> newList = new ArrayList<>();
        Matrix matrix = ds.getMatrix(chromosome, chromosome);
        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, resolution));

        for (Feature2D loop : feature2DList) {
            long start1 = loop.getStart1() / resolution - 1;
            long start2 = loop.getStart2() / resolution - 1;
            long end1 = loop.getEnd1() / resolution + 1;
            long end2 = loop.getEnd2() / resolution + 1;

            try {
                double[][] data = HiCFileTools.extractLocalBoundedRegion(zd, start1, end1,
                        start2, end2, (int) (end1 - start1), (int) (end2 - start2),
                        norm, false).getData();
                System.out.println(chromosome.getName());
                System.out.println(start1 * resolution + " " + end1 * resolution);
                System.out.println(start2 * resolution + " " + end2 * resolution);
                System.out.println();

                System.exit(6);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return newList;
    }

}
