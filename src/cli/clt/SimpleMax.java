package cli.clt;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class SimpleMax {

    private static final NormalizationType VC = NormalizationHandler.VC;
    private static final NormalizationType NONE = NormalizationHandler.NONE;
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
                        // load matrix in "resolution" resolution of that chromosome in sparse format
                        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                        if (zd != null) { // if it's not empty
                            try {
                                // newLoops is the output looplist that has the calibrated centers
                                List<Feature2D> newLoops = new ArrayList<>();
                                for (Feature2D loop : loops) {
                                    // loop through the looplist
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

        // Initialize variables
        long binXStart = loop.getStart2() / resolution;
        long binYStart = loop.getStart1() / resolution;
        long binXEnd = loop.getEnd2() / resolution;
        long binYEnd = loop.getEnd1() / resolution;

        /*
        PREVIOUSLY:
        long binXStart = loop.getStart1() / resolution;
        long binYStart = loop.getEnd1() / resolution;
        long binXEnd = loop.getStart2() / resolution;
        long binYEnd = loop.getEnd2() / resolution;
         */
        boolean getDataUnderDiagonal = true;

        // VC normalization
        List<Block> blocksVC = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, VC, getDataUnderDiagonal);
        float counts = 0;
        int[] binCoordsVC = new int[2];
        for (Block b : blocksVC) {
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    if (rec.getCounts() > 0) { // will skip NaNs
                        // can choose to use the BIN coordinates
                        // int binX = rec.getBinX();
                        // int binY = rec.getBinY();
                        if (rec.getCounts() > counts) {
                            counts = rec.getCounts();
                            binCoordsVC[0] = rec.getBinX();
                            binCoordsVC[1] = rec.getBinY();
                        }
                    }
                }
            }
        }

        // NONE normalization
        List<Block> blocksNONE = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, NONE, getDataUnderDiagonal);
        // reset counts back to 0
        // we can reuse variable since we only care about bin #
        counts = 0;
        int[] binCoordsNONE = new int[2];
        for (Block b : blocksNONE) {
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    if (rec.getCounts() > 0) { // will skip NaNs
                        if (rec.getCounts() > counts) {
                            counts = rec.getCounts();
                            binCoordsNONE[0] = rec.getBinX();
                            binCoordsNONE[1] = rec.getBinY();
                        }
                    }
                }
            }
        }

        if (binCoordsNONE[0] != binCoordsVC[0] || binCoordsNONE[1] != binCoordsVC[1]) {
            System.out.println("ERROR, NONE AND VC REPORT DIFFERENT MAX PIXELS");
            // debugging below time
            long startX_VC = (long) binCoordsVC[0] * resolution;
            long endX_VC = startX_VC + resolution;
            long startY_VC = (long) binCoordsVC[1] * resolution;
            long endY_VC = startY_VC + resolution;
            // this is just for debugging:
            long startX_NONE = (long) binCoordsNONE[0] * resolution;
            long endX_NONE = startX_NONE + resolution;
            long startY_NONE = (long) binCoordsNONE[1] * resolution;
            long endY_NONE = startY_NONE + resolution;
            Feature2D feature = new Feature2D(loop.getFeatureType(), loop.getChr1(), startX_VC, endX_VC, loop.getChr2(), startY_VC, endY_VC, Color.BLACK, loop.getAttributes());

            /*
            System.out.println("VC Normalization: " +
                    String.valueOf(startX_VC) + " " +
                    String.valueOf(endX_VC)) + " " +
                    String.valueOf(startY_VC) + " " +
                    String.valueOf(endY_VC));


             */
            return feature;
        } else {
            // arbitrarily chose binCoordsNONE, since NONE and VC have same coordinates
            long startX = (long) binCoordsNONE[0] * resolution;
            long endX = startX + resolution;
            long startY = (long) binCoordsNONE[1] * resolution;
            long endY = startY + resolution;
            Feature2D feature = new Feature2D(loop.getFeatureType(), loop.getChr1(), startX, endX, loop.getChr2(), startY, endY, Color.BLACK, loop.getAttributes());
            return feature;
        }
    }

        /*
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


         */
}
