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
import javastraw.reader.norm.NormalizationPicker;
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
                                // newLoops is the output looplist that has the calibrated centers
                                List<Feature2D> newLoops = new ArrayList<>();
                                for (Feature2D loop : loops) {
                                    // loop through the looplist
                                    getTheMaxPixel(loop, zd, resolution, newLoops, norm);
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

    private static void getTheMaxPixel(Feature2D loop, MatrixZoomData zd, int resolution, List<Feature2D> newLoops,
                                       NormalizationType norm) {
        // Initialize variables
        // can assume loop from top right corner of diagonal
        long binXStart = loop.getStart1() / resolution;
        long binYStart = loop.getStart2() / resolution;
        long binXEnd = (loop.getEnd1() / resolution) + 1;
        long binYEnd = (loop.getEnd2() / resolution) + 1;

        boolean getDataUnderDiagonal = true;

        List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, norm, getDataUnderDiagonal);
        int[] binCoords = getMaxPixelInBox(blocks, binXStart, binYStart, binXEnd, binYEnd);
        blocks.clear();

        // arbitrarily chose binCoordsNONE, since NONE and VC have same coordinates
        long startX = (long) binCoords[0] * resolution;
        long endX = startX + resolution;
        long startY = (long) binCoords[1] * resolution;
        long endY = startY + resolution;
        Feature2D feature = new Feature2D(loop.getFeatureType(), loop.getChr1(), startX, endX,
                loop.getChr2(), startY, endY, Color.BLACK, loop.getAttributes());
        newLoops.add(feature);

    }

    private static int[] getMaxPixelInBox(List<Block> blocks, long binXStart, long binYStart, long binXEnd, long binYEnd) {
        float maxCounts = 0;
        int[] coordinates = new int[2];
        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    // try checking if bounds fixed bug
                    if (inBounds(rec, binXStart, binYStart, binXEnd, binYEnd)) {
                        if (rec.getCounts() > maxCounts) { // will skip NaNs
                            maxCounts = rec.getCounts();
                            coordinates[0] = rec.getBinX();
                            coordinates[1] = rec.getBinY();
                        }
                    }
                }
            }
        }
        return coordinates;
    }

    private static boolean inBounds(ContactRecord rec, long binXStart, long binYStart, long binXEnd, long binYEnd) {
        int x = rec.getBinX();
        int y = rec.getBinY();
        return x >= binXStart && x < binXEnd && y >= binYStart && y < binYEnd;
    }
}
