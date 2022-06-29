package cli.clt;

import cli.utils.StandardDevUtils;
import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.util.List;

import static cli.utils.StandardDevUtils.StdDeviationFinder;

public class HotSpot {

    private static List<float[][]> getMatrices(int res, int window, String strNorm, String[] files) {
        /*
        Accepts parameters of the HOTSPOT command line tool and returns a list of 2D float arrays that represent the Hi-C maps corresponding to the files. The 2D float
        arrays will have pixel sizes equal window argument
         */
        boolean useCache = false;
        List<float[][]> mtxList = null;
        for (String file : files) {
            // create a hic dataset object
            Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, useCache, true);
            // choose norm: we know the datasets we're using will have SCALE available
            NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{strNorm, "NONE"});

            // choose to use 2Kb resolution
            int resolution = res;

            // Get a list of the chromosome objects
            Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();

            // Instantiates the 2D float array that will be used to represent the Hi-C map of individual tissue genomes
            int genomeLength = 0;
            for (Chromosome chromosome : chromosomes) {
                genomeLength += chromosome.getLength();
            }
            int genomeBinsLength = (int) Math.ceil(genomeLength / (double) resolution);
            float[][] singleMatrix = new float[genomeBinsLength][genomeBinsLength];

            boolean getDataUnderTheDiagonal = true;

            // Iterates through each possible chromosome pair
            // ASSUMPTION: (for prevAbsYBinCoord) that chromosomes is a list in order from chromosome 0 to 13
            int prevAbsYBinCoord = 0;
            for (int i = 0; i < chromosomes.length; i++) {
                int prevAbsXBinCoord = 0;
                for (int j = i; i < chromosomes.length; i++) {
                    Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[j]);
                    if (matrix == null) continue;
                    MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                    if (zd == null) continue;

                    long binXStart = 0, binYStart = 0;
                    // TODO: deal with edge cases (length doesn't divide evenly by resolution)
                    long binXEnd = chromosomes[i].getLength() / resolution;
                    long binYEnd = chromosomes[j].getLength() / resolution;

                    // Iterates through the blocks of the chromosome pair, eventually grabbing the bin coordinates and count to record them in singleMatrix
                    List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, norm, getDataUnderTheDiagonal);
                    for (Block b : blocks) {
                        if (b != null) {
                            for (ContactRecord rec : b.getContactRecords()) {
                                if (rec.getCounts() > 0) { // will skip NaNs
                                    int binX = rec.getBinX();
                                    int binY = rec.getBinY();
                                    singleMatrix[binX + prevAbsXBinCoord][binY + prevAbsYBinCoord] = rec.getCounts();
                                }
                            }
                        }

                    }
                }
                // adds the singleMatrix to the list of overall matrices, mtxList
                mtxList.add(singleMatrix);
            }
        }
        return mtxList;
    }

    private static void printUsageAndExit() {
        /* example print: ("apa [--min-dist minval] [--max-dist max_val] [--window window] [-r resolution]" +
                " [-k NONE/VC/VC_SQRT/KR] [--corner-width corner_width] [--include-inter include_inter_chr] [--ag-norm]" +
                " <input.hic> <loops.bedpe> <outfolder>"); */
        System.out.println("hotspot [--res resolution] [--window window] [--norm normalization] <file1.hic,file2.hic,...> <name1,name2,...> <out_folder>");
        System.exit(19);
    }

    public static void run(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            printUsageAndExit();
        }
        final int DEFAULT_RES = 2000;
        final int DEFAULT_WINDOW = 0;
        // hotspot [--res int] [--window int] [--norm string] <file1.hic,file2.hic,...> <name1,name2,...> <out_folder>
        String[] files = args[2].split(",");
        List<float[][]> mtxList = getMatrices(parser.getResolutionOption(DEFAULT_RES), parser.getWindowSizeOption(DEFAULT_WINDOW), parser.getNormalizationStringOption(), files);
        StandardDevUtils stdDevObj = new StandardDevUtils();
        float[][] stdDevArray = new float[mtxList.get(0).length][mtxList.get(0).length];
        stdDevArray = StdDeviationFinder(mtxList);

        // saves 2D float array as a NumpyMatrix to outfolder
        // ASSUMPTION: in command line, outfolder argument is inputted as an entire path
        MatrixTools.saveMatrixTextNumpy(args[3], stdDevArray);
    }
}