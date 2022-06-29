package cli.clt;

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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static cli.utils.StandardDevUtils.StdDeviationFinder;

public class HotSpot {

    private static void getMatrices(int resolution, int window, String strNorm, String file,
                                    Map<Integer, float[][]> results) {
        /*
        Accepts parameters of the HOTSPOT command line tool and returns a list of 2D float arrays that represent the Hi-C maps corresponding to the files. The 2D float
        arrays will have pixel sizes equal window argument
         */
        boolean useCache = false;
        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(file, false, useCache, true);
        // choose norm: we know the datasets we're using will have SCALE available
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{strNorm, "NONE"});

        // Get a list of the chromosome objects
        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();

        // Instantiates the 2D float array that will be used to represent the Hi-C map of individual tissue genomes

        boolean getDataUnderTheDiagonal = true;

        // Iterates through each possible chromosome pair
        // ASSUMPTION: (for prevAbsYBinCoord) that chromosomes is a list in order from chromosome 0 to 13
        //int prevAbsYBinCoord = 0;


        for (int i = 0; i < chromosomes.length; i++) {
            //int prevAbsXBinCoord = 0;

            Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[i]);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) continue;

            int binXStart = 3000000 / resolution; // todo make this slide across the map
            // TODO: deal with edge cases (length doesn't divide evenly by resolution)
            long binXEnd = binXStart + window;

            float[][] singleMatrix = results.get(chromosomes[i].getIndex());

            // Iterates through the blocks of the chromosome pair, eventually grabbing the bin coordinates and count to record them in singleMatrix
            List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binXStart, binXEnd, binXEnd, norm, getDataUnderTheDiagonal);
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) { // will skip NaNs
                            int binX = rec.getBinX() - binXStart;
                            int binY = rec.getBinY() - binXStart;
                            singleMatrix[binX][binY] = rec.getCounts();
                        }
                    }
                }
            }
        }
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
        final int DEFAULT_WINDOW = 1000;
        // hotspot [--res int] [--window int] [--norm string] <file1.hic,file2.hic,...> <name1,name2,...> <out_folder>
        String[] files = args[2].split(",");

        Map<Integer, float[][]> results = new HashMap<>();

        int window = parser.getWindowSizeOption(DEFAULT_WINDOW);

        Dataset ds = HiCFileTools.extractDatasetForCLT(files[0], false, false, true);
        for (Chromosome chromosome : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            results.put(chromosome.getIndex(), new float[window][window]);
        }

        for (String file : files) {
            HotSpot.getMatrices(parser.getResolutionOption(DEFAULT_RES),
                    window, parser.getNormalizationStringOption(), file,
                    results);
        }
        //StandardDevUtils stdDevObj = new StandardDevUtils();
        //float[][] stdDevArray = new float[mtxList.get(0).length][mtxList.get(0).length];

        float[][] stdDevArray = StdDeviationFinder(mtxList);

        // saves 2D float array as a NumpyMatrix to outfolder
        // ASSUMPTION: in command line, outfolder argument is inputted as an entire path
        MatrixTools.saveMatrixTextNumpy(args[3], stdDevArray);
    }
}