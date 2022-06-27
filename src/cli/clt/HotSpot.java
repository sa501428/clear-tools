package cli.clt;

import cli.utils.StandardDevUtils;
import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.io.IOException;
import java.util.List;

public class HotSpot {
    /*
    private static void tissueTypeToNumpyFile(String tissue) {

        boolean useCache = false;

        //QUESTION:
        // https://s3.us-east-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter1/TISSUE_nd.hic
        String filename = "https://s3.us-east-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter1/" + tissue + "_nd.hic"; //insert filename here. I don't think URL is acceptable

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, useCache, true);
        // choose norm: we know the datasets we're using will have SCALE available
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "NONE"});

        // choose to use 2Kb resolution
        int resolution = 2000;

        // Get a list of the chromosome objects
        Chromosome chr5 = ds.getChromosomeHandler().getChromosomeFromName("chr5");

        Matrix matrix = ds.getMatrix(chr5, chr5);
        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));

        boolean getDataUnderTheDiagonal = true;

        // our bounds will be binXStart, binYStart, binXEnd, binYEnd
        // these are in BIN coordinates, not genome coordinates
        int binXStart = 119848237 / resolution;
        int binYStart = 119836267 / resolution;
        int binXEnd = 121000236 / resolution;
        int binYEnd = 120988666 / resolution;

        int numRows = binXEnd - binXStart + 1; // replace with actual number later
        int numCols = binYEnd - binYStart + 1; // replace later

        float[][] float2DArray = new float[0][];
        try {
            float2DArray = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, binXStart, binXEnd, binYStart, binYEnd, numRows, numCols, norm, getDataUnderTheDiagonal);
            MatrixTools.saveMatrixTextNumpy(tissue + "numpymatrix.to.output.npy", float2DArray);
        } catch (IOException e) {
            //throw new RuntimeException(e);
            System.err.println(e.getLocalizedMessage());
            //System.exit(10);
        }
    }
    */
    private static List<float[][]> getMatrices(int res, int window, String strNorm, String[] files, String[] names) {
        /*
        Accepts parameters of the HOTSPOT command line tool and returns a list of 2D float arrays that represent the Hi-C maps corresponding to the files. The 2D float
        arrays will have pixel sizes equal window argument
         */
        boolean useCache = false;

        for (String filename : files) {

            // create a hic dataset object
            Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, useCache, true);
            // choose norm: we know the datasets we're using will have SCALE available
            NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{strNorm, "NONE"});

            // choose to use 2Kb resolution
            int resolution = res;

            // Get a list of the chromosome objects
            Chromosome chr5 = ds.getChromosomeHandler().getChromosomeFromName("chr5");

            Matrix matrix = ds.getMatrix(chr5, chr5);
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));

            boolean getDataUnderTheDiagonal = true;

            // our bounds will be binXStart, binYStart, binXEnd, binYEnd
            // these are in BIN coordinates, not genome coordinates
            int binXStart = 119848237 / resolution;
            int binYStart = 119836267 / resolution;
            int binXEnd = 121000236 / resolution;
            int binYEnd = 120988666 / resolution;

            int numRows = binXEnd - binXStart + 1; // replace with actual number later
            int numCols = binYEnd - binYStart + 1; // replace later

            float[][] float2DArray = new float[0][];
            try {
                float2DArray = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, binXStart, binXEnd, binYStart, binYEnd, numRows, numCols, norm, getDataUnderTheDiagonal);
                MatrixTools.saveMatrixTextNumpy(tissue + "numpymatrix.to.output.npy", float2DArray);
            } catch (IOException e) {
                //throw new RuntimeException(e);
                System.err.println(e.getLocalizedMessage());
                //System.exit(10);
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
        if (args.length != 6 || args[3].length() != args[4].length()) {
            printUsageAndExit();
        }
        final int DEFAULT_RES = 0;
        final int DEFAULT_WINDOW = 0;
        // hotspot [--res int] [--window int] [--norm string] <file1.hic,file2.hic,...> <name1,name2,...> <out_folder>
        List<float[][]> mtxList = getMatrices(parser.getResolutionOption(DEFAULT_RES), parser.getWindowSizeOption(DEFAULT_WINDOW), parser.getNormalizationStringOption(), args[3], args[4]);
        StandardDevUtils stdDevObj = new StandardDevUtils();
        // return float[][] stdDevObj.StdDeviationFinder(mtxList);
    }
}