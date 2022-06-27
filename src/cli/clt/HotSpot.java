package cli.clt;

import cli.utils.StandardDevUtils;
import javastraw.reader.Dataset;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.util.Iterator;
import java.io.IOException;
import java.util.List;

public class HotSpot {
    private static List<float[][]> getMatrices(int res, int window, String strNorm, String filename, String[] names) {
        boolean useCache = false;
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
        int chrLength = (int) chr5.getLength();

        int binXStart = 1 / resolution;
        int binXEnd = 1 / resolution;
        int binYStart = chrLength / resolution;
        int binYEnd = chrLength / resolution;
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

    /*
    private static List<float[][]> getMatrices(int res, int window, String strNorm, String[] files, String[] names) {
        boolean useCache = false;

        for (String filename : files) {

            // create a hic dataset object
            Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, useCache, true);
            // choose norm: we know the datasets we're using will have SCALE available
            NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{strNorm, "NONE"});

            // choose to use 2Kb resolution
            int resolution = res;

            // Get a list of the chromosome objects
            // Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();

            for (Chromosome chromosome : chromosomes) {

                Matrix matrix = ds.getMatrix(chromosome, chromosome);
                if (matrix == null) continue;
                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                if (zd == null) continue;

                boolean getDataUnderTheDiagonal = false;
                // zd is now a data structure that contains pointers to the data
                // *** Let's show 2 different ways to access data ***

                // OPTION 1
                // iterate on all the data for the whole chromosome in sparse format
                Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
                while (iterator.hasNext()) {
                    ContactRecord record = iterator.next();
                    // now do whatever you want with the contact record
                    int binX = record.getBinX();
                    int binY = record.getBinY();
                    float counts = record.getCounts();

                    // binX and binY are in BIN coordinates, not genome coordinates
                    // to switch, we can just multiply by the resolution
                    int genomeX = binX * resolution;
                    int genomeY = binY * resolution;

                    if (counts > 0) { // will skip NaNs

                        // do task

                        // the iterator only iterates above the diagonal
                        // to also fill in data below the diagonal, flip it
                        if (binX != binY) {
                            binX = record.getBinY();
                            binY = record.getBinX();
                            counts = record.getCounts();

                        }
                    }
                }
            }

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
    */

    private static void printUsageAndExit() {
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
        print(stdDevObj.StdDeviationFinder(mtxList));
    }
}}