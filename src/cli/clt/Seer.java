package cli.clt;

import javastraw.reader.Dataset;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

import java.util.Iterator;

public class Seer {
    /* takes in one file currently (for ease of testing: can change later to a list of files and easily iterate over).
    resolution also input but currently set to 50 manually.
    */

    // run file with chromosome file, create main class, output as numpy (desktop)
    public static int[] rowSum(String filename) {

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, true);

        int resolution = 50;

        // hiC zoom objects, type & resolution fields (want the smallest number which is the resolution field). Iterate over.
        ds.getBpZooms();
        // first find the genome length
        /*
        int genomeLength = 0;

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; i < chromosomes.length; i++) {
                genomeLength++;
            }
        }
        */

        // iterate over a chromosome for now (chromosome 10)

        Chromosome chr10 = ds.getChromosomeHandler().getChromosomeFromName("chr10");
        Matrix matrix = ds.getMatrix(chr10, chr10);
        if (matrix == null) return null;
        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
        if (zd == null) return null;

        int[] rowSummation = new int[(int) (chr10.getLength() / resolution + 1)];

        Iterator<ContactRecord> iterator = zd.getDirectIterator();
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            float counts = record.getCounts();

            if (counts > 0) { // will skip NaNs
                int binX = record.getBinX();
                int binY = record.getBinY();
                rowSummation[binX] += record.getCounts();

                if (binX != binY) {
                    rowSummation[binY] += record.getCounts();
                    // do task
                }
            }
        }

        // returns the row summation matrix
        return rowSummation;
    }
    private static void printUsageAndExit() {
        System.out.println("seer <file> <out_folder>");
        System.exit(19);
    }

    public static void run(String[] args) {
        // check length of arguments equal to 3

        int[] results = rowSum(args[1]);
        MatrixTools.saveMatrixTextNumpy(args[2], results);
    }
}
