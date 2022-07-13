package cli.clt;

import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Seer {
    /* takes in one file currently (for ease of testing: can change later to a list of files and easily iterate over).
    resolution also input but currently set to 50 manually.
    */

    // run file with chromosome file, create main class, output as numpy (desktop)
    public static void calculateRowSums(String filename, Map<Chromosome, int[]> chromToRowSumsMap) {

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, true);

        int resolution = 50;

        // todo hiC zoom objects, type & resolution fields
        //  (want the smallest number which is the resolution field). Iterate over.
        ds.getBpZooms();

        // iterate over a chromosome for now (chromosome 10)
        for (Chromosome chromosome : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) continue;

            int[] rowSummation = new int[(int) (chromosome.getLength() / resolution + 1)];

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

            chromToRowSumsMap.put(chromosome, rowSummation);
        }
    }
    private static void printUsageAndExit() {
        System.out.println("seer <file> <out_folder>");
        System.exit(19);
    }

    public static void run(String[] args) {
        // check length of arguments equal to 3

        Map<Chromosome, int[]> chromToRowSumsMap = new HashMap<>();
        calculateRowSums(args[1], chromToRowSumsMap);
        UNIXTools.makeDir(args[2]);
        exportRowSumsToBedgraph(chromToRowSumsMap, args[2]);
        // MatrixTools.saveMatrixTextNumpy(outputFileName, results);
    }

    private static void exportRowSumsToBedgraph(Map<Chromosome, int[]> chromToRowSumsMap, String arg) {

        // String outputFileName = new File(args[2], "rowSums.bedgraph").getAbsolutePath();
        // todo first you need to make a bufferedfilewritere / filewriter

        for (Chromosome chromosome : chromToRowSumsMap.keySet()) {
            int[] sums = chromToRowSumsMap.get(chromosome);

            // todo for every bin, you will write a line to the files
            // position = bin * resolution
            // chromosome " " startPosition + " " + endPosition + " " actual sum


        }

        // close the writer

    }
}
