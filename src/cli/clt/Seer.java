package cli.clt;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Seer {
    /* takes in one file currently (for ease of testing: can change later to a list of files and easily iterate over).
    resolution also input but currently set to 50 manually.
    */

    // run file with chromosome file, create main class, output as numpy (desktop)

    public static void calculateRowSums(String filename, Map<Chromosome, int[]> chromToRowSumsMap, int lowResolution, int highResolution) {

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, true);

        // todo hiC zoom objects, type & resolution fields
        //  (want the smallest number which is the resolution field). Iterate over.
        /*
        List<HiCZoom> objList = ds.getBpZooms();
        resolution = Integer.MAX_VALUE;
        for (HiCZoom obj : objList){
            if (obj. < resolution) {
                resolution = obj.
            }
        }
         */

        // iterate over a chromosome for now (chromosome 10)
        for (Chromosome chromosome : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(highResolution));
            if (zd == null) continue;

            int[] rowSummation = new int[(int) (chromosome.getLength() / highResolution + 1)];

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

    public static void run(String[] args) throws IOException {
        // check length of arguments equal to 3

        Map<Chromosome, int[]> chromToRowSumsMap = new HashMap<>();
        int lowResolution = 100;
        int highResolution = 50;
        calculateRowSums(args[1], chromToRowSumsMap, lowResolution, highResolution);
        UNIXTools.makeDir(args[2]);
        exportRowSumsToBedgraph(chromToRowSumsMap, args[2], highResolution);
        // MatrixTools.saveMatrixTextNumpy(outputFileName, results);
    }

    private static void exportRowSumsToBedgraph(Map<Chromosome, int[]> chromToRowSumsMap, String arg, int resolution) throws IOException {

        // todo first you need to make a bufferedfilewritere / filewriter

        // why is file giving an error? why do we need path in filename
        File outputFileName = new File(arg, "rowSums.bedgraph");
        outputFileName.createNewFile();
        // String outputFilePath = new File(arg, "rowSums.bedgraph").getAbsolutePath();

        FileWriter fw = new FileWriter(outputFileName);
        BufferedWriter bw = new BufferedWriter(fw);

        for (Chromosome chromosome : chromToRowSumsMap.keySet()) {
            int[] sums = chromToRowSumsMap.get(chromosome);

            for (int i = 0; i < sums.length; i++) {
                int startPosition = i * resolution;
                int endPosition = startPosition + resolution;
                int value = sums[i];
                if (value > 0) {
                    bw.write(chromosome.getName() + " " + startPosition + " " + endPosition + " " + value);
                    bw.newLine();
                }
            }

            // bin values --> sums[chromosome]
            // todo for every bin, you will write a line to the files
            // position = bin * resolution
            // chromosome " " startPosition + " " + endPosition + " " actual sum
        }

        // close the writer
        bw.close();

    }
}