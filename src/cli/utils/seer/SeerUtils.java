package cli.utils.seer;

import cli.utils.sift.SimpleLocation;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

public class SeerUtils {

    // create bedgraph file from a rowSums map
    public static void exportRowSumsToBedgraph(Map<Chromosome, int[]> chromToRowSumsMap,
                                               String path, int resolution) throws IOException {
        File outputFileName = new File(path);
        outputFileName.createNewFile();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
        // write for every chromosome
        for (Chromosome chromosome : chromToRowSumsMap.keySet()) {
            int[] sums = chromToRowSumsMap.get(chromosome);
            // iterate through bins, write in the format <chr> <start_pos> <end_pos> <value>
            for (int i = 0; i < sums.length; i++) {
                int startPosition = i * resolution;
                int endPosition = startPosition + resolution;
                if (sums[i] > 0) {
                    bw.write(chromosome.getName() + " " + startPosition + " " + endPosition + " " + sums[i]);
                    bw.newLine();
                }
            }
        }
        // close the writer
        bw.close();
    }

    public static void exportRowFloatsToBedgraph(Map<Chromosome, float[]> chromToRowSumsMap,
                                                 String path, int resolution) throws IOException {
        File outputFileName = new File(path);
        outputFileName.createNewFile();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
        // write for every chromosome
        for (Chromosome chromosome : chromToRowSumsMap.keySet()) {
            float[] sums = chromToRowSumsMap.get(chromosome);
            // iterate through bins, write in the format <chr> <start_pos> <end_pos> <value>
            for (int i = 0; i < sums.length; i++) {
                int startPosition = i * resolution;
                int endPosition = startPosition + resolution;
                if (sums[i] > 0) {
                    bw.write(chromosome.getName() + " " + startPosition + " " + endPosition + " " + sums[i]);
                    bw.newLine();
                }
            }
        }
        // close the writer
        bw.close();
    }

    public static void exportRowIntsToBedgraph(Map<Chromosome, int[]> chromToRowSumsMap, String path, int resolution) throws IOException {
        File outputFileName = new File(path);
        outputFileName.createNewFile();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
        // write for every chromosome
        for (Chromosome chromosome : chromToRowSumsMap.keySet()) {
            int[] sums = chromToRowSumsMap.get(chromosome);
            // iterate through bins, write in the format <chr> <start_pos> <end_pos> <value>
            for (int i = 0; i < sums.length; i++) {
                int startPosition = i * resolution;
                int endPosition = startPosition + resolution;
                if (sums[i] > 0) {
                    bw.write(chromosome.getName() + " " + startPosition + " " + endPosition + " " + sums[i]);
                    bw.newLine();
                }
            }
        }
        // close the writer
        bw.close();
    }

    // calculates rowSums for a chromosome
    public static int[] getRowSumsForZD(Chromosome chromosome, int highResolution, Iterator<ContactRecord> iterator) {
        int[] rowSums = new int[(int) (chromosome.getLength() / highResolution + 1)];
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            float counts = record.getCounts();
            if (counts > 0) { // will skip NaNs
                int binX = record.getBinX();
                int binY = record.getBinY();
                rowSums[binX] += record.getCounts();
                if (binX != binY) {
                    rowSums[binY] += record.getCounts();
                }
            }
        }
        return rowSums;
    }

    public static SimpleLocation updateToHigherResPosition(SimpleLocation genomePosition, double[] hiResCDF,
                                                           int lowResolution, int highResolution, Random rand) {
        int window = lowResolution / highResolution; // e.g. 100
        int startBinX = genomePosition.getBinX() / highResolution;
        int startBinY = genomePosition.getBinY() / highResolution;

        int genomeX = getHigherQualityIndex(startBinX, window, hiResCDF, rand) * highResolution;
        int genomeY = getHigherQualityIndex(startBinY, window, hiResCDF, rand) * highResolution;

        return new SimpleLocation(genomeX, genomeY);
    }

    private static int getHigherQualityIndex(int startBin, int window, double[] hiResCDF, Random rand) {
        double r = rand.nextDouble();

        int realStartBin = startBin;
        int realEndBin = startBin + window;

        if (realEndBin >= hiResCDF.length) {
            realEndBin = hiResCDF.length - 1;
        }

        if (startBin > 0) {
            realStartBin = startBin - 1;
        }

        double x0 = hiResCDF[realStartBin];
        double xF = hiResCDF[realEndBin];

        double range = xF - x0;
        double target = r * range + x0;

        return BinarySearch.runBinarySearchIteratively(hiResCDF, target, realStartBin, realEndBin);
    }

    public static double[] convertToCDF(int[] numbers) {
        double[] cumSums = new double[numbers.length];
        // matrix of cumulative sums of the first array
        cumSums[0] = numbers[0];
        for (int i = 1; i < numbers.length; i++) {
            cumSums[i] = cumSums[i - 1] + numbers[i];
        }
        double sum = cumSums[cumSums.length - 1];
        // normalizes every term by having last value = 1
        for (int i = 0; i < cumSums.length; i++) {
            cumSums[i] /= sum;
        }
        return cumSums;
    }
}
