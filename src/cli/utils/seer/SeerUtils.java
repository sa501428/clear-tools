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

public class SeerUtils {

    public static void exportRowSumsToBedgraph(Map<Chromosome, int[]> chromToRowSumsMap, String arg, int resolution) throws IOException {
        File outputFileName = new File(arg, "rowSums.bedgraph");
        outputFileName.createNewFile();

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
        }

        // close the writer
        bw.close();

    }

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

    public static SimpleLocation updateToHigherResPosition(SimpleLocation genomePosition, int[] hiResRowSums,
                                                           int lowResolution, int highResolution) {

        return null;
    }
}
