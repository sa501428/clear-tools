package cli.utils.anchors;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GlobalAnchorCounter {

    public static Map<Chromosome, int[]> createFromPeaks(Chromosome[] chromosomes, Map<Chromosome,
            Map<Chromosome, List<ContactRecord>>> globalPeakMap, int resolution) {
        Map<Chromosome, int[]> countsForAnchor = new HashMap<>();
        for (Chromosome chromosome : chromosomes) {
            int[] counts = new int[(int) (chromosome.getLength() / resolution) + 1];
            countsForAnchor.put(chromosome, counts);
        }

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i + 1; j < chromosomes.length; j++) {
                populate(countsForAnchor.get(chromosomes[i]),
                        countsForAnchor.get(chromosomes[j]),
                        globalPeakMap.get(chromosomes[i]).get(chromosomes[j]));
            }
        }
        return countsForAnchor;
    }

    private static void populate(int[] counts1, int[] counts2,
                                 List<ContactRecord> records) {
        for (ContactRecord record : records) {
            counts1[record.getBinX()]++;
            counts2[record.getBinY()]++;
        }
    }
}
