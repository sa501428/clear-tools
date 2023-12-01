package cli.utils;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;

public class StrawUtils {
    public static Chromosome[] getChromosomes(ChromosomeHandler handler, String chrom) {
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        if (chrom != null && chrom.length() > 0) {
            String[] chroms = chrom.split(",");
            chromosomes = new Chromosome[chroms.length];
            for (int i = 0; i < chroms.length; i++) {
                chromosomes[i] = handler.getChromosomeFromName(chroms[i]);
            }
        }
        return chromosomes;
    }
}
