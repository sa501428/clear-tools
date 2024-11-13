package cli.utils.general;

import javastraw.reader.basics.Chromosome;

import java.util.HashSet;
import java.util.Set;

public class TranslocationSet {
    private final Set<InterChromosomeRegion> regions = new HashSet<>();

    public void add(Chromosome chr1, Chromosome chr2) {
        regions.add(new InterChromosomeRegion(chr1, chr2));
    }

    public boolean contains(Chromosome chr1, Chromosome chr2) {
        return regions.contains(new InterChromosomeRegion(chr1, chr2));
    }
}
