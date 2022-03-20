package cli.apa;

import javastraw.reader.basics.Chromosome;

public class RegionConfiguration {
    private final Chromosome chrom1, chrom2;

    public RegionConfiguration(Chromosome chrom1, Chromosome chrom2) {
        this.chrom1 = chrom1;
        this.chrom2 = chrom2;
    }

    public Chromosome getChr1() {
        return chrom1;
    }

    public Chromosome getChr2() {
        return chrom2;
    }
}
