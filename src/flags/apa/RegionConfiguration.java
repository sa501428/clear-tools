package flags.apa;

import javastraw.reader.basics.Chromosome;

public class RegionConfiguration {
    private final Chromosome chrom1, chrom2;
    private final int distIndex;

    public RegionConfiguration(Chromosome chrom1, Chromosome chrom2, int distIndex) {
        this.chrom1 = chrom1;
        this.chrom2 = chrom2;
        this.distIndex = distIndex;
    }

    public Chromosome getChr1() {
        return chrom1;
    }

    public Chromosome getChr2() {
        return chrom2;
    }

    public int getDistIndex() {
        return distIndex;
    }
}
