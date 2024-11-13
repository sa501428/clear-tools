package cli.utils.general;

import javastraw.reader.basics.Chromosome;

import java.util.Objects;

public class InterChromosomeRegion {
    Chromosome c1, c2;

    public InterChromosomeRegion(Chromosome c1, Chromosome c2) {
        if (c1.getIndex() <= c2.getIndex()) {
            this.c1 = c1;
            this.c2 = c2;
        } else {
            this.c1 = c2;
            this.c2 = c1;
        }
    }

    public boolean is(Chromosome a, Chromosome b) {
        if (c1.getIndex() == a.getIndex() && c2.getIndex() == b.getIndex()) return true;
        return c2.getIndex() == a.getIndex() && c1.getIndex() == b.getIndex();
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj instanceof InterChromosomeRegion) {
            InterChromosomeRegion o = (InterChromosomeRegion) obj;
            return o.is(c1, c2);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(c1, c2);
    }
}
