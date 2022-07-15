package cli.utils.sift;

import javastraw.reader.block.ContactRecord;

import java.util.Objects;

public class SimpleLocation {
    private final int binX, binY;

    public SimpleLocation(int binX, int binY) {
        this.binX = binX;
        this.binY = binY;
    }

    public SimpleLocation(ContactRecord cr) {
        this.binX = cr.getBinX();
        this.binY = cr.getBinY();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        SimpleLocation o2 = (SimpleLocation) o;
        return binX == o2.binX && binY == o2.binY;
    }

    @Override
    public int hashCode() {
        return Objects.hash(binX, binY);
    }

    public int getBinX() {
        return this.binX;
    }

    public int getBinY() {
        return this.binY;
    }
}
