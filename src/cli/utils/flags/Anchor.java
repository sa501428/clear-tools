package cli.utils.flags;

import javastraw.feature1D.Feature1D;

import java.util.Objects;

public class Anchor extends Feature1D {

    protected final String chromName;
    protected final long start, end, mid;
    protected final int chrIndex;

    public Anchor(String chromName, long start1, long end1, int chrIndex) {
        if (chromName.startsWith("chr")) {
            this.chromName = chromName;
        } else {
            this.chromName = "chr" + chromName;
        }
        this.start = start1;
        this.end = end1;
        this.mid = (start1 + end1) / 2;
        this.chrIndex = chrIndex;
    }

    @Override
    public String getKey() {
        return "" + chrIndex;
    }

    public long getStart() {
        return start;
    }

    public long getEnd() {
        return end;
    }

    public long getMid() {
        return mid;
    }

    @Override
    public Feature1D deepClone() {
        return new Anchor(chromName, start, end, chrIndex);
    }

    @Override
    public String toString() {
        return chromName + "\t" + start + "\t" + end;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Anchor) {
            Anchor a = (Anchor) obj;
            return a.chromName.equals(chromName) && a.start == start && a.end == end;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromName, start, end);
    }

    public int getWidth() {
        return (int) (end - start);
    }
}
