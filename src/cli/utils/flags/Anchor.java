package cli.utils.flags;

import javastraw.feature1D.Feature1D;

public class Anchor extends Feature1D {

    private final String chromName;
    private final int chromIndex, start, end, mid;

    public Anchor(String chromName, int chromIndex, int start1, int end1) {
        this.chromName = chromName;
        this.chromIndex = chromIndex;
        this.start = start1;
        this.end = end1;
        this.mid = (start1 + end1) / 2;
    }

    @Override
    public String getKey() {
        return "" + chromIndex;
    }

    @Override
    public Feature1D deepClone() {
        return new Anchor(chromName, chromIndex, start, end);
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getMid() {
        return mid;
    }
}
