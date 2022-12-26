package cli.utils.data;

public class BoundsInfo {

    private final int binYStart;
    private final int binYEnd;
    private final int binXStart;

    public BoundsInfo(int binYStart, int binYEnd, int binXStart) {
        this.binYStart = binYStart;
        this.binYEnd = binYEnd;
        this.binXStart = binXStart;
    }

    public int getBinYStart() {
        return binYStart;
    }

    public int getBinYEnd() {
        return binYEnd;
    }

    public int getBinXStart() {
        return binXStart;
    }

    public boolean contains(int binY) {
        return binY >= binYStart && binY < binYEnd;
    }
}
