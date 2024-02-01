package cli.utils.data;

public class Bounds2DInfo {

    private final int binXStart;
    private final int binXEnd;
    private final int binYStart;
    private final int binYEnd;


    public Bounds2DInfo(int binXStart, int binXEnd, int binYStart, int binYEnd) {
        this.binXStart = binXStart;
        this.binXEnd = binXEnd;
        this.binYStart = binYStart;
        this.binYEnd = binYEnd;
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

    public boolean contains(int binX, int binY) {
        return binX >= binXStart && binX < binXEnd && binY >= binYStart && binY < binYEnd;
    }
}
