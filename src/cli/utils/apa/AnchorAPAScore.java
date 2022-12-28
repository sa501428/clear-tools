package cli.utils.apa;

import javastraw.reader.basics.Chromosome;

public class AnchorAPAScore {

    private final Chromosome chromosome;
    private final int gx1, gx2;
    private final String id;
    private final float score;
    private final boolean isUpstream;

    public AnchorAPAScore(Chromosome chromosome, int resolution,
                          int x, int width, String id, float score, boolean isUpstream) {
        this.chromosome = chromosome;
        this.gx1 = (x * resolution) - width;
        this.gx2 = gx1 + (2 * width);
        this.id = id;
        this.score = score;
        this.isUpstream = isUpstream;
    }

    public String getLineForBedgraphFile() {
        return chromosome.getName() + "\t" + gx1 + "\t" + gx2 + "\t" + score;
    }

    public String getLineForBEDFile() {
        return chromosome.getName() + "\t" + gx1 + "\t" + gx2 + "\t" +
                id + "\t" + score + "\t" + (isUpstream ? "+" : "-") +
                gx1 + "\t" + gx2 + "\t" + getColor();
    }

    public String getColor() {
        if (isUpstream) {
            return "0,255,0";
        } else {
            return "255,0,0";
        }
    }
}
