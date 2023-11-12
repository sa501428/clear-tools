package cli.utils.anchors;

import cli.utils.flags.Anchor;
import javastraw.reader.basics.Chromosome;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class AnchorPeakFinder {

    public static Set<Anchor> getPeaks(int resolution, Map<Chromosome, double[]> allUpStreamOEProd,
                                       Map<Chromosome, double[]> allDownStreamOEProd,
                                       Map<Chromosome, double[]> allBothStreamOEProd) {
        Set<Anchor> peaks = new HashSet<>();
        for (Chromosome chrom : allBothStreamOEProd.keySet()) {
            double[] both = allBothStreamOEProd.get(chrom);
            double[] up = allUpStreamOEProd.get(chrom);
            double[] down = allDownStreamOEProd.get(chrom);

            double[] bothSmooth = ConvolutionTools.smooth(both);
            double[] upSmooth = ConvolutionTools.smooth(up);
            double[] downSmooth = ConvolutionTools.smooth(down);

            for (int i = 2; i < both.length - 2; i++) {
                if (elevated5(both, i) && elevated5(bothSmooth, i)) {
                    if (elevated5(upSmooth, i) || elevated5(downSmooth, i)) {
                        if (elevated3(up, i) || elevated3(down, i)) {
                            peaks.add(new Anchor(chrom.getName(),
                                    (long) i * resolution, (long) (i + 1) * resolution,
                                    chrom.getIndex()));
                        }
                    }
                }
            }
        }
        return peaks;
    }

    private static boolean elevated5(double[] data, int i) {
        return data[i] > 1.1 * data[i - 1]
                && data[i] > 1.1 * data[i + 1]
                && data[i - 1] > 1.1 * data[i - 2]
                && data[i + 1] > 1.1 * data[i + 2]
                && data[i + 2] > 1
                && data[i - 2] > 1;
    }

    private static boolean elevated3(double[] data, int i) {
        return data[i] > 1.1 * data[i - 1]
                && data[i] > 1.1 * data[i + 1]
                && data[i + 1] > 1
                && data[i - 1] > 1;
    }
}
