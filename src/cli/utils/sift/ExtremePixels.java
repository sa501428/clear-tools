package cli.utils.sift;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.expected.LogExpectedModel;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class ExtremePixels {

    private static final int HIRES_ZSCORE_CUTOFF = 2;
    private static final float LOWRES_ZSCORE_CUTOFF = 1.5f;
    private static final NormalizationType SCALE = NormalizationHandler.SCALE;
    private static final int MAX_DIST = 10000000;
    private static final int MIN_DIST = 10000;


    public static Set<ContactRecord> getExtremePixelsForResolution(Set<ContactRecord> initialPoints, Dataset ds,
                                                                   MatrixZoomData zdLow, Chromosome chrom,
                                                                   int hiRes, int lowRes) {
        Set<ContactRecord> points = filterInStagesFrom(initialPoints, hiRes, lowRes);
        if (Main.printVerboseComments) {
            System.out.println("Num loops after pre filter (overlaps) " + points.size() +
                    "\nStart LowRes pass (" + lowRes + ")");
        }

        Set<SimpleLocation> enrichedRegions = ExtremePixels.getExtremeLocations(ds, chrom, lowRes,
                zdLow, MAX_DIST / lowRes, MIN_DIST / lowRes); // 2000, 8

        NMSUtils.filterOutByOverlap(points, enrichedRegions, lowRes / hiRes);
        if (Main.printVerboseComments)
            System.out.println("Num loops after low res global filter " + points.size());

        enrichedRegions.clear();

        EnrichmentChecker.filterOutIfNotLocalMax(zdLow, points, lowRes / hiRes, SCALE);
        if (Main.printVerboseComments) {
            System.out.println("Num loops after low res local filter " + points.size() +
                    "\nLowRes pass done (" + lowRes + ")\n");
        }
        return points;
    }

    public static Set<ContactRecord> getHiResExtremePixels(MatrixZoomData zd, int hiRes) {

        int maxBin = MAX_DIST / hiRes;
        int minBin = MIN_DIST / hiRes;

        LogExpectedModel model = new LogExpectedModel(zd, SCALE, maxBin, true, 1);
        ZScores zScores = model.getZscores();

        Set<ContactRecord> records = new HashSet<>();
        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = ExpectedUtils.getDist(cr);
                if (dist > minBin && dist < maxBin) {
                    dist = model.logp1i(dist);
                    float zscore = zScores.getZscore(dist, LogExpectedModel.logp1(cr.getCounts()));
                    if (zscore > HIRES_ZSCORE_CUTOFF) {
                        records.add(cr);
                    }
                }
            }
        }

        return records;
    }

    public static Set<SimpleLocation> getExtremeLocations(Dataset ds, Chromosome chromosome, int resolution,
                                                          MatrixZoomData zd, int maxBin, int minBin) {
        int chrIdx = chromosome.getIndex();
        double[] nvSCALE;
        try {
            nvSCALE = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), SCALE).getData().getValues().get(0);
        } catch (Exception e) {
            System.err.println("No norm vector found for " + chromosome.getName() + " resolution " + resolution);
            System.exit(8);
            return new HashSet<>();
        }

        LogExpectedModel model = new LogExpectedModel(zd, SCALE, maxBin, false, 1);
        ZScores zScores = model.getZscores();

        Set<SimpleLocation> records = new HashSet<>();
        Iterator<ContactRecord> it = zd.getDirectIterator();

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = ExpectedUtils.getDist(cr);
                if (dist > minBin && dist < maxBin) {
                    double percentContact = model.getPercentContact(dist, cr.getCounts());
                    if (isReasonableEnrichment(percentContact)) {
                        double nv1 = nvSCALE[cr.getBinX()];
                        double nv2 = nvSCALE[cr.getBinY()];
                        double denomScale = nv1 * nv2;
                        if (denomScale > 1 && nv1 > 1 && nv2 > 1) {
                            double valScale = (cr.getCounts() / denomScale);
                            if (valScale > 1) {
                                dist = model.logp1i(dist);
                                valScale = LogExpectedModel.logp1(valScale);
                                if (zScores.getZscore(dist, valScale) > LOWRES_ZSCORE_CUTOFF) {
                                    records.add(new SimpleLocation(cr));
                                }
                            }
                        }
                    }
                }
            }
        }

        return records;
    }

    public static boolean isReasonableEnrichment(double val) {
        return val > 0.01 && val < 0.5;
    }

    private static Set<ContactRecord> filterInStagesFrom(Set<ContactRecord> initialPoints, int hiRes, int lowRes) {
        Set<ContactRecord> points = new HashSet<>(initialPoints);
        int factor = lowRes / 500;
        for (int k = 1; k <= factor; k++) {
            int res = 500 * k;
            NMSUtils.filterOutByOverlap(points, res / hiRes);
        }
        return points;
    }
}
