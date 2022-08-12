package cli.utils.sift;

import cli.utils.ExpectedUtils;
import cli.utils.expected.LogExpectedModel;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class ExtremePixels {

    private static final int CONTACT_ZSCORE_CUTOFF = 2;
    private static final int MAX_DIST = 10000000;
    private static final int MIN_DIST = 10000;


    public static Set<SimpleLocation> getExtremePixelsForResolution(Dataset ds, MatrixZoomData zd, Chromosome chrom,
                                                                    int res, NormalizationType norm) {
        Set<ContactRecord> enrichedRegions = ExtremePixels.getExtremeLocations(ds, chrom, res,
                zd, MAX_DIST / res, MIN_DIST / res, norm);

        EnrichmentChecker.filterOutIfNotLocalMax(zd, enrichedRegions, norm);

        Set<SimpleLocation> locations = new HashSet<>();
        for (ContactRecord record : enrichedRegions) {
            locations.add(new SimpleLocation(record));
        }
        enrichedRegions.clear();

        return locations;
    }

    public static Set<ContactRecord> getExtremeLocations(Dataset ds, Chromosome chromosome, int resolution,
                                                         MatrixZoomData zd, int maxBin, int minBin,
                                                         NormalizationType norm) {
        int chrIdx = chromosome.getIndex();
        double[] nvSCALE;
        try {
            nvSCALE = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), norm).getData().getValues().get(0);
        } catch (Exception e) {
            System.err.println("No norm vector found for " + chromosome.getName() + " resolution " + resolution);
            System.exit(8);
            return new HashSet<>();
        }

        LogExpectedModel model = new LogExpectedModel(zd, norm, maxBin, 1);
        ZScores zScores = model.getZscores();

        Set<ContactRecord> records = new HashSet<>();
        Iterator<ContactRecord> it = getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = ExpectedUtils.getDist(cr);
                if (dist > minBin && dist < maxBin) {
                    double percentContact = model.getPercentContact(dist, cr.getCounts());
                    if (isReasonableEnrichment(percentContact)) {
                        double nv1 = nvSCALE[cr.getBinX()];
                        double nv2 = nvSCALE[cr.getBinY()];
                        if (nv1 > 1 && nv2 > 1) {
                            double valScale = (cr.getCounts() / (nv1 * nv2));
                            if (valScale > 1) {
                                dist = model.logp1i(dist);
                                valScale = LogExpectedModel.logp1(valScale);
                                if (zScores.getZscore(dist, valScale) > CONTACT_ZSCORE_CUTOFF) {
                                    records.add(cr);
                                }
                            }
                        }
                    }
                }
            }
        }

        return records;
    }

    public static Iterator<ContactRecord> getIterator(MatrixZoomData zd, NormalizationType norm) {
        if (norm.getLabel().equalsIgnoreCase("none")) {
            return zd.getDirectIterator();
        } else {
            return zd.getNormalizedIterator(norm);
        }
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
