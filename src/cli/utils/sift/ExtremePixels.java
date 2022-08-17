package cli.utils.sift;

import cli.utils.ExpectedUtils;
import cli.utils.expected.LogExpectedModel;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class ExtremePixels {

    private static final float CONTACT_ZSCORE_CUTOFF = 1.5f;
    private static final int MAX_DIST = 10000000;
    private static final int MIN_DIST = 10000;


    public static Set<SimpleLocation> getExtremePixelsForResolution(Dataset ds, MatrixZoomData zd, Chromosome chrom,
                                                                    int res, NormalizationType norm) {
        Set<ContactRecord> enrichedRegions = ExtremePixels.getExtremeLocations(ds, chrom, res,
                zd, MAX_DIST / res, MIN_DIST / res, norm);

        return FeatureUtils.toLocationsAndClear(coalescePixelsToCentroid(enrichedRegions));
    }

    public static Set<ContactRecord> getExtremeLocations(Dataset ds, Chromosome chromosome, int resolution,
                                                         MatrixZoomData zd, int maxBin, int minBin,
                                                         NormalizationType norm) {
        int chrIdx = chromosome.getIndex();
        double[] nv;
        try {
            nv = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), norm).getData().getValues().get(0);
        } catch (Exception e) {
            System.err.println("No norm vector found for " + chromosome.getName() + " resolution " + resolution);
            System.exit(8);
            return new HashSet<>();
        }

        List<ContactRecord> records = populateRecordsInRange(zd, norm, minBin, maxBin, nv, 1);
        System.out.println(resolution + " - num records: " + records.size());

        LogExpectedModel model = new LogExpectedModel(records, maxBin);
        ZScores zScores = model.getZscores();

        Set<ContactRecord> extremes = new HashSet<>();
        for (ContactRecord cr : records) {
            double percentContact = model.getPercentContact(cr);
            if (isReasonableEnrichment(percentContact)) {
                int dist = model.logp1i(ExpectedUtils.getDist(cr));
                double val = LogExpectedModel.logp1(cr.getCounts());
                if (zScores.getZscore(dist, val) > CONTACT_ZSCORE_CUTOFF) {
                    extremes.add(cr);
                }
            }
        }
        System.out.println(resolution + " - num extremes: " + extremes.size());

        records.clear();

        return extremes;
    }

    private static List<ContactRecord> populateRecordsInRange(MatrixZoomData zd, NormalizationType norm,
                                                              int minBin, int maxBin, double[] nv, int minVal) {
        List<ContactRecord> records = new LinkedList<>();
        Iterator<ContactRecord> it = getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > minVal) {
                int dist = ExpectedUtils.getDist(cr);
                if (dist > minBin && dist < maxBin) {
                    if (nv[cr.getBinX()] > 1 && nv[cr.getBinY()] > 1) {
                        records.add(cr);
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
        return val > 0.005 && val < 0.5;
    }

    public static Set<ContactRecord> coalescePixelsToCentroid(Set<ContactRecord> regions) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        LinkedList<ContactRecord> featureLL = new LinkedList<>(regions);
        featureLL.sort((o1, o2) -> Float.compare(-o1.getCounts(), -o2.getCounts()));

        Set<ContactRecord> coalesced = new HashSet<>();
        while (!featureLL.isEmpty()) {
            ContactRecord pixel = featureLL.pollFirst();
            if (pixel != null) {
                coalesced.add(pixel);
                featureLL.remove(pixel);

                int buffer = 1;

                int binX0 = pixel.getBinX() - buffer;
                int binY0 = pixel.getBinY() - buffer;
                int binX1 = pixel.getBinX() + buffer + 1;
                int binY1 = pixel.getBinY() + buffer + 1;

                int prevSize = 0;
                Set<ContactRecord> pixelList = new HashSet<>();
                pixelList.add(pixel);

                while (prevSize != pixelList.size()) {
                    prevSize = pixelList.size();
                    for (ContactRecord px : featureLL) {
                        if (contains(px, binX0, binY0, binX1, binY1)) {
                            pixelList.add(px);
                            binX0 = Math.min(binX0, px.getBinX() - buffer);
                            binY0 = Math.min(binY0, px.getBinY() - buffer);
                            binX1 = Math.max(binX1, px.getBinX() + buffer + 1);
                            binY1 = Math.max(binY1, px.getBinY() + buffer + 1);
                        }
                    }
                    featureLL.removeAll(pixelList);
                }
            }
        }
        return coalesced;
    }

    private static boolean contains(ContactRecord px, int binX0, int binY0, int binX1, int binY1) {
        return binX0 <= px.getBinX() && binX1 > px.getBinX() && binY0 <= px.getBinY() && binY1 > px.getBinY();
    }
}
