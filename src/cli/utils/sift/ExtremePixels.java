package cli.utils.sift;

import cli.utils.expected.ExpectedModel;
import cli.utils.expected.ExpectedUtils;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class ExtremePixels {

    public static Set<SimpleLocation> getExtremePixelsForResolution(Dataset ds, MatrixZoomData zd, Chromosome chrom,
                                                                    int res, NormalizationType norm,
                                                                    int maxBin, int minBin, ExpectedModel poly) {
        Set<ContactRecord> enrichedRegions = ExtremePixels.getExtremeLocations(ds, chrom, res,
                zd, maxBin, minBin, norm, poly);

        return FeatureUtils.toLocationsAndClear(coalescePixelsToCentroid(enrichedRegions));
    }

    public static Set<ContactRecord> getExtremeLocations(Dataset ds, Chromosome chromosome, int resolution,
                                                         MatrixZoomData zd, int maxBin, int minBin,
                                                         NormalizationType norm, ExpectedModel poly) {
        int chrIdx = chromosome.getIndex();
        double[] nv;
        try {
            nv = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), norm).getData().getValues().get(0);
        } catch (Exception e) {
            System.err.println("No norm vector found for " + chromosome.getName() + " resolution " + resolution);
            System.exit(8);
            return new HashSet<>();
        }

        int minVal = 1;
        //List<ContactRecord> records = populateRecordsInRange(zd, norm, maxBin, 1, minBin, nv);

        Set<ContactRecord> extremes = new HashSet<>();
        Iterator<ContactRecord> it = getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > minVal && isReasonableNorm(cr, nv)) {
                int dist = ExpectedUtils.getDist(cr);
                if (dist > minBin && dist < maxBin) {
                    if (poly.isReasonableEnrichment(cr) && poly.isReasonablePercentContact(cr)) {
                        extremes.add(cr);
                    }
                }
            }
        }

        System.out.println(resolution + " - num extremes: " + extremes.size());
        return extremes;
    }

    private static boolean isReasonableNorm(ContactRecord cr, double[] nv) {
        return nv[cr.getBinX()] > 1 && nv[cr.getBinY()] > 1;
    }

    private static List<ContactRecord> populateRecordsInRange(MatrixZoomData zd, NormalizationType norm,
                                                              int maxBin, int minVal, int minBin, double[] nv) {
        List<ContactRecord> records = new LinkedList<>();
        Iterator<ContactRecord> it = getIterator(zd, norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > minVal) {
                if (isReasonableNorm(cr, nv)) {
                    int dist = ExpectedUtils.getDist(cr);
                    if (dist > minBin && dist < maxBin) {
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

    public static Set<ContactRecord> coalescePixelsToCentroid(Set<ContactRecord> regions) {
        // HashSet intermediate for removing duplicates
        // LinkedList used so that we can pop out highest obs values
        Map<SimpleLocation, LinkedList<ContactRecord>> map = NMSUtils.groupNearbyRecords(regions, 250);
        Set<ContactRecord> coalesced = new HashSet<>();

        for (LinkedList<ContactRecord> records : map.values()) {
            records.sort((o1, o2) -> Float.compare(-o1.getCounts(), -o2.getCounts()));

            while (!records.isEmpty()) {
                ContactRecord pixel = records.pollFirst();
                if (pixel != null) {
                    records.remove(pixel);

                    int buffer = 2;

                    int binX0 = pixel.getBinX() - buffer;
                    int binY0 = pixel.getBinY() - buffer;
                    int binX1 = pixel.getBinX() + buffer + 1;
                    int binY1 = pixel.getBinY() + buffer + 1;

                    int prevSize = 0;
                    Set<ContactRecord> pixelList = new HashSet<>();
                    pixelList.add(pixel);

                    while (prevSize != pixelList.size()) {
                        prevSize = pixelList.size();
                        for (ContactRecord px : records) {
                            if (contains(px, binX0, binY0, binX1, binY1)) {
                                pixelList.add(px);
                                binX0 = Math.min(binX0, px.getBinX() - buffer);
                                binY0 = Math.min(binY0, px.getBinY() - buffer);
                                binX1 = Math.max(binX1, px.getBinX() + buffer + 1);
                                binY1 = Math.max(binY1, px.getBinY() + buffer + 1);
                            }
                        }
                        records.removeAll(pixelList);
                    }

                    if (pixelList.size() > buffer) {
                        coalesced.add(pixel);
                    }
                }
            }
        }
        return coalesced;
    }

    private static boolean contains(ContactRecord px, int binX0, int binY0, int binX1, int binY1) {
        return binX0 <= px.getBinX() && binX1 > px.getBinX() && binY0 <= px.getBinY() && binY1 > px.getBinY();
    }
}
