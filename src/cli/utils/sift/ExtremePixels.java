package cli.utils.sift;

import cli.clt.Sift;
import cli.utils.expected.ExpectedModel;
import cli.utils.expected.ExpectedUtils;
import cli.utils.sift.collapse.CentroidCollapser;
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

    public static Set<ContactRecord> getExtremePixelsForResolution(Dataset ds, MatrixZoomData zd, Chromosome chrom,
                                                                   int res, NormalizationType norm,
                                                                   int maxBin, int minBin, ExpectedModel poly) {
        Set<ContactRecord> enrichedRegions = ExtremePixels.getExtremeLocations(ds, chrom, res,
                zd, maxBin, minBin, norm, poly);
        int radius = Math.max(Sift.MIN_RADIUS_0 / res, 2);
        return CentroidCollapser.coalesce(enrichedRegions, radius, radius);
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

        Set<ContactRecord> extremes = new HashSet<>();
        Iterator<ContactRecord> it = ExpectedUtils.getIterator(zd, norm);
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
        return nv[cr.getBinX()] >= Sift.MIN_NORM && nv[cr.getBinY()] >= Sift.MIN_NORM;
    }
}
