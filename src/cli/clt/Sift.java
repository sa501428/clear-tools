package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import cli.utils.expected.LogExpectedModel;
import cli.utils.sift.*;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

import java.awt.*;
import java.io.File;
import java.util.List;
import java.util.*;


public class Sift {
    private final int MAX_DIST = 10000000;

    private static final int HIRES_ZSCORE_CUTOFF = 2;
    private static final int LOWRES_ZSCORE_CUTOFF = 2;
    private static final NormalizationType SCALE = NormalizationHandler.SCALE;
    private static final NormalizationType VC = NormalizationHandler.VC;
    private static final NormalizationType VC_SQRT = NormalizationHandler.VC_SQRT;
    private final int MIN_DIST = 10000;
    private final int window = 5;

    public Sift(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            Main.printGeneralUsageAndExit(5);
        }
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1],
                false, false, false);
        //File outFolder = UNIXTools.makeDir(new File(args[2]));
        int hires = parser.getResolutionOption(200);
        int lowres = parser.getLowResolutionOption(5000);
        BoundingBoxWithContacts.width = parser.getWindowSizeOption(5);
        BoundingBoxWithContacts.buffer = 2 * BoundingBoxWithContacts.width;

        BoundingBoxWithContacts.minEnrichment = parser.getMinOption(1.2);
        BoundingBoxWithContacts.maxEnrichment = parser.getMaxOption(50);

        Feature2DList refinedLoops = siftThroughCalls(ds, hires, lowres);
        refinedLoops.exportFeatureList(new File(args[2] + ".sift.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("sift complete");
    }

    private static Set<ContactRecord> getHiResExtremePixels(MatrixZoomData zd, int maxBin, int minBin) {

        int maxCompressedBin = LogExpectedModel.logp1i(maxBin) + 1;
        int minCompressedBin = LogExpectedModel.logp1i(minBin);

        ZScores zScores = getZscores(zd, maxCompressedBin, true);

        Set<ContactRecord> records = new HashSet<>();
        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = LogExpectedModel.logp1i(ExpectedUtils.getDist(cr));
                if (dist > minCompressedBin && dist < maxCompressedBin) {
                    float zscore = zScores.getZscore(dist, LogExpectedModel.logp1(cr.getCounts()));
                    if (zscore > HIRES_ZSCORE_CUTOFF) {
                        records.add(cr);
                    }
                }
            }
        }

        return records;
    }

    private static Set<SimpleLocation> getExtremeLocations(Dataset ds, int chrIdx, int resolution,
                                                           MatrixZoomData zd, int maxBin, int minBin) {

        int maxCompressedBin = LogExpectedModel.logp1i(maxBin) + 1;
        int minCompressedBin = LogExpectedModel.logp1i(minBin);

        //double[] nvVC = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), VC).getData().getValues().get(0);
        //double[] nvVCSqrt = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), VC_SQRT).getData().getValues().get(0);
        double[] nvSCALE = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), SCALE).getData().getValues().get(0);

        //Z4Scores zScores = getZ4scores(zd, maxCompressedBin, nvVC, nvVCSqrt, nvSCALE);

        ZScores zScores = getZscores(zd, maxCompressedBin, false);

        Set<SimpleLocation> records = new HashSet<>();
        Iterator<ContactRecord> it = zd.getDirectIterator();

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = LogExpectedModel.logp1i(ExpectedUtils.getDist(cr));
                if (dist > minCompressedBin && dist < maxCompressedBin) {
                    //double denomVC = nvVC[cr.getBinX()] * nvVC[cr.getBinY()];
                    //double denomVCSqrt = nvVCSqrt[cr.getBinX()] * nvVCSqrt[cr.getBinY()];
                    double denomScale = nvSCALE[cr.getBinX()] * nvSCALE[cr.getBinY()];

                    if (denomScale > 0.75) { // denomVC > 1 && denomVCSqrt > 1 &&
                        //double valVC = (cr.getCounts() / denomVC);
                        //double valVCSqrt = (cr.getCounts() / denomVCSqrt);
                        double valScale = (cr.getCounts() / denomScale);
                        if (valScale > 1) { // valVC > 1 && valVCSqrt > 1 &&
                            //double raw = logp1(cr.getCounts());
                            //valVC = logp1(valVC);
                            //valVCSqrt = logp1(valVCSqrt);
                            valScale = LogExpectedModel.logp1(valScale);

                            // .passesAllZscores(dist, LOWRES_ZSCORE_CUTOFF,
                            //                                    raw, valVC, valVCSqrt, valScale)
                            if (zScores.getZscore(dist, valScale) > LOWRES_ZSCORE_CUTOFF) {
                                records.add(new SimpleLocation(cr));
                            }
                        }
                    }
                }
            }
        }

        return records;
    }

    public static ZScores getZscores(MatrixZoomData zd, int length, boolean useNone) {
        WelfordStats stats = LogExpectedModel.getSummaryStats(zd, length, useNone, 1, SCALE);
        return stats.getZscores();
    }

    public static List<Feature2D> convertToFeature2Ds(Set<ContactRecord> records,
                                                      Chromosome c1, Chromosome c2, int resolution) {
        List<Feature2D> features = new ArrayList<>();
        for (ContactRecord record : records) {
            long start1 = (long) record.getBinX() * resolution;
            long end1 = start1 + resolution;
            long start2 = (long) record.getBinY() * resolution;
            long end2 = start2 + resolution;
            features.add(new Feature2D(Feature2D.FeatureType.PEAK, c1.getName(), start1, end1,
                    c2.getName(), start2, end2, Color.BLACK, new HashMap<>()));
        }
        return features;
    }

    /*
     * Iteration notes
     * 1 - original pass
     * 2 - similar, fix stuff
     * 3 zscore at 3, diff filtering
     * 4 zscore t0 2.5
     * 5 - change zscore to 2
     * 6 - change hires from 100 to 200, and add 500 and 2000 to the resolution runs (already have 1000, 5000 )
     * 7 - remove 500 bp res
     * 8 - coalesce pixels radius 13kb
     * 9 - coalesce pixels radius 5kb
     * 10 - use multiple normalizations *** seems best so far, but pretty conservative
     * 11 - only validate at 5000, not 1000/2000
     * 12 - higher zscore cutoff for 5kb
     * 13 - just the hires calls, collapse at 2k
     * 14 - scale/vc vector filtering, collapse at 5k
     * 15 - restore global max filtering
     * 16 - linear distance for 5k lowres expected
     * 17 - incorporate local enrichment; only use 1 norm for global enrichment at lowres,
     *      remove vector filtering but keep denom filtering
     * 18 - no max neighbor
     * 19 - change to 20% enrichment, not 25%; also put prefiltering earlier
     * 20 - change low res zscore to 1
     * 21 - change low res zscore to 1.5
     * 22 - same as 19 (low res zscore 2)
     * 23 - 2x min local enrichment
     */
    private Feature2DList siftThroughCalls(Dataset ds, int hiRes, int lowRes) {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList output = new Feature2DList();
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chrom, chrom);

            if (matrix != null) {
                MatrixZoomData zdHigh = matrix.getZoomData(new HiCZoom(hiRes));
                System.out.println("Start HiRes pass (" + hiRes + ")");
                Set<ContactRecord> initialPoints = getHiResExtremePixels(zdHigh, MAX_DIST / hiRes, MIN_DIST / hiRes);
                System.out.println("HiRes pass done (" + hiRes + ")");
                matrix.clearCacheForZoom(new HiCZoom(hiRes));

                System.out.println("Num initial loops " + initialPoints.size());

                NMSUtils.filterOutByOverlap(initialPoints, lowRes / hiRes);
                System.out.println("Num loops after pre filter (overlaps) " + initialPoints.size());

                MatrixZoomData zdLow = matrix.getZoomData(new HiCZoom(lowRes));
                System.out.println("Start LowRes pass (" + lowRes + ")");
                Set<SimpleLocation> enrichedRegions = getExtremeLocations(ds, chrom.getIndex(), lowRes,
                        zdLow, MAX_DIST / lowRes, MIN_DIST / lowRes);

                NMSUtils.filterOutByOverlap(initialPoints, enrichedRegions, lowRes / hiRes);
                enrichedRegions.clear();
                System.out.println("LowRes pass done (" + lowRes + ")");

                System.out.println("Num loops after low res global filter " + initialPoints.size());

                EnrichmentChecker.filterOutIfNotLocalMax(zdLow, initialPoints, lowRes / hiRes, SCALE);
                matrix.clearCache();

                System.out.println("Num loops after low res local filter " + initialPoints.size());

                SiftUtils.coalesceAndRetainCentroids(initialPoints, hiRes, lowRes);
                System.out.println("Num loops after filter3 " + initialPoints.size());

                output.addByKey(Feature2DList.getKey(chrom, chrom), convertToFeature2Ds(initialPoints,
                        chrom, chrom, hiRes));
            }
        }

        return output;
    }
}
