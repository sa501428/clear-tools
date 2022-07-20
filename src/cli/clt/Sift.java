package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
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

        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);

        Feature2DList[] refinedLoops = siftThroughCalls(ds);
        //refinedLoops[0].exportFeatureList(new File(args[2] + ".c1.sift.bedpe"), false, Feature2DList.ListFormat.NA);
        refinedLoops[1].exportFeatureList(new File(args[2] + ".sift.bedpe"), false, Feature2DList.ListFormat.NA);
        //refinedLoops[2].exportFeatureList(new File(args[2] + ".c3.sift.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("sift complete");
    }

    private static Set<ContactRecord> getHiResExtremePixels(MatrixZoomData zd, int maxBin, int minBin) {

        LogExpectedModel model = new LogExpectedModel(zd, SCALE, maxBin, true);
        ZScores zScores = model.getZscores();

        Set<ContactRecord> records = new HashSet<>();
        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = ExpectedUtils.getDist(cr);
                if (dist > minBin && dist < maxBin) {
                    dist = LogExpectedModel.logp1i(dist);
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

        double[] nvSCALE = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), SCALE).getData().getValues().get(0);

        LogExpectedModel model = new LogExpectedModel(zd, SCALE, maxBin, false);
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
                        double denomScale = nvSCALE[cr.getBinX()] * nvSCALE[cr.getBinY()];
                        if (denomScale > 1) { // denomVC > 1 && denomVCSqrt > 1 &&
                            double valScale = (cr.getCounts() / denomScale);
                            if (valScale > 1) {
                                dist = LogExpectedModel.logp1i(dist);
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

    private static boolean isReasonableEnrichment(double val) {
        return val > 0.01 && val < 0.5;
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
     *
     * 30 - try with 1k and 5k, with pc filter
     * 31 - restrict 1k partially
     * 32 - remove 1k
     * 33 - simplify local filtering (medium status)
     * 34 - add back 1kb
     * 40 - adjust neighbor windows to 1 and 2 (good)
     * 50 - restructure (worse??)
     * 51 - include 2kb, 10kb
     * 52 - higher cutoff
     *
     * 55 - try stepwise filtering, full range
     */
    private Feature2DList[] siftThroughCalls(Dataset ds) {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList[] outputs = new Feature2DList[]{new Feature2DList(), new Feature2DList(), new Feature2DList()};

        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chrom, chrom);

            if (matrix != null) {
                int hiRes = 200;
                MatrixZoomData zdHigh = matrix.getZoomData(new HiCZoom(hiRes));
                System.out.println("Start HiRes pass (" + hiRes + ")");
                Set<ContactRecord> initialPoints = getHiResExtremePixels(zdHigh, MAX_DIST / hiRes, MIN_DIST / hiRes);
                System.out.println("HiRes pass done (" + hiRes + ")");
                matrix.clearCacheForZoom(new HiCZoom(hiRes));
                System.out.println("Num initial loops " + initialPoints.size());

                Map<ContactRecord, Short> countsForRecord = new HashMap<>();
                for (int lowRes : new int[]{1000, 2000, 5000, 10000}) { // 1000,
                    Set<ContactRecord> points = filterInStagesFrom(initialPoints, hiRes, lowRes);
                    System.out.println("Num loops after pre filter (overlaps) " + points.size() +
                            "\nStart LowRes pass (" + lowRes + ")");

                    MatrixZoomData zdLow = matrix.getZoomData(new HiCZoom(lowRes));

                    Set<SimpleLocation> enrichedRegions = getExtremeLocations(ds, chrom.getIndex(), lowRes,
                            zdLow, MAX_DIST / lowRes, MIN_DIST / lowRes); // 2000, 8

                    NMSUtils.filterOutByOverlap(points, enrichedRegions, lowRes / hiRes);
                    System.out.println("Num loops after low res global filter " + points.size());

                    enrichedRegions.clear();

                    EnrichmentChecker.filterOutIfNotLocalMax(zdLow, points, lowRes / hiRes, SCALE);
                    System.out.println("Num loops after low res local filter " + points.size() +
                            "\nLowRes pass done (" + lowRes + ")\n");

                    matrix.clearCacheForZoom(new HiCZoom(lowRes));

                    for (ContactRecord record : points) {
                        if (countsForRecord.containsKey(record)) {
                            countsForRecord.put(record, (short) (countsForRecord.get(record) + 1));
                        } else {
                            countsForRecord.put(record, (short) 1);
                        }
                    }
                    points.clear();
                }
                matrix.clearCache();
                initialPoints.clear();

                for (int c = 0; c < outputs.length; c++) {
                    Set<ContactRecord> finalPoints = getPointWithMoreThan(countsForRecord, c + 1);
                    SiftUtils.coalesceAndRetainCentroids(finalPoints, hiRes, 5000);
                    System.out.println("Num loops after final filter[" + c + "] " + finalPoints.size());
                    outputs[c].addByKey(Feature2DList.getKey(chrom, chrom), convertToFeature2Ds(finalPoints,
                            chrom, chrom, hiRes));
                    finalPoints.clear();
                }
                countsForRecord.clear();
            }
        }

        return outputs;
    }

    private Set<ContactRecord> filterInStagesFrom(Set<ContactRecord> initialPoints, int hiRes, int lowRes) {
        Set<ContactRecord> points = new HashSet<>(initialPoints);
        int factor = lowRes / 500;
        for (int k = 1; k <= factor; k++) {
            int res = 500 * k;
            NMSUtils.filterOutByOverlap(points, res / hiRes);
        }
        return points;
    }

    private Set<ContactRecord> getPointWithMoreThan(Map<ContactRecord, Short> countsForRecord, int cutoff) {
        Set<ContactRecord> finalSet = new HashSet<>();
        for (ContactRecord record : countsForRecord.keySet()) {
            if (countsForRecord.get(record) > cutoff) {
                finalSet.add(record);
            }
        }
        return finalSet;
    }
}
