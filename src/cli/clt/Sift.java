package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import cli.utils.sift.SiftUtils;
import cli.utils.sift.SimpleLocation;
import cli.utils.sift.Z4Scores;
import cli.utils.sift.ZScores;
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
        Feature2DList refinedLoops = siftThroughCalls(ds);
        refinedLoops.exportFeatureList(new File(args[2] + ".sift.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("sift complete");
    }

    private static Set<ContactRecord> getHiResExtremePixels(MatrixZoomData zd, int maxBin, int minBin) {

        int maxCompressedBin = logp1i(maxBin) + 1;
        int minCompressedBin = logp1i(minBin);

        ZScores zScores = getZscores(zd, maxCompressedBin, true);

        Set<ContactRecord> records = new HashSet<>();
        for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist > minCompressedBin && dist < maxCompressedBin) {
                    float zscore = zScores.getZscore(dist, logp1(cr.getCounts()));
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

        int maxCompressedBin = logp1i(maxBin) + 1;
        int minCompressedBin = logp1i(minBin);

        double[] nvVC = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), VC).getData().getValues().get(0);
        double[] nvVCSqrt = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), VC_SQRT).getData().getValues().get(0);
        double[] nvSCALE = ds.getNormalizationVector(chrIdx, new HiCZoom(resolution), SCALE).getData().getValues().get(0);

        Z4Scores zScores = getZ4scores(zd, maxCompressedBin, nvVC, nvVCSqrt, nvSCALE);
        Set<SimpleLocation> records = new HashSet<>();
        Iterator<ContactRecord> it = zd.getDirectIterator();

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist > minCompressedBin && dist < maxCompressedBin) {
                    double denomVC = nvVC[cr.getBinX()] * nvVC[cr.getBinY()];
                    double denomVCSqrt = nvVCSqrt[cr.getBinX()] * nvVCSqrt[cr.getBinY()];
                    double denomScale = nvSCALE[cr.getBinX()] * nvSCALE[cr.getBinY()];

                    if (denomVC > 0 && denomVCSqrt > 0 && denomScale > 0) {
                        double valVC = (cr.getCounts() / denomVC);
                        double valVCSqrt = (cr.getCounts() / denomVCSqrt);
                        double valScale = (cr.getCounts() / denomScale);
                        if (valVC > 1 && valVCSqrt > 1 && valScale > 1) {
                            double raw = logp1(cr.getCounts());
                            valVC = logp1(valVC);
                            valVCSqrt = logp1(valVCSqrt);
                            valScale = logp1(valScale);

                            if (zScores.passesAllZscores(dist, LOWRES_ZSCORE_CUTOFF,
                                    raw, valVC, valVCSqrt, valScale)) {
                                records.add(new SimpleLocation(cr));
                            }
                        }
                    }
                }
            }
        }

        return records;
    }

    private static ZScores getZscores(MatrixZoomData zd, int length, boolean useNone) {
        WelfordStats stats = new WelfordStats(length);

        Iterator<ContactRecord> it;
        if (useNone) {
            it = zd.getDirectIterator();
        } else {
            it = zd.getNormalizedIterator(SCALE);
        }

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < length) {
                    stats.addValue(dist, logp1(cr.getCounts()));
                }
            }
        }
        return stats.getZscores();
    }

    private static Z4Scores getZ4scores(MatrixZoomData zd, int length,
                                        double[] nvVC, double[] nvVCSqrt, double[] nvSCALE) {
        WelfordStats rawStats = new WelfordStats(length);
        WelfordStats vcStats = new WelfordStats(length);
        WelfordStats vcSqrtStats = new WelfordStats(length);
        WelfordStats scaleStats = new WelfordStats(length);
        Iterator<ContactRecord> it = zd.getDirectIterator();

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist < length) {
                    rawStats.addValue(dist, logp1(cr.getCounts()));
                    populateNormedStats(cr, nvVC, vcStats, dist);
                    populateNormedStats(cr, nvVCSqrt, vcSqrtStats, dist);
                    populateNormedStats(cr, nvSCALE, scaleStats, dist);
                }
            }
        }
        return new Z4Scores(rawStats, vcStats, vcSqrtStats, scaleStats);
    }

    private static void populateNormedStats(ContactRecord cr, double[] norm, WelfordStats stats, int dist) {
        double denom = norm[cr.getBinX()] * norm[cr.getBinY()];
        if (denom > 0) {
            double val = cr.getCounts() / denom;
            if (val > 1) {
                stats.addValue(dist, logp1(val));
            }
        }
    }

    private static double logp1(double x) {
        return Math.log(1 + x);
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
     * 13 - just the hires calls
     */
    private Feature2DList siftThroughCalls(Dataset ds) {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList output = new Feature2DList();
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chrom, chrom);

            if (matrix != null) {
                int hires = 200;
                MatrixZoomData zd2 = matrix.getZoomData(new HiCZoom(hires));
                System.out.println("Start HiRes pass (" + hires + ")");
                Set<ContactRecord> initialPoints = getHiResExtremePixels(zd2, MAX_DIST / hires, MIN_DIST / hires);
                System.out.println("HiRes pass done (" + hires + ")");

                System.out.println("Num initial loops " + initialPoints.size());
                /*
                for (int lowRes : new int[]{5000}) { // 1000, 2000,
                    MatrixZoomData zd1 = matrix.getZoomData(new HiCZoom(lowRes));
                    System.out.println("Start LowRes pass (" + lowRes + ")");
                    Set<SimpleLocation> enrichedRegions = getExtremeLocations(ds, chrom.getIndex(), lowRes,
                            zd1, MAX_DIST / lowRes, MIN_DIST / lowRes);

                    filterOutByOverlap(initialPoints, enrichedRegions, lowRes / hires);
                    enrichedRegions.clear();
                    System.out.println("LowRes pass done (" + lowRes + ")");
                }
                */

                filterOutByOverlap(initialPoints, 2000 / hires);
                System.out.println("Num initial loops after filter1 " + initialPoints.size());

                matrix.clearCache();

                SiftUtils.coalesceAndRetainCentroids(initialPoints, hires, 5000);
                System.out.println("Num initial loops after filter2 " + initialPoints.size());

                output.addByKey(Feature2DList.getKey(chrom, chrom), convertToFeature2Ds(initialPoints,
                        chrom, chrom, hires));
            }
        }

        return output;
    }

    private void filterOutByOverlap(Set<ContactRecord> initialPoints, int window) {
        Set<ContactRecord> toRemove = new HashSet<>();
        Map<SimpleLocation, List<ContactRecord>> locationMap = new HashMap<>();
        for (ContactRecord cr : initialPoints) {
            SimpleLocation region = new SimpleLocation(cr.getBinX() / window, cr.getBinY() / window);
            if (locationMap.containsKey(region)) {
                locationMap.get(region).add(cr);
            } else {
                List<ContactRecord> values = new ArrayList<>();
                values.add(cr);
                locationMap.put(region, values);
            }
        }

        nonMaxSuppressionInGroup(initialPoints, toRemove, locationMap);
    }

    private void nonMaxSuppressionInGroup(Set<ContactRecord> initialPoints, Set<ContactRecord> toRemove, Map<SimpleLocation, List<ContactRecord>> locationMap) {
        for (List<ContactRecord> overlaps : locationMap.values()) {
            if (overlaps.size() > 1) {
                float max = getMax(overlaps);
                for (ContactRecord cr2 : overlaps) {
                    if (cr2.getCounts() < max) {
                        toRemove.add(cr2);
                    }
                }
            }
        }
        initialPoints.removeAll(toRemove);
        locationMap.clear();
    }

    private float getMax(List<ContactRecord> records) {
        float maxVal = records.get(0).getCounts();
        for (ContactRecord cr : records) {
            if (cr.getCounts() > maxVal) {
                maxVal = cr.getCounts();
            }
        }
        return maxVal;
    }

    private boolean inRegions(ContactRecord cr, Set<SimpleLocation> regions, int scalar) {
        SimpleLocation region = new SimpleLocation(cr.getBinX() / scalar, cr.getBinY() / scalar);
        return regions.contains(region);
    }

    private static int logp1i(int x) {
        return (int) Math.floor(Math.log(1 + x));
    }

    private static double logp1(float x) {
        return Math.log(1 + x);
    }

    private void filterOutByOverlap(Set<ContactRecord> initialPoints, Set<SimpleLocation> regions, int scalar) {
        Set<ContactRecord> toRemove = new HashSet<>();
        Map<SimpleLocation, List<ContactRecord>> locationMap = new HashMap<>();
        for (ContactRecord cr : initialPoints) {
            SimpleLocation region = new SimpleLocation(cr.getBinX() / scalar, cr.getBinY() / scalar);
            if (regions.contains(region)) {
                if (locationMap.containsKey(region)) {
                    locationMap.get(region).add(cr);
                } else {
                    List<ContactRecord> values = new ArrayList<>();
                    values.add(cr);
                    locationMap.put(region, values);
                }
            } else {
                toRemove.add(cr);
            }
        }

        nonMaxSuppressionInGroup(initialPoints, toRemove, locationMap);
    }
}
