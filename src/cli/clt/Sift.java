package cli.clt;

import cli.Main;
import cli.utils.ExpectedUtils;
import cli.utils.WelfordStats;
import cli.utils.sift.SiftUtils;
import cli.utils.sift.SimpleLocation;
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
import javastraw.tools.UNIXTools;

import java.awt.*;
import java.io.File;
import java.util.List;
import java.util.*;


public class Sift {
    private final int MAX_DIST = 10000000;

    private static final int ZSCORE_CUTOFF = 2;
    private static final NormalizationType scale = NormalizationHandler.SCALE;
    private final int MIN_DIST = 10000;
    private final int window = 5;

    public Sift(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            Main.printGeneralUsageAndExit(5);
        }
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1],
                false, false, false);
        File outFolder = UNIXTools.makeDir(new File(args[2]));
        Feature2DList refinedLoops = siftThroughCalls(ds);
        refinedLoops.exportFeatureList(new File(outFolder, "sift.bedpe"), false, Feature2DList.ListFormat.NA);
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
                    if (zscore > ZSCORE_CUTOFF) {
                        records.add(cr);
                    }
                }
            }
        }

        return records;
    }

    private static Set<SimpleLocation> getExtremeLocations(MatrixZoomData zd, int maxBin, int minBin, boolean useNone) {

        int maxCompressedBin = logp1i(maxBin) + 1;
        int minCompressedBin = logp1i(minBin);

        ZScores zScores = getZscores(zd, maxCompressedBin, useNone);

        Set<SimpleLocation> records = new HashSet<>();

        Iterator<ContactRecord> it;
        if (useNone) {
            it = zd.getDirectIterator();
        } else {
            it = zd.getNormalizedIterator(scale);
        }

        while (it.hasNext()) {
            ContactRecord cr = it.next();
            if (cr.getCounts() > 1) {
                int dist = logp1i(ExpectedUtils.getDist(cr));
                if (dist > minCompressedBin && dist < maxCompressedBin) {
                    float zscore = zScores.getZscore(dist, logp1(cr.getCounts()));
                    if (zscore > ZSCORE_CUTOFF) {
                        records.add(new SimpleLocation(cr));
                    }
                }
            }
        }

        return records;
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
     * 8 - coalesce pixels
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

                for (int lowres : new int[]{1000, 2000, 5000}) {
                    MatrixZoomData zd1 = matrix.getZoomData(new HiCZoom(lowres));
                    System.out.println("Start LowRes pass (" + lowres + ")");
                    Set<SimpleLocation> enrichedRegions = getExtremeLocations(zd1, MAX_DIST / lowres, MIN_DIST / lowres,
                            false);
                    filterOutByOverlap(initialPoints, enrichedRegions, lowres / hires);
                    enrichedRegions.clear();
                    System.out.println("LowRes pass done (" + lowres + ")");
                }
                matrix.clearCache();

                SiftUtils.coalescePixelsToCentroid(initialPoints, hires, 13000);

                output.addByKey(Feature2DList.getKey(chrom, chrom), convertToFeature2Ds(initialPoints,
                        chrom, chrom, hires));
            }
        }

        return output;
    }

    private static ZScores getZscores(MatrixZoomData zd, int length, boolean useNone) {
        WelfordStats stats = new WelfordStats(length);

        Iterator<ContactRecord> it;
        if (useNone) {
            it = zd.getDirectIterator();
        } else {
            it = zd.getNormalizedIterator(scale);
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
}