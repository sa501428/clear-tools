package cli.clt;

import cli.Main;
import cli.utils.expected.ExpectedModel;
import cli.utils.expected.LogExpectedPolynomial;
import cli.utils.sift.ContactRecordBox;
import cli.utils.sift.ExtremePixels;
import cli.utils.sift.FeatureUtils;
import cli.utils.sift.SimpleLocation;
import cli.utils.sift.collapse.MultiResCentroidCollapser;
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
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;


public class Sift {

    // 410 defaults
    // 411 require only 2 resolutions (instead of 3)
    // 412 norm 1->0.9, EVE 1.75->1.5, EVN 1.2 -> 1.1
    // 413 expand bounds in second pass
    // 414 MPC .005 -> .001, MPC 0.4 -> 0.5
    // 415 MPC restored, norm 1, EVE 1.5, EVN 1.5, require 2 resolutions
    // 416 EVN 1.25
    // 417 EVN 1.5,

    /* hardcoded variables for sifting through pixels */
    public static final float MIN_PC = 0.005f;
    public static final float MAX_PC = 0.4f;
    public static final int MIN_RADIUS_0 = 1000;
    public static final float MIN_NORM = 1;
    public static final float ENRICMENT_VS_EXPECTED = 1.5f;
    public static final float ENRICHMENT_VS_NEIGHBORS = 1.5f;
    private static final int NUM_RES_TO_PASS = 2;


    private static final int[] resolutions = new int[]{100, 200, 500, 1000, 2000, 5000, 10000}; //  10000
    private NormalizationType norm = NormalizationHandler.NONE;
    private static final int MAX_DIST = 10000000;
    private static final int MIN_DIST = 10000;

    public Sift(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            Main.printGeneralUsageAndExit(5);
        }

        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);

        String normString = parser.getNormalizationStringOption();
        if (normString != null && normString.length() > 1) {
            norm = ds.getNormalizationHandler().getNormTypeFromString(normString);
        }
        System.out.println("Using norm: " + norm.getLabel());

        Feature2DList refinedLoops = siftThroughCalls(ds, args[2]);
        refinedLoops.exportFeatureList(new File(args[2] + ".sift.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("sift complete");
    }

    private Feature2DList siftThroughCalls(Dataset ds, String outname) {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Feature2DList output = new Feature2DList();

        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();

        AtomicInteger cIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chromosome = chromosomes[currIndex];
                List<Feature2D> sharpLoops = findSiftedFeatures(chromosome, ds, outname);
                if (sharpLoops.size() > 0) {
                    synchronized (output) {
                        output.addByKey(Feature2DList.getKey(chromosome, chromosome), sharpLoops);
                    }
                }
                currIndex = cIndex.getAndIncrement();
            }
        });

        return output;
    }

    private List<Feature2D> findSiftedFeatures(Chromosome chromosome, Dataset ds, String outname) {
        Matrix matrix = ds.getMatrix(chromosome, chromosome);
        if (matrix == null) return new ArrayList<>();

        final Map<Integer, Set<ContactRecord>> pixelsForResolutions = new HashMap<>();

        AtomicInteger rIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currResIndex = rIndex.getAndIncrement();
            while (currResIndex < resolutions.length) {
                int lowRes = resolutions[currResIndex];

                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(lowRes));
                if (zd != null) {
                    ExpectedModel poly = new LogExpectedPolynomial(zd, norm, MAX_DIST / lowRes);
                    Set<ContactRecord> points = ExtremePixels.getExtremePixelsForResolution(ds, zd,
                            chromosome, lowRes, norm, MAX_DIST / lowRes, MIN_DIST / lowRes, poly);
                    matrix.clearCacheForZoom(new HiCZoom(lowRes));

                    synchronized (pixelsForResolutions) {
                        pixelsForResolutions.put(lowRes, points);
                        System.out.println(lowRes + " found (" + points.size() + ")");
                    }
                }
                currResIndex = rIndex.getAndIncrement();
            }
        });

        matrix.clearCache();

        Set<ContactRecordBox> allRecords = collapse(pixelsForResolutions);
        pixelsForResolutions.clear();

        Set<ContactRecordBox> finalBoxes = MultiResCentroidCollapser.coalesce(allRecords, NUM_RES_TO_PASS);
        allRecords.clear();

        return FeatureUtils.convertToFeature2Ds(finalBoxes, chromosome);
    }

    private Set<ContactRecordBox> collapse(Map<Integer, Set<ContactRecord>> recordMap) {
        HashSet<ContactRecordBox> records = new HashSet<>();
        for (Integer res : recordMap.keySet()) {
            for (ContactRecord record : recordMap.get(res)) {
                records.add(new ContactRecordBox(record, res));
            }
        }
        return records;
    }

    private Feature2DList convert(Set<SimpleLocation> locations, Chromosome chromosome, int res) {
        Feature2DList list = new Feature2DList();
        List<Feature2D> features = new ArrayList<>();
        for (SimpleLocation location : locations) {
            features.add(location.toRegion(res).toFeature2D(chromosome));
        }
        list.addByKey(Feature2DList.getKey(chromosome, chromosome), features);
        return list;
    }
}
