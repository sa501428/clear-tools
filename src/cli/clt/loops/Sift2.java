package cli.clt.loops;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.anchors.GlobalAnchorCounter;
import cli.utils.data.InterLogZscoreExpectedModel;
import cli.utils.seer.SeerUtils;
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

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;


public class Sift2 {

    /*

    TODO
    of the (n(n - 1)/2) interactions, > 70% should be valid hits
    call the clique
    compare with loci cliques if shifted by 1 or 10 positions upstream / down stream, should still be enriched
    inter-chromosomal or regular APA on the clique
    test with hct 116 rad21 treated?
     */

    public static String usage = "inter-sift [--res int] <file.hic> <outfile.stem>";
    private static int resolution = 50000;
    private NormalizationType norm = NormalizationHandler.VC;

    public Sift2(String[] args, CommandLineParser parser) {
        if (args.length != 3) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        Dataset ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, false);

        resolution = parser.getResolutionOption(resolution);

        String normString = parser.getNormalizationStringOption();
        if (normString != null && normString.length() > 1) {
            norm = ds.getNormalizationHandler().getNormTypeFromString(normString);
        }
        System.out.println("Using norm: " + norm.getLabel());

        siftThroughCalls(ds, args[2]);

        //interClique.exportFeatureList(new File(args[2] + ".cliques.bedpe"), false, Feature2DList.ListFormat.NA);
        System.out.println("inter sift complete");
    }

    private void siftThroughCalls(Dataset ds, String outStem) {
        ChromosomeHandler handler = ds.getChromosomeHandler();

        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        Map<Chromosome, Map<Chromosome, List<ContactRecord>>> globalPeakMap = initialize(chromosomes);

        AtomicInteger cIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int currIndex = cIndex.getAndIncrement();
            while (currIndex < chromosomes.length) {
                Chromosome chromosome1 = chromosomes[currIndex];
                for (int j = currIndex + 1; j < chromosomes.length; j++) {
                    Chromosome chromosome2 = chromosomes[j];
                    List<ContactRecord> interLoops = findSiftedFeatures(chromosome1, chromosome2, ds);

                    if (interLoops.size() > 0) {
                        synchronized (globalPeakMap) {
                            populate(globalPeakMap, chromosome1, chromosome2, interLoops);
                        }
                    }
                }
                currIndex = cIndex.getAndIncrement();
            }
        });

        // check for shared / overlapping rows, retain loci > 3 hits
        Map<Chromosome, int[]> globalAnchors = GlobalAnchorCounter.createFromPeaks(chromosomes, globalPeakMap, resolution);

        try {
            SeerUtils.exportRowSumsToBedgraph(globalAnchors, outStem, resolution);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private Map<Chromosome, Map<Chromosome, List<ContactRecord>>> initialize(Chromosome[] chromosomes) {
        Map<Chromosome, Map<Chromosome, List<ContactRecord>>> map = new HashMap<>();
        for (int i = 0; i < chromosomes.length; i++) {
            map.put(chromosomes[i], new HashMap<>());
            for (int j = i + 1; j < chromosomes.length; j++) {
                map.get(chromosomes[i]).put(chromosomes[j], new ArrayList<>());
            }
        }
        return map;
    }

    private void populate(Map<Chromosome, Map<Chromosome, List<ContactRecord>>> output,
                          Chromosome chromosome1,
                          Chromosome chromosome2, List<ContactRecord> interLoops) {
        output.get(chromosome1).put(chromosome2, interLoops);
    }

    private List<ContactRecord> findSiftedFeatures(Chromosome chromosomeI, Chromosome chromosomeJ, Dataset ds) {
        Matrix matrix = ds.getMatrix(chromosomeI, chromosomeJ);
        if (matrix == null) return new ArrayList<>();


        MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
        if (zd != null) {
            InterLogZscoreExpectedModel expected = new InterLogZscoreExpectedModel(zd, norm);
            // identify all z score > 2 loci interchromosomally after log
            List<ContactRecord> records = expected.getAllContactsZscore(zd, norm, 2);
            Map<Integer, Map<Integer, Float>> recordsAsMap = convert(records);
            // non-max suppression / simple check for all neighbors are less than the pixel)
            List<ContactRecord> recordsNMS = nonMaxSuppression(records, recordsAsMap);
            records.clear();
            recordsAsMap.clear();
            return recordsNMS;
        }
        return new ArrayList<>();
    }

    private List<ContactRecord> nonMaxSuppression(List<ContactRecord> records,
                                                  Map<Integer, Map<Integer, Float>> recordsAsMap) {
        List<ContactRecord> survivedNMS = new LinkedList<>();
        for (ContactRecord record : records) {
            boolean isMax = true;
            int i = record.getBinX();
            for (int di = -1; di < 2; di++) {
                int j = record.getBinY();
                for (int dj = -1; dj < 2; dj++) {
                    if (di == 0 && dj == 0) continue;
                    if (recordsAsMap.containsKey(i + di) &&
                            recordsAsMap.get(i + di).containsKey(j + dj)) {
                        if (recordsAsMap.get(i + di).get(j + dj) > record.getCounts()) {
                            isMax = false;
                            break;
                        }
                    }
                }
            }
            if (isMax) {
                survivedNMS.add(record);
            }
        }
        return survivedNMS;
    }

    private Map<Integer, Map<Integer, Float>> convert(List<ContactRecord> records) {
        Map<Integer, Map<Integer, Float>> topPixels = new HashMap<>();
        for (ContactRecord record : records) {
            if (!topPixels.containsKey(record.getBinX())) {
                topPixels.put(record.getBinX(), new HashMap<>());
            }
            topPixels.get(record.getBinX()).put(record.getBinY(), record.getCounts());
        }
        return topPixels;
    }
}
