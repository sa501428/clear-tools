package cli.clt.flat.file;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.*;

public class LoopDiffFlatFileMaker {
    public static String usage = "create-diff-flat-file <genomeID> <output.bedpe> " +
            "<_positive_suffix.bedpe> <_negative_suffix.bedpe> <stem1,stem2,...,stemN>\n" +
            "creates a flat file, with a 1 if loop is present, -1 if not present," +
            "0 if indeterminate\n";

    public static void run(String[] args, String command, CommandLineParser parser) {
        if (args.length != 6) {
            Main.printGeneralUsageAndExit(58, usage);
        }

        String genomeID = args[1];
        String outPath = args[2];

        String posSuffix = args[3];
        String negSuffix = args[4];
        String[] stems = args[5].split(",");

        createFlatFile(genomeID, outPath, posSuffix, negSuffix, stems);

        System.out.println("flat file complete");
    }

    private static void createFlatFile(String genomeID, String outPath, String posSuffix, String negSuffix, String[] stems) {
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);

        Map<String, Map<String, Integer>> codeToStemToStatus = new HashMap<>();
        Map<String, Set<Feature2D>> chromKeyToLoops = new HashMap<>();

        for (String stem : stems) {
            Feature2DList loops = Feature2DParser.loadFeatures(stem + posSuffix, handler,
                    false, null, false);
            populateForEachLoopInList(codeToStemToStatus, chromKeyToLoops, stem, loops, 1);

            loops = Feature2DParser.loadFeatures(stem + negSuffix, handler,
                    false, null, false);
            populateForEachLoopInList(codeToStemToStatus, chromKeyToLoops, stem, loops, -1);
        }
        fillInIndeterminates(codeToStemToStatus, chromKeyToLoops, stems);

        Feature2DList flat = makeFromMapping(chromKeyToLoops, codeToStemToStatus);
        flat.exportFeatureList(new File(outPath), false, Feature2DList.ListFormat.NA);
    }

    private static Feature2DList makeFromMapping(Map<String, Set<Feature2D>> chromKeyToLoops, Map<String, Map<String, Integer>> codeToStemToStatus) {
        Feature2DList flat = new Feature2DList();

        for (String key : chromKeyToLoops.keySet()) {
            Set<Feature2D> loops = chromKeyToLoops.get(key);
            for (Feature2D loop : loops) {
                String code = getLoopCode(loop);
                Map<String, Integer> stemToStatus = codeToStemToStatus.get(code);
                setAttributes(loop, stemToStatus);
            }
            flat.addByKey(key, new LinkedList<>(loops));
        }
        return flat;
    }

    private static void setAttributes(Feature2D loop, Map<String, Integer> stemToStatus) {
        for (String stem : stemToStatus.keySet()) {
            loop.addIntAttribute(stem, stemToStatus.get(stem));
        }
    }

    private static void fillInIndeterminates(Map<String, Map<String, Integer>> codeToStemToStatus,
                                             Map<String, Set<Feature2D>> chromKeyToLoops, String[] stems) {
        for (String chromKey : chromKeyToLoops.keySet()) {
            Set<Feature2D> loops = chromKeyToLoops.get(chromKey);
            for (Feature2D loop : loops) {
                String code = getLoopCode(loop);
                if (!codeToStemToStatus.containsKey(code)) {
                    codeToStemToStatus.put(code, new HashMap<>());
                }
                for (String stem : stems) {
                    if (!codeToStemToStatus.get(code).containsKey(stem)) {
                        codeToStemToStatus.get(code).put(stem, 0);
                    }
                }
            }
        }
    }

    private static void populateForEachLoopInList(Map<String, Map<String, Integer>> codeToStemToStatus,
                                                  Map<String, Set<Feature2D>> chromKeyToLoops, String stem,
                                                  Feature2DList loops, int val) {
        loops.processLists((s, list) -> {
            for (Feature2D loop : list) {
                String code = getLoopCode(loop);
                if (!codeToStemToStatus.containsKey(code)) {
                    codeToStemToStatus.put(code, new HashMap<>());
                }
                codeToStemToStatus.get(code).put(stem, val);
            }
            if (!chromKeyToLoops.containsKey(s)) {
                chromKeyToLoops.put(s, new HashSet<>());
            }
            chromKeyToLoops.get(s).addAll(list);
        });
    }

    private static String getLoopCode(Feature2D loop) {
        return loop.getChr1() + "_" + loop.getStart1() + "_" + loop.getEnd1() + "_" +
                loop.getChr2() + "_" + loop.getStart2() + "_" + loop.getEnd2();
    }
}
