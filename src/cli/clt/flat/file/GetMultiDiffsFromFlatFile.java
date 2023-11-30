package cli.clt.flat.file;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.*;

import static cli.clt.sieve.RetainOverlap.getLoopCode;

public class GetMultiDiffsFromFlatFile {
    public static String usage = "get-multi-diffs-from-flat-file <genomeID> <flat.file.bedpe> <outfolder> " +
            "<stem1,stem2,...> <stemA,stemB,...> ...\n" +
            "creates differential loop lists from the flat file\n";

    public static void run(String[] args, String command, CommandLineParser parser) {
        if (args.length < 6) {
            Main.printGeneralUsageAndExit(58, usage);
        }

        String genomeID = args[1];
        String flatFile = args[2];
        String outFolder = args[3];

        String[][] stems = new String[args.length - 4][];
        for (int i = 0; i < stems.length; i++) {
            stems[i] = args[i + 4].split(",");
        }

        File folder = UNIXTools.makeDir(new File(outFolder));

        getFromFlatFile(genomeID, flatFile, stems, folder);

        System.out.println("extraction complete");
    }

    private static void getFromFlatFile(String genomeID, String flatFile, String[][] stems, File folder) {
        Feature2DList flat = Feature2DParser.loadFeatures(flatFile, ChromosomeTools.loadChromosomes(genomeID),
                true, null, false);

        Feature2DList[] onlyStemN = new Feature2DList[stems.length];
        Feature2DList[] onlyNotStemN = new Feature2DList[stems.length];

        for (int i = 0; i < stems.length; i++) {
            onlyStemN[i] = new Feature2DList();
            onlyNotStemN[i] = new Feature2DList();
        }

        flat.processLists((s, list) -> {
            Map<String, Feature2D> keyToLoop = makeMapping(list);

            List<List<Feature2D>> present = new ArrayList<>(stems.length);
            List<List<Feature2D>> absent = new ArrayList<>(stems.length);
            for (int i = 0; i < stems.length; i++) {
                present.add(new LinkedList<>());
                absent.add(new LinkedList<>());
            }

            System.out.println("processing " + s);

            Map<String, Map<Integer, Float>> keyToStemCount = makeEmptyMap(keyToLoop.keySet());
            for (Feature2D loop : list) {
                for (int i = 0; i < stems.length; i++) {
                    float counter = 0;
                    for (String stem : stems[i]) {
                        counter += Integer.parseInt(loop.getAttribute(stem));
                    }
                    keyToStemCount.get(getLoopCode(loop)).put(i, counter / stems[i].length);
                }
            }

            for (String loopKey : keyToStemCount.keySet()) {
                int whichList = getUniqueLoop(keyToStemCount.get(loopKey), true, stems.length);
                if (whichList > -1) {
                    present.get(whichList).add(keyToLoop.get(loopKey));
                }
                whichList = getUniqueLoop(keyToStemCount.get(loopKey), false, stems.length);
                if (whichList > -1) {
                    absent.get(whichList).add(keyToLoop.get(loopKey));
                }
            }

            for (int i = 0; i < stems.length; i++) {
                if (present.get(i).size() > 0) {
                    onlyStemN[i].addByKey(s, present.get(i));
                }
                if (absent.get(i).size() > 0) {
                    onlyNotStemN[i].addByKey(s, absent.get(i));
                }
            }
        });

        for (int i = 0; i < stems.length; i++) {
            onlyStemN[i].exportFeatureList(new File(folder, "only." + String.join(".", stems[i]) + ".bedpe"), false, Feature2DList.ListFormat.NA);
            onlyNotStemN[i].exportFeatureList(new File(folder, "not." + String.join(".", stems[i]) + ".bedpe"), false, Feature2DList.ListFormat.NA);
        }
    }

    private static Map<String, Map<Integer, Float>> makeEmptyMap(Set<String> keySet) {
        Map<String, Map<Integer, Float>> map = new HashMap<>();
        for (String key : keySet) {
            map.put(key, new HashMap<>());
        }
        return map;
    }

    private static Map<String, Feature2D> makeMapping(List<Feature2D> list) {
        Map<String, Feature2D> mapping = new HashMap<>();
        for (Feature2D loop : list) {
            mapping.put(getLoopCode(loop), loop);
        }
        return mapping;
    }

    private static int getUniqueLoop(Map<Integer, Float> integerFloatMap, boolean isPresent, int n) {
        if (isPresent) {
            int numEqualOne = 0;
            int numNotEqualOne = 0;
            for (int i : integerFloatMap.keySet()) {
                if (integerFloatMap.get(i) > 0.99) {
                    numEqualOne++;
                } else if (integerFloatMap.get(i) < 0) {
                    numNotEqualOne++;
                }
            }
            if (numEqualOne == 1 && numNotEqualOne > n / 2) {
                for (int i : integerFloatMap.keySet()) {
                    if (integerFloatMap.get(i) > 0.99) {
                        return i;
                    }
                }
            }
        } else {
            int numEqualNegOne = 0;
            int numNotEqualNegOne = 0;
            for (int i : integerFloatMap.keySet()) {
                if (integerFloatMap.get(i) < -0.99) {
                    numEqualNegOne++;
                } else if (integerFloatMap.get(i) > 0) {
                    numNotEqualNegOne++;
                }
            }
            if (numEqualNegOne == 1 && numNotEqualNegOne > n / 2) {
                for (int i : integerFloatMap.keySet()) {
                    if (integerFloatMap.get(i) < -0.99) {
                        return i;
                    }
                }
            }
        }
        return -1;
    }
}
