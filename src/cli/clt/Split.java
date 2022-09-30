package cli.clt;

import cli.Main;
import cli.utils.clean.LoopTools;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class Split {
    public static String usage = "split[-clean] <genomeID> <number> <input.bedpe> <output_>";

    public static void run(String[] args, String command) {
        if (args.length < 5) {
            Main.printGeneralUsageAndExit(5);
        }
        // split[-clean] <genomeID> <number> <input.bedpe> <output>
        String genomeID = args[1];
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        boolean noAttributes = command.contains("clean");

        int numberOfSplits = Integer.parseInt(args[2]);

        String path = args[3];
        Feature2DList loopList = LoopTools.loadFilteredBedpe(path, handler, !noAttributes);

        String outStem = args[4];

        Feature2DList[] outputs = split(loopList, numberOfSplits);

        for (int i = 0; i < outputs.length; i++) {
            outputs[i].exportFeatureList(new File(outStem + "split" + i + ".bedpe"),
                    false, Feature2DList.ListFormat.NA);
        }

        System.out.println("split complete");
    }

    private static Feature2DList[] split(Feature2DList loopList, int numberOfSplits) {
        if (numberOfSplits < 2) {
            return splitListByChromosomes(loopList);
        } else {
            return splitRandomly(loopList, numberOfSplits);
        }
    }

    private static Feature2DList[] splitListByChromosomes(Feature2DList loopList) {
        List<Feature2DList> finalLists = new ArrayList<>();
        loopList.processLists((s, list) -> {
            Feature2DList newList = new Feature2DList();
            newList.addByKey(s, list);
            synchronized (finalLists) {
                finalLists.add(newList);
            }
        });
        return finalLists.toArray(new Feature2DList[0]);
    }

    private static Feature2DList[] splitRandomly(Feature2DList loopList, int numberOfSplits) {
        Feature2DList[] lists = new Feature2DList[numberOfSplits];
        for (int z = 0; z < lists.length; z++) {
            lists[z] = new Feature2DList();
        }
        loopList.processLists((s, list) -> {
            List<List<Feature2D>> breakups = new ArrayList<>();
            for (int z = 0; z < numberOfSplits; z++) {
                breakups.add(new ArrayList<>());
            }

            int counter = 0;
            for (Feature2D feature : list) {
                breakups.get(counter).add(feature);
                counter = (counter + 1) % numberOfSplits;
            }

            synchronized (lists) {
                for (int z = 0; z < numberOfSplits; z++) {
                    lists[z].addByKey(s, breakups.get(z));
                }
            }
        });
        return lists;
    }
}
