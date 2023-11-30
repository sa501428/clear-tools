package cli.clt.flat.file;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

public class GetDiffsFromFlatFile {
    public static String usage = "get-diffs-from-flat-file <genomeID> <flat.file.bedpe> <stem1> <stem2> " +
            "<outfolder>\n" +
            "creates shared and differential loop lists from the flat file\n";

    public static void run(String[] args, String command, CommandLineParser parser) {
        if (args.length != 6) {
            Main.printGeneralUsageAndExit(58, usage);
        }

        String genomeID = args[1];
        String flatFile = args[2];

        String stem1 = args[3];
        String stem2 = args[4];

        String outFolder = args[5];
        File folder = UNIXTools.makeDir(new File(outFolder));

        getFromFlatFile(genomeID, flatFile, stem1, stem2, folder);

        System.out.println("extraction complete");
    }

    private static void getFromFlatFile(String genomeID, String flatFile, String stem1, String stem2, File folder) {
        Feature2DList flat = Feature2DParser.loadFeatures(flatFile, ChromosomeTools.loadChromosomes(genomeID),
                true, null, false);

        Feature2DList onlyStem1 = new Feature2DList();
        Feature2DList onlyStem2 = new Feature2DList();
        Feature2DList bothStems = new Feature2DList();
        Feature2DList neitherStems = new Feature2DList();

        flat.processLists((s, list) -> {
            List<Feature2D> only1 = new LinkedList<>();
            List<Feature2D> only2 = new LinkedList<>();
            List<Feature2D> both = new LinkedList<>();
            List<Feature2D> neither = new LinkedList<>();
            System.out.println("processing " + s);

            for (Feature2D loop : list) {
                int status1 = Integer.parseInt(loop.getAttribute(stem1));
                int status2 = Integer.parseInt(loop.getAttribute(stem2));

                if (status1 == 1 && status2 == 1) {
                    both.add(loop);
                } else if (status1 == 1 && status2 == -1) {
                    only1.add(loop);
                } else if (status2 == 1 && status1 == -1) {
                    only2.add(loop);
                } else if (status1 == -1 && status2 == -1) {
                    neither.add(loop);
                }
            }
            onlyStem1.addByKey(s, only1);
            onlyStem2.addByKey(s, only2);
            bothStems.addByKey(s, both);
            neitherStems.addByKey(s, neither);
        });

        onlyStem1.exportFeatureList(new File(folder, "only." + stem1 + ".bedpe"), false, Feature2DList.ListFormat.NA);
        onlyStem2.exportFeatureList(new File(folder, "only." + stem2 + ".bedpe"), false, Feature2DList.ListFormat.NA);
        bothStems.exportFeatureList(new File(folder, "both." + stem1 + "." + stem2 + ".bedpe"), false, Feature2DList.ListFormat.NA);
        neitherStems.exportFeatureList(new File(folder, "neither." + stem1 + "." + stem2 + ".bedpe"), false, Feature2DList.ListFormat.NA);
    }
}
