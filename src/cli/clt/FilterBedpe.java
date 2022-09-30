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

public class FilterBedpe {

    // overlap can be adjusted; exact means exact indices; default will use any overlap
    // clean means don't save attributes
    public static String usage = "filter[-contain][-clean] <genomeID> <input.bedpe> <output.bedpe>";

    public static void run(String[] args, String command) {

        boolean doCheckContain = checkForContainment(command);

        if (args.length != 4) {
            Main.printGeneralUsageAndExit(15);
        }
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        Feature2DList features = LoopTools.loadFilteredBedpe(args[2], handler, true);

        Feature2DList output = filter(features, handler, doCheckContain);
        output.exportFeatureList(new File(args[3]), false, Feature2DList.ListFormat.NA);
        System.out.println("filtering complete");

    }

    private static Feature2DList filter(Feature2DList features, ChromosomeHandler handler,
                                        boolean checkContain) {
        features.filterLists((s, list) -> filterForContainment(list, checkContain));

        return features;
    }

    private static List<Feature2D> filterForContainment(List<Feature2D> features, boolean checkContain) {
        List<Feature2D> newList = new ArrayList<>();
        for (Feature2D feature : features) {
            if (checkContain) {
                if (properlyContainedLocalization(feature)) {
                    newList.add(feature);
                }
            }
        }
        return newList;
    }

    private static boolean properlyContainedLocalization(Feature2D feature) {
        return localizationContained(feature.getStart1(), feature.getEnd1(), feature.getAttribute("localX"))
                && localizationContained(feature.getStart2(), feature.getEnd2(), feature.getAttribute("localY"));
    }

    private static boolean localizationContained(long start, long end, String local) {
        try {
            long position = Long.parseLong(local);
            return start <= position && position <= end;
        } catch (Exception e) {
            return false;
        }
    }


    private static boolean checkForContainment(String command) {
        if (command.contains("contain")) {
            if (Main.printVerboseComments) System.out.println("Checking containment");
            return true;
        } else {
            if (Main.printVerboseComments) System.out.println("Not implemented!");
            System.exit(33);
            return false;
        }
    }
}
