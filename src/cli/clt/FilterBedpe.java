package cli.clt;

import cli.Main;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class FilterBedpe {

    // overlap can be adjusted; exact means exact indices; default will use any overlap
    //
    public static String usage = "filter-contain[-clean] <genomeID> <input.bedpe> <output.bedpe>\n" +
            "\t\tcontain checks if localX and localY are within the bounds of the feature\n" +
            "\t\tclean means don't save attributes";

    public static void run(String[] args, String command) {

        boolean doCheckContain = checkForContainment(command);

        if (args.length != 4) {
            Main.printGeneralUsageAndExit(15, usage);
        }
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        boolean noAttributes = command.contains("clean");
        Feature2DList features = Feature2DParser.loadFeatures(args[2], handler, !noAttributes, null, false);

        Feature2DList output = filter(features, doCheckContain);
        output.exportFeatureList(new File(args[3]), false, Feature2DList.ListFormat.NA);
        System.out.println("filtering complete");

    }

    private static Feature2DList filter(Feature2DList features, boolean checkContain) {
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
