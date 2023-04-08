package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.general.FusionTools;
import cli.utils.general.OverlapTools;
import cli.utils.general.QuickGrouping;
import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class IntersectBedpe {

    // overlap can be adjusted; exact means exact indices; default will use any overlap
    // clean means don't save attributes
    public static String usage = "intersect[-subtract][-clean][-exact][-bounding-box] [-w window] <genomeID> " +
            "<fileA.bedpe> <fileB.bedpe> <output.bedpe>\n" +
            "\t\tsubtract retains features in A with no overlap in B, default requires an overlap\n" +
            "\t\texact requires exact matches, default is any overlap";

    public static void run(String[] args, String command, CommandLineParser parser) {

        // subtract will keep things in A that don't overlap with B
        // otherwise retain things in A that have overlap with B
        boolean doSubtraction = checkForSubtraction(command);
        boolean useExactMatch = checkForExact(command);
        boolean doBoundingBox = checkForBounds(command);

        int window = parser.getWindowSizeOption(0);
        if (window > 0 && Main.printVerboseComments) {
            System.out.println("Features will be expanded by " + window);
        }

        if (args.length != 5) {
            Main.printGeneralUsageAndExit(15, usage);
        }
        // intersect <genomeID> <fileA.bedpe> <fileB.bedpe> <output.bedpe>
        boolean noAttributes = command.contains("clean");
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        Feature2DList featuresA = Feature2DParser.loadFeatures(args[2], handler, !noAttributes, null, false);
        Feature2DList featuresB = Feature2DParser.loadFeatures(args[3], handler, !noAttributes, null, false);

        Feature2DList output = coalesceFeatures(featuresA, featuresB, handler,
                doSubtraction, useExactMatch, window, doBoundingBox);
        output.exportFeatureList(new File(args[4]), false, Feature2DList.ListFormat.NA);
        System.out.println("fusion complete");

    }

    private static Feature2DList coalesceFeatures(Feature2DList featuresA, Feature2DList featuresB,
                                                  ChromosomeHandler handler, boolean doSubtraction,
                                                  boolean useExactMatch, int window, boolean doBoundingBox) {
        Feature2DList result = new Feature2DList();
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                String key = Feature2DList.getKey(chromosomes[i], chromosomes[j]);
                List<Feature2D> features = intersect(featuresA.get(key), featuresB.get(key),
                        doSubtraction, useExactMatch, window, doBoundingBox);
                result.addByKey(key, features);
            }
        }
        return result;
    }

    private static List<Feature2D> intersect(List<Feature2D> listA, List<Feature2D> listB,
                                             boolean doSubtraction, boolean useExactMatch, int window,
                                             boolean doBoundingBox) {
        Map<SimpleLocation, List<Feature2D>> mapA = QuickGrouping.groupNearbyRecords(listA, 1000000);
        Map<SimpleLocation, List<Feature2D>> mapB = QuickGrouping.groupNearbyRecordsWithOverlap(listB, 1000000);

        List<SimpleLocation> keys = new ArrayList<>(mapA.keySet());
        Set<Feature2D> allCoalesced = new HashSet<>();

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {

            Set<Feature2D> coalesced = new HashSet<>();
            int i = index.getAndIncrement();
            while (i < keys.size()) {

                SimpleLocation location = keys.get(i);
                List<Feature2D> featuresA = mapA.get(location);
                List<Feature2D> featuresB = new ArrayList<>();
                if (mapB.containsKey(location)) {
                    featuresB = mapB.get(location);
                }

                processFeatures(coalesced, featuresA, featuresB, doSubtraction, useExactMatch, window, doBoundingBox);

                i = index.getAndIncrement();
            }

            synchronized (allCoalesced) {
                allCoalesced.addAll(coalesced);
            }
        });

        mapA.clear();
        mapB.clear();

        return new ArrayList<>(allCoalesced);
    }

    private static void processFeatures(Set<Feature2D> coalesced, List<Feature2D> featuresA, List<Feature2D> featuresB,
                                        boolean doSubtraction, boolean useExactMatch, int window,
                                        boolean doBoundingBox) {
        if (useExactMatch) {
            Set<Feature2D> setA = new HashSet<>(featuresA);
            if (doSubtraction) {
                featuresB.forEach(setA::remove);
            } else {
                setA.retainAll(featuresB);
            }
            coalesced.addAll(setA);
            setA.clear();
            return;
        }

        for (Feature2D pixelA : featuresA) {
            List<Feature2D> pixelList = OverlapTools.getMatchesWithOverlap(pixelA, featuresB, window);

            if (doSubtraction) {
                if (pixelList.isEmpty()) coalesced.add(pixelA);
            } else if (pixelList.size() > 0) {
                if (doBoundingBox) {
                    pixelList.add(pixelA);
                    coalesced.add(FusionTools.getFeatureFromBounds(pixelList));
                } else {
                    coalesced.add(pixelA);
                }
            }
        }
    }

    private static boolean checkForSubtraction(String command) {
        if (command.contains("subtract")) {
            if (Main.printVerboseComments) System.out.println("Doing subtraction");
            return true;
        } else {
            if (Main.printVerboseComments) System.out.println("Doing intersection");
            return false;
        }
    }

    private static boolean checkForExact(String command) {
        if (command.contains("exact")) {
            if (Main.printVerboseComments) System.out.println("Will use exact matches");
            return true;
        } else {
            if (Main.printVerboseComments) System.out.println("Will use any amount of overlap");
            return false;
        }
    }

    private static boolean checkForBounds(String command) {
        if (command.contains("bound")) {
            if (Main.printVerboseComments) System.out.println("Will use bounding box of overlapping features");
            return true;
        } else {
            return false;
        }
    }
}
