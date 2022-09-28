package cli.clt;

import cli.Main;
import cli.utils.clean.LoopTools;
import cli.utils.general.QuickGrouping;
import cli.utils.sift.SimpleLocation;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class IntersectBedpe {

    public static String usage = "intersect[-subtract][-clean][-exact] [-w window] <genomeID> <fileA.bedpe> <fileB.bedpe> <output.bedpe>";

    public static void run(String[] args, String command, CommandLineParser parser) {

        // subtract will keep things in A that don't overlap with B
        // otherwise retain things in A that have overlap with B

        boolean doSubtraction = checkForSubtraction(command);
        boolean useExactMatch = checkForExact(command);

        int window = parser.getWindowSizeOption(0);
        if (window > 0 && Main.printVerboseComments) {
            System.out.println("Features will be expanded by " + window);
        }

        // overlap can be adjusted
        // exact means exact
        // default will use any overlap

        // clean means don't save attributes

        if (args.length != 5) {
            Main.printGeneralUsageAndExit(15);
        }
        // intersect <genomeID> <fileA.bedpe> <fileB.bedpe> <output.bedpe>
        String genomeID = args[1];
        String inputA = args[2];
        String inputB = args[3];
        String outFile = args[4];

        boolean noAttributes = command.contains("clean");

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList featuresA = LoopTools.loadFilteredBedpe(inputA, handler, !noAttributes);
        Feature2DList featuresB = LoopTools.loadFilteredBedpe(inputB, handler, !noAttributes);


        Feature2DList output = coalesceFeatures(featuresA, featuresB, handler,
                doSubtraction, useExactMatch, window);
        output.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
        System.out.println("fusion complete");

    }

    private static Feature2DList coalesceFeatures(Feature2DList featuresA, Feature2DList featuresB,
                                                  ChromosomeHandler handler, boolean doSubtraction,
                                                  boolean useExactMatch, int window) {
        Feature2DList result = new Feature2DList();
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                String key = Feature2DList.getKey(chromosomes[i], chromosomes[j]);
                List<Feature2D> features = intersect(featuresA.get(key), featuresB.get(key),
                        doSubtraction, useExactMatch, window);
                result.addByKey(key, features);
            }
        }
        return result;
    }

    private static List<Feature2D> intersect(List<Feature2D> listA, List<Feature2D> listB,
                                             boolean doSubtraction, boolean useExactMatch, int window) {
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

                processFeatures(coalesced, featuresA, featuresB, doSubtraction, useExactMatch, window);

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

    private static void processFeatures(Set<Feature2D> coalesced, List<Feature2D> featuresA, List<Feature2D> featuresB, boolean doSubtraction, boolean useExactMatch, int window) {

        // todo
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
}
