package cli.clt;

import cli.Main;
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

public class SubtractSharedAnchors {

    private static final int resolution = 100;

    // overlap can be adjusted; exact means exact indices; default will use any overlap
    // clean means don't save attributes
    public static String usage = "subtract-shared-anchors[-clean] [-w window] <genomeID> " +
            "<fileA.bedpe> <fileB1.bedpe> <fileB2.bedpe> ...\n" +
            "\t\treturn only loops that have no shared anchors with anything that loops in files B1, B2, etc.\n";

    public static void run(String[] args, String command, CommandLineParser parser) {

        int window = parser.getWindowSizeOption(0);
        if (window > 0 && Main.printVerboseComments) {
            System.out.println("Anchors will be expanded by " + window);
        }

        if (args.length < 5) {
            Main.printGeneralUsageAndExit(15, usage);
        }

        boolean noAttributes = command.contains("clean");
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        String inputName = args[2];
        Feature2DList featuresA = Feature2DParser.loadFeatures(inputName, handler, !noAttributes, null, false);

        Map<Integer, BitSet> upstreamAnchors = new HashMap<>();
        Map<Integer, BitSet> downstreamAnchors = new HashMap<>();

        populateFilteringAnchors(upstreamAnchors, downstreamAnchors, handler, args, window, 3);

        Feature2DList output = filterFeatures(featuresA, handler, upstreamAnchors, downstreamAnchors);
        output.exportFeatureList(new File(getOutputName(inputName)), false, Feature2DList.ListFormat.NA);
        System.out.println("anchor filtering complete");

    }

    private static String getOutputName(String inputName) {
        return inputName.replace(".bedpe", "") + "_after_anchor_filter.bedpe";
    }

    private static Feature2DList filterFeatures(Feature2DList featuresA, ChromosomeHandler handler,
                                                Map<Integer, BitSet> upstreamAnchors,
                                                Map<Integer, BitSet> downstreamAnchors) {
        final Feature2DList result = new Feature2DList();
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < chromosomes.length) {
                Chromosome chromosome = chromosomes[i];
                String key = Feature2DList.getKey(chromosome, chromosome);
                List<Feature2D> features = intersect(featuresA.get(key),
                        upstreamAnchors.get(chromosome.getIndex()),
                        downstreamAnchors.get(chromosome.getIndex()));
                synchronized (result) {
                    result.addByKey(key, features);
                }
                i = index.getAndIncrement();
            }
        });
        return result;
    }

    private static List<Feature2D> intersect(List<Feature2D> features, BitSet upStream, BitSet downStream) {
        final Set<Feature2D> toSaveFinal = new HashSet<>();
        for (Feature2D feature : features) {
            boolean overlapsAnchor = hasAnchor(upStream, feature.getStart1(), feature.getEnd1()) ||
                    hasAnchor(downStream, feature.getStart2(), feature.getEnd2());
            if (!overlapsAnchor) {
                toSaveFinal.add(feature);
            }
        }
        return new ArrayList<>(toSaveFinal);
    }

    private static boolean hasAnchor(BitSet stream, long start1, long end1) {
        for (int x = (int) (start1 / resolution); x < end1 / resolution; x++) {
            if (stream.get(x)) return true;
        }
        return false;
    }

    private static void populateFilteringAnchors(Map<Integer, BitSet> upstreamAnchors,
                                                 Map<Integer, BitSet> downstreamAnchors,
                                                 ChromosomeHandler handler, String[] args, int window,
                                                 int startIndex) {
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            int length = (int) ((chrom.getLength() / resolution) + 1);
            upstreamAnchors.put(chrom.getIndex(), new BitSet(length));
            downstreamAnchors.put(chrom.getIndex(), new BitSet(length));
        }

        for (int k = startIndex; k < args.length; k++) {
            Feature2DList featuresB = Feature2DParser.loadFeatures(args[k], handler, false, null, false);
            for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
                BitSet[] bitSets = getUpStreamDownStreamAnchors(featuresB.get(chrom.getIndex(),
                        chrom.getIndex()), chrom, window);
                upstreamAnchors.get(chrom.getIndex()).or(bitSets[0]);
                downstreamAnchors.get(chrom.getIndex()).or(bitSets[1]);
            }
        }
    }

    private static BitSet[] getUpStreamDownStreamAnchors(List<Feature2D> features, Chromosome chrom, int window) {
        int length = (int) ((chrom.getLength() / resolution) + 1);
        BitSet bitSetUpStream = new BitSet(length);
        BitSet bitSetDownStream = new BitSet(length);
        for (Feature2D feature : features) {
            populate(bitSetUpStream, feature.getStart1() - window, feature.getEnd1() + window);
            populate(bitSetDownStream, feature.getStart2() - window, feature.getEnd2() + window);
        }
        return new BitSet[]{bitSetUpStream, bitSetDownStream};
    }

    private static void populate(BitSet anchorStream, long start1, long end1) {
        int realStart = (int) Math.max(start1 / resolution, 0);
        int realEnd = (int) Math.max(realStart + 1, end1 / resolution);
        anchorStream.set(realStart, realEnd);
    }
}
