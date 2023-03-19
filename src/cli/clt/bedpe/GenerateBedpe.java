package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.flags.Anchor;
import cli.utils.general.BedFileParser;
import cli.utils.loops.AnchorTools;
import cli.utils.loops.DomainTools;
import cli.utils.loops.OffsetTools;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;

public class GenerateBedpe {

    public static String usage = "generate-from-anchors <forward.bed> <reverse.bed> " +
            "<min_genome_dist> <max_genome_dist> <genomeID> <output.bedpe>\n" +
            "\t\tcreate potential loop locations using the anchors\n\n" +
            "\t\tgenerate-from-domains <genomeID> <domains.bedpe> <output_>\n\n" +
            "\t\tgenerate-from-offsets <genomeID> <loops.bedpe> <output.bedpe> <+-offset>";

    public static void run(String[] args, CommandLineParser parser, String command) {
        if (command.contains("anchors")) {
            buildFromAnchors(args, parser);
        } else if (command.contains("domains")) {
            buildFromDomains(args, parser);
        } else if (command.contains("offset")) {
            buildFromOffsets(args, parser);
        }
        System.out.println("generation complete");
    }

    public static void buildFromAnchors(String[] args, CommandLineParser parser) {

        if (args.length != 7) {
            Main.printGeneralUsageAndExit(4, usage);
        }
        String forwardMotifFile = args[1];
        String reverseMotifFile = args[2];
        long minDist = Long.parseLong(args[3]);
        long maxDist = Long.parseLong(args[4]);
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[5]);
        String outname = args[6];

        int resolution = parser.getResolutionOption(0);
        int percentile = parser.getPercentileOption(-1);

        GenomeWide1DList<Anchor> forwardAnchors = BedFileParser.loadFromBEDFile(handler, forwardMotifFile, percentile, true, true);
        GenomeWide1DList<Anchor> reverseAnchors = BedFileParser.loadFromBEDFile(handler, reverseMotifFile, percentile, true, true);
        Feature2DList output = AnchorTools.createLoops(handler, forwardAnchors, reverseAnchors, minDist, maxDist, resolution);
        output.exportFeatureList(new File(outname), false, Feature2DList.ListFormat.NA);
    }

    private static void buildFromDomains(String[] args, CommandLineParser parser) {
        // generate-from-domains <genomeID> <domains.bedpe> <output_>
        if (args.length != 4) {
            Main.printGeneralUsageAndExit(4, usage);
        }
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        Feature2DList domains = Feature2DParser.loadFeatures(args[2], handler, false, null, false);
        String outStem = args[3];
        Feature2DList[] output = DomainTools.createLoops(domains);
        output[0].exportFeatureList(new File(outStem + ".corners.bedpe"), false, Feature2DList.ListFormat.NA);
        output[1].exportFeatureList(new File(outStem + ".centers.bedpe"), false, Feature2DList.ListFormat.NA);
    }

    private static void buildFromOffsets(String[] args, CommandLineParser parser) {
        // generate-from-offsets <genomeID> <loops.bedpe> <output.bedpe> <+-offset>
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(4, usage);
        }
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        Feature2DList loops = Feature2DParser.loadFeatures(args[2], handler, false, null, false);
        String outFile = args[3];
        boolean addOffset = args[4].contains("+");
        boolean subtractOffset = args[4].contains("-");
        int offset = Integer.parseInt(args[4].replaceAll("\\+", "").replaceAll("-", ""));
        Feature2DList output = OffsetTools.createLoops(loops, offset, addOffset, subtractOffset);
        output.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }
}
