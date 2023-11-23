package cli.clt.bedpe;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.general.FusionTools;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.*;

public class Fusion {
    public static String usage = "fuse[-nms][-clean][-exact][-file-id] [--round val] " +
            "[--attributes a,b,c,...] <genomeID> <output.bedpe> <file1.bedpe> " +
            "<file2.bedpe> [...files.bedpe]\n" +
            "\t\tdefault behavior combines and coalesces overlapping features into a bounding box\n" +
            "\t\tnms retains high res feature and removes anything overlapping it\n" +
            "\t\tclean avoids saving old attributes\n" +
            "\t\texact uses precise bounds and requires exact matches, essentially a duplicate remover\n" +
            "\t\tfile-id saves an index for each input file\n" +
            "\t\t--attributes allows for the saving of specific attributes in the final list\\n\" +\n" +
            "\t\t--round rounds down to the nearest multiple of val";

    public static void run(String[] args, String command, CommandLineParser parser) {
        if (args.length < 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }
        // fusion <genomeID> <output.bedpe> <file1.bedpe> <file2.bedpe> [...files.bedpe]
        String genomeID = args[1];
        String outFile = args[2];
        String[] bedpeFiles = new String[args.length - 3];
        System.arraycopy(args, 3, bedpeFiles, 0, bedpeFiles.length);
        String[] attributes = parser.getAttributesOption();
        int val = parser.getRoundOption();

        if (command.contains("clean") && command.contains("exact")
                && (!command.contains("nms")) && (!command.contains("file-id"))) {

            simpleDedup(bedpeFiles, genomeID, outFile);

        } else {
            FusionTools.coalesceFeatures(bedpeFiles, genomeID, outFile,
                    command.contains("nms"), command.contains("clean"), command.contains("exact"),
                    command.contains("file-id"), attributes, val);
        }
        System.out.println("fusion complete");
    }

    private static void simpleDedup(String[] bedpeFiles, String genomeID, String outFile) {
        Map<String, Set<Feature2D>> dedupedMapping = new HashMap<>();
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);

        for (String path : bedpeFiles) {
            Feature2DList loopList = Feature2DParser.loadFeatures(path, handler,
                    false, null, false);
            loopList.processLists((key, list) -> {
                if (!dedupedMapping.containsKey(key)) {
                    dedupedMapping.put(key, new HashSet<>());
                }
                dedupedMapping.get(key).addAll(list);
            });
        }
        Feature2DList deduped = new Feature2DList();
        for (String key : dedupedMapping.keySet()) {
            deduped.addByKey(key, new ArrayList<>(dedupedMapping.get(key)));
        }
        deduped.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }
}
