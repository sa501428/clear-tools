package cli.clt.bedpe;

import cli.Main;
import cli.utils.general.FusionTools;

public class Fusion {
    public static String usage = "fuse[-nms][-clean][-exact][-file-id] <genomeID> <output.bedpe> <file1.bedpe> " +
            "<file2.bedpe> [...files.bedpe]\n" +
            "\t\tdefault behavior combines and coalesces overlapping features into a bounding box\n" +
            "\t\tnms retains high res feature and removes anything overlapping it\n" +
            "\t\tclean avoids saving old attributes\n" +
            "\t\texact uses precise bounds and requires exact matches, essentially a duplicate remover\n" +
            "\t\tfile-id saves an index for each input file";

    public static void run(String[] args, String command) {
        if (args.length < 4) {
            Main.printGeneralUsageAndExit(5, usage);
        }
        // fusion <genomeID> <output.bedpe> <file1.bedpe> <file2.bedpe> [...files.bedpe]
        String genomeID = args[1];
        String outFile = args[2];
        String[] bedpeFiles = new String[args.length - 3];
        System.arraycopy(args, 3, bedpeFiles, 0, bedpeFiles.length);
        FusionTools.coalesceFeatures(bedpeFiles, genomeID, outFile,
                command.contains("nms"), command.contains("clean"), command.contains("exact"),
                command.contains("file-id"));
        System.out.println("fusion complete");
    }
}
