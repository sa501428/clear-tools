package cli.clt;

import cli.Main;
import cli.utils.general.FusionTools;

public class Fusion {
    public static String usage = "fuse[-nms] <genomeID> <output.bedpe> <file1.bedpe> <file2.bedpe> [...files.bedpe]";

    public static void run(String[] args, String command) {
        if (args.length < 5) {
            Main.printGeneralUsageAndExit(5);
        }
        // fusion <genomeID> <output.bedpe> <file1.bedpe> <file2.bedpe> [...files.bedpe]
        String genomeID = args[1];
        String outFile = args[2];
        String[] bedpeFiles = new String[args.length - 3];
        System.arraycopy(args, 3, bedpeFiles, 0, bedpeFiles.length);
        FusionTools.coalesceFeatures(bedpeFiles, genomeID, outFile, command.contains("nms"));
        System.out.println("fusion complete");
    }
}
