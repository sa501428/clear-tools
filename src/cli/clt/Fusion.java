package cli.clt;

import cli.Main;
import cli.utils.general.FusionTools;

public class Fusion {
    public Fusion(String[] args, CommandLineParser parser, boolean useNMS) {
        if (args.length < 5) {
            Main.printGeneralUsageAndExit(5);
        }
        // fusion <genomeID> <output.bedpe> <file1.bedpe> <file2.bedpe> [...files.bedpe]
        String genomeID = args[1];
        String outFile = args[2];
        String[] bedpeFiles = new String[args.length - 3];
        System.arraycopy(args, 3, bedpeFiles, 0, bedpeFiles.length);
        FusionTools.coalesceFeatures(bedpeFiles, genomeID, outFile, useNMS);
        System.out.println("fusion complete");
    }
}
