package cli;

import cli.clt.*;
import jargs.gnu.CmdLineParser;

public class Main {

    public static final String VERSION_NUM = "0.26.0";
    public static final int DEFAULT_RESOLUTION = 5000;
    public static final int DEFAULT_CUTOFF = 500;
    public static boolean printVerboseComments = false;

    public static void printGeneralUsageAndExit(int exitCode) {
        System.out.println("CLEAR Tools Version " + VERSION_NUM);
        System.out.println("Usage:");
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Commands: \n" +
                "flags [--cutoff int] [--res int] [--norm string] <input.hic> <loops.bedpe> <out_folder>\n" +
                "enhance [--res int] [--norm string] <out_folder> <loops.bedpe> <hic_files>\n" +
                "probability [--res int] <input.hic> <loops.bedpe> <out_folder>\n" +
                "pinpoint [--res int] <input.hic> <loops.bedpe> <output.bedpe>\n" +
                "clean <input.hic> <loops.bedpe> <output.bedpe>\n" +
                "apa [options] <input.hic> <loops.bedpe> <outfolder>\n" +
                "ata [--res int] <signal.bw> <peaks.bed> <outfile> <genome>\n" +
                "recap [--loop] <loops.bedpe> <outfolder> <file1.hic,file2.hic,...> <name1,name2,...>\n" +
                "hotspot [--res int] [--norm string] <file1.hic,file2.hic,...> <outfile>\n" +
                "fuse <genomeID> <output.bedpe> <file1.bedpe> <file2.bedpe> [...files.bedpe]\n" +
                "sift [--window int] [--min double] [--max double] [--res int] [--low-res int] <file.hic> <outfile>\n" +
                "seer [--res int] [--low-res int] [-k norm] [--seed seed] <file> <out_folder> <num_contacts>\n" +
                "hack [--res int] <out_folder> <file1.hic,file2.hic,...> <name1,name2,...>\n" +
                "simple-max [-r resolution] [-k norm] <file.hic> <loops.bedpe> <output.bedpe>");
        System.out.println("Exit code " + exitCode);
        System.exit(exitCode);
    }

    public static void main(String[] argv) throws CmdLineParser.UnknownOptionException, CmdLineParser.IllegalOptionValueException {
        if (argv.length == 0 || argv[0].equals("-h") || argv[0].equals("--help") || argv[0].equals("-V") || argv[0].equals("--version")) {
            printGeneralUsageAndExit(1);
        }

        CommandLineParser parser = new CommandLineParser();
        parser.parse(argv);
        boolean help = parser.getHelpOption();
        boolean version = parser.getVersionOption();
        printVerboseComments = printVerboseComments || parser.getVerboseOption();

        String[] args = parser.getRemainingArgs();
        if(help || version){
            printGeneralUsageAndExit(2);
        }

        String command = args[0].toLowerCase();
        if(command.equals("flags")){
            Flags.run(args, parser.getResolutionOption(Main.DEFAULT_RESOLUTION), parser.getCutoffOption(), parser.getNormalizationStringOption());
        } else if (command.equals("enhance") || command.equals("amplifi") || command.equals("amplify")) {
            Enhance.run(args, parser.getResolutionOption(Main.DEFAULT_RESOLUTION), parser.getNpyOption());
        } else if (command.equals("pinpoint")) {
            Pinpoint.run(args, parser);
        } else if (command.startsWith("clean")) {
            Cleaner.run(args);
        } else if (command.startsWith("prob")) {
            Probability.run(args, parser.getResolutionOption(Main.DEFAULT_RESOLUTION), parser.getLogOption());
        } else if (command.startsWith("apa")) {
            APA apa = new APA(args, parser);
            apa.run();
        } else if (command.startsWith("ata")) {
            ATA ata = new ATA(args, parser);
            ata.run();
        } else if (command.startsWith("recap") || command.startsWith("compile")) {
            new Recap(args, parser);
        } else if (command.startsWith("hotspot")) {
            HotSpot.run(args, parser);
        } else if (command.startsWith("sift")) {
            new Sift(args, parser);
        } else if (command.startsWith("fuse") || command.startsWith("fusion")) {
            new Fusion(args, parser);
        } else if (command.startsWith("seer")) {
            Seer.run(args, parser);
        } else if (command.startsWith("hack")) {
            NormHack.run(args, parser);
        } else if (command.startsWith("random")) {
            RandomLoops.run(args, parser);
        } else if (command.startsWith("generate")) {
            GenerateBedpe.run(args, parser);
        } else if (command.startsWith("simple")) {
            SimpleMax.run(args, parser);
        } else {
            printGeneralUsageAndExit(3);
        }
    }
}
