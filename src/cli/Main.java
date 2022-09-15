package cli;

import cli.clt.*;
import jargs.gnu.CmdLineParser;

public class Main {

    public static final String VERSION_NUM = "0.35.0";
    public static final int DEFAULT_RESOLUTION = 5000;
    public static final int DEFAULT_CUTOFF = 500;
    public static boolean printVerboseComments = false;

    public static void printGeneralUsageAndExit(int exitCode) {
        System.out.println("CLEAR Tools Version " + VERSION_NUM);
        System.out.println("Usage:");
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Commands:");
        System.out.println(Flags.usage);
        System.out.println(Pinpoint.usage);
        System.out.println(Cleaner.usage);
        System.out.println(APA.usage);
        System.out.println(ATA.usage);
        System.out.println(Recap.usage);
        System.out.println(Sieve.usage);
        System.out.println(HotSpot.usage);
        System.out.println(Fusion.usage);
        System.out.println(Sift.usage);
        System.out.println(NormHack.usage);
        System.out.println(SimplePeak.usage);
        System.out.println(SimpleMax.usage);
        System.out.println(GenerateBedpe.usage);
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
        } else if (command.startsWith("split") || command.startsWith("join")) {
            new SplitOrJoin(command, args);
        } else if (command.startsWith("sieve")) {
            new Sieve(args, parser, command);
        } else if (command.startsWith("hotspot")) {
            HotSpot.run(args, parser);
        } else if (command.startsWith("sift")) {
            new Sift(args, parser);
        } else if (command.startsWith("fuse") || command.startsWith("fusion")) {
            new Fusion(args, command);
        } else if (command.startsWith("seer")) {
            Seer.run(args, parser);
        } else if (command.startsWith("hack")) {
            NormHack.run(args, parser);
        } else if (command.startsWith("random")) {
            RandomLoops.run(args, parser);
        } else if (command.startsWith("generate")) {
            GenerateBedpe.run(args, parser);
        } else if (command.startsWith("simple")) {
            if (command.contains("max")) {
                SimpleMax.run(args, parser);
            } else {
                SimplePeak.run(args, parser);
            }
        } else {
            printGeneralUsageAndExit(3);
        }
    }
}
