package cli;

import jargs.gnu.CmdLineParser;

public class Main {

    public static final String VERSION_NUM = "0.3.1";
    public static final int DEFAULT_RESOLUTION = 5000;
    public static final int DEFAULT_CUTOFF = 500;
    public static final String DEFAULT_NORMALIZATION = "SCALE";
    public static boolean printVerboseComments = false;

    public static void printGeneralUsageAndExit(int exitCode) {
        System.out.println("Hi-C FLAGS Version " + VERSION_NUM);
        System.out.println("Usage:");
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Commands: \n" +
                "flags [--cutoff int] [--res int] [--norm string] <hic_file> <bed_file> <out_folder>\n" +
                "amplifi [--res int] [--norm string] <out_folder> <bedpe_file> <hic_files>");
        System.out.println("Exit code "+exitCode);
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
            Flags.run(args, parser.getResolutionOption(), parser.getCutoffOption(), parser.getNormalizationStringOption());
        } else if (command.equals("enhance") || command.equals("amplifi") || command.equals("amplify")){
            Enhance.run(args, parser.getResolutionOption(), parser.getNormalizationStringOption());
        } else {
            printGeneralUsageAndExit(3);
        }
    }


}
