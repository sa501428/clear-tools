package flags;

import flags.apa.APA;
import flags.apa.Anchor;
import jargs.gnu.CmdLineParser;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

public class Main {

    public static final String VERSION_NUM = "0.2.0";
    public static final int DEFAULT_RESOLUTION = 5000;
    public static final int DEFAULT_CUTOFF = 500;
    public static final String DEFAULT_NORMALIZATION = "SCALE";
    public static boolean printVerboseComments = false;

    private static void printGeneralUsageAndExit() {
        System.out.println("Hi-C FLAGS Version " + VERSION_NUM);
        System.out.println("Usage:");
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Usage: flags.jar [--cutoff int] [--res int] [--norm string] " +
                "<hic_file> <bed_file> <out_folder>");
        System.exit(0);
    }

    public static void main(String[] argv) throws CmdLineParser.UnknownOptionException, CmdLineParser.IllegalOptionValueException {
        if (argv.length == 0 || argv[0].equals("-h") || argv[0].equals("--help") || argv[0].equals("-V") || argv[0].equals("--version")) {
            printGeneralUsageAndExit();
        }

        CommandLineParser parser = new CommandLineParser();
        parser.parse(argv);
        boolean help = parser.getHelpOption();
        boolean version = parser.getVersionOption();
        printVerboseComments = printVerboseComments || parser.getVerboseOption();

        String[] args = parser.getRemainingArgs();
        if(args.length != 3 || help || version){
            printGeneralUsageAndExit();
        }

        actualRun(args, parser.getResolutionOption(), parser.getCutoffOption(), parser.getNormalizationStringOption());
    }

    public static void actualRun(String[] args, int resolution, int cutoff, String normString) {

        String hicfile = args[0];
        String bedfile = args[1];
        String outfolder = args[2];

        Dataset ds = HiCFileTools.extractDatasetForCLT(hicfile, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();
        
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normString, "KR", "SCALE", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        GenomeWide1DList<Anchor> anchors = GenericLocusParser.loadFromBEDFile(handler, bedfile, cutoff);
        System.out.println("Number of anchors: " + anchors.size());
        APA apa = new APA(ds, outfolder, norm, anchors, resolution);
        apa.run();
    }
}
