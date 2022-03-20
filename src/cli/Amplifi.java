package cli;

public class Amplifi {
    public static void run(String[] args, int resolution, String normString) {
        if(args.length != 4){
            Main.printGeneralUsageAndExit();
        }

        String hicFiles = args[1];
        String bedpeFile = args[2];
        String outFolder = args[3];

        /*
        Dataset ds = HiCFileTools.extractDatasetForCLT(hicFile, false, true);
        ChromosomeHandler handler = ds.getChromosomeHandler();

        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{normString, "KR", "SCALE", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        GenomeWide1DList<Anchor> anchors = GenericLocusParser.loadFromBEDFile(handler, bedpeFile, cutoff);
        System.out.println("Number of anchors: " + anchors.size());
        APA apa = new APA(ds, outFolder, norm, anchors, resolution);
        apa.run();
        */

    }
}
