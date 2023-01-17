package cli.clt.loops;

import cli.clt.CommandLineParser;
import cli.utils.motifs.IndexedBedFile;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class MotifAssignment {

    public static String usage = "assign-motifs [--window val] <genomeID> <loops.bedpe> " +
            "<upstream.motifs.bed> <downstream.motifs.bed> <output.bedpe>";
    protected final int window;
    private final ChromosomeHandler handler;
    private final Feature2DList loopList;
    private final int binSize;
    private final String outFile;
    private Map<Integer, Map<Integer, List<int[]>>> upBed;
    private Map<Integer, Map<Integer, List<int[]>>> downBed;

    public MotifAssignment(String[] args, CommandLineParser parser) {
        if (args.length != 6) {
            System.out.println(usage);
            System.exit(6);
        }

        window = parser.getWindowSizeOption(250);
        binSize = 3 * window;
        handler = ChromosomeTools.loadChromosomes(args[1]);
        loopList = Feature2DParser.loadFeatures(args[2], handler,
                true, null, false);
        outFile = args[5];
        try {
            upBed = IndexedBedFile.index(args[3], handler, binSize);
            downBed = IndexedBedFile.index(args[4], handler, binSize);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(7);
        }
    }

    public void run() {
        Feature2DList result = new Feature2DList();
        for (Chromosome chromosome : handler.getChromosomeArrayWithoutAllByAll()) {
            List<Feature2D> chrLoops = loopList.get(chromosome.getIndex(), chromosome.getIndex());
            List<Feature2D> loopsToSave = new ArrayList<>();
            Map<Integer, List<int[]>> upMotifs = upBed.get(chromosome.getIndex());
            Map<Integer, List<int[]>> downMotifs = downBed.get(chromosome.getIndex());

            for (Feature2D loop : chrLoops) {
                int[] upMotif = IndexedBedFile.getUniqueMotif(loop.getAttribute("localX"), upMotifs, binSize, window);
                int[] downMotif = IndexedBedFile.getUniqueMotif(loop.getAttribute("localY"), downMotifs, binSize, window);
                if (upMotif != null && downMotif != null) {
                    IndexedBedFile.setMotifAttributes(loop, upMotif, true);
                    IndexedBedFile.setMotifAttributes(loop, downMotif, false);
                    loopsToSave.add(loop);
                }
            }

            result.addByKey(Feature2DList.getKey(chromosome, chromosome), loopsToSave);
        }
        result.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }
}
