package cli.clt.loops;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.general.BedGraphParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.Map;

public class SieveBedgraph {

    public static String usage = "sieve-bedgraph [-r resolution] <loops.bedpe> <out.stem> <file.bedgraph> <genomeID>";

    public int resolution = 1000;

    public SieveBedgraph(String[] args, CommandLineParser parser, String command) {
        // sieve <loops.bedpe> <output.bedpe> <file1.hic> <res1,res2,...>
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        String loopListPath = args[1];
        String outStem = args[2];
        String bedgraphPath = args[3];
        String genomeID = args[4];

        resolution = parser.getResolutionOption(resolution);

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);

        Map<String, float[]> bedgraphData = BedGraphParser.parse(handler, bedgraphPath, resolution);
        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler, true, null, false);

        Feature2DList result = intersect(loopList, bedgraphData, resolution, outStem);
        result.exportFeatureList(new File(outStem + ".anchor.bedgraph.bedpe"), false, Feature2DList.ListFormat.NA);

        System.out.println("sieve-bedgraph complete");
    }

    private Feature2DList intersect(Feature2DList loopList, Map<String, float[]> bedgraphData, int resolution, String outStem) {
        Feature2DList result = new Feature2DList();
        loopList.processLists((chromKey, list) -> {

            for (Feature2D loop : list) {

                float upVal = getValAtMidpoint(bedgraphData.get(loop.getChr1()), loop.getStart1(), loop.getEnd1(), resolution);
                float downVal = getValAtMidpoint(bedgraphData.get(loop.getChr2()), loop.getStart2(), loop.getEnd2(), resolution);

                loop.addFloatAttribute(outStem + "_upstream_anchor_score", upVal);
                loop.addFloatAttribute(outStem + "_downstream_anchor_score", downVal);
            }

        });

        return result;
    }

    private float getValAtMidpoint(float[] data, long start, long end, int resolution) {
        if (data != null && data.length > 0) {
            return data[(int) (((start + end) / 2) / resolution)];
        }
        return 0;
    }
}
