package cli.clt.sieve;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.general.BedGraphParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class SieveBedgraph {

    public static String usage = "sieve-bedgraph [-r resolution] <loops.bedpe> <out.stem> <genomeID> <file.bedgraph> [file2.bedgraph]";

    public int resolution = 1000;

    public SieveBedgraph(String[] args, CommandLineParser parser, String command) {
        // sieve <loops.bedpe> <output.bedpe> <file1.hic> <res1,res2,...>
        if (args.length != 5 && args.length != 6) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        String loopListPath = args[1];
        String outStem = args[2];
        String genomeID = args[3];
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        resolution = parser.getResolutionOption(resolution);

        Map<String, float[]> bedgraphData = BedGraphParser.parse(handler, args[4], resolution);

        Feature2DList loopList = Feature2DParser.loadFeatures(loopListPath, handler, true, null, false);
        Feature2DList[] result;

        if (args.length == 6) {
            Map<String, float[]> bedgraphData2 = BedGraphParser.parse(handler, args[5], resolution);
            result = intersect(loopList, bedgraphData, bedgraphData2, resolution, outStem);
        } else {
            result = intersect(loopList, bedgraphData, bedgraphData, resolution, outStem);
        }

        String outPath = loopListPath.replaceAll(".bedpe", outStem + ".anchor.bedgraph.bedpe");
        result[0].exportFeatureList(new File(outPath), false, Feature2DList.ListFormat.NA);
        outPath = loopListPath.replaceAll(".bedpe", outStem + ".good.anchor.bedgraph.bedpe");
        result[1].exportFeatureList(new File(outPath), false, Feature2DList.ListFormat.NA);
        outPath = loopListPath.replaceAll(".bedpe", outStem + ".bad.anchor.bedgraph.bedpe");
        result[2].exportFeatureList(new File(outPath), false, Feature2DList.ListFormat.NA);

        System.out.println("sieve-bedgraph complete");
    }

    private Feature2DList[] intersect(Feature2DList loopList, Map<String, float[]> bedgraphData,
                                      Map<String, float[]> bedgraphData2, int resolution, String outStem) {
        Feature2DList result = new Feature2DList();
        Feature2DList good = new Feature2DList();
        Feature2DList bad = new Feature2DList();


        loopList.processLists((chromKey, list) -> {

            List<Feature2D> outList = new LinkedList<>();
            List<Feature2D> goodList = new LinkedList<>();
            List<Feature2D> badList = new LinkedList<>();

            for (Feature2D loop : list) {

                float upVal = getValAtMidpoint(bedgraphData.get(loop.getChr1()), loop.getStart1(), loop.getEnd1(), resolution);
                float downVal = getValAtMidpoint(bedgraphData2.get(loop.getChr2()), loop.getStart2(), loop.getEnd2(), resolution);

                loop.addFloatAttribute(outStem + "_upstream_anchor_score", upVal);
                loop.addFloatAttribute(outStem + "_downstream_anchor_score", downVal);

                outList.add(loop);
                if (upVal > 2 && downVal > 2
                        && isEnriched(bedgraphData.get(loop.getChr1()), loop.getStart1(), loop.getEnd1(), resolution)
                        && isEnriched(bedgraphData2.get(loop.getChr2()), loop.getStart2(), loop.getEnd2(), resolution)) {
                    goodList.add(loop);
                } else {
                    badList.add(loop);
                }
            }

            result.addByKey(chromKey, outList);
            good.addByKey(chromKey, goodList);
            bad.addByKey(chromKey, badList);
        });

        return new Feature2DList[]{result, good, bad};
    }

    private boolean isEnriched(float[] data, long start, long end, int resolution) {
        if (data != null && data.length > 0) {
            int index = (int) (((start + end) / 2) / resolution);
            if (index > 0 && index < data.length - 1) {
                return data[index] > data[index - 1] && data[index] > data[index + 1];
            }
        }
        return false;
    }

    private float getValAtMidpoint(float[] data, long start, long end, int resolution) {
        if (data != null && data.length > 0) {
            return data[(int) (((start + end) / 2) / resolution)];
        }
        return 0;
    }
}
