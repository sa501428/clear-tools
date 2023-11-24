package cli.clt.misc;

import cli.Main;
import cli.clt.CommandLineParser;
import cli.utils.general.BedFileParser;
import cli.utils.general.BedGraphParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class IntersectBedWithBedgraph {

    public static String usage = "intersect-bed-bedgraph [-r resolution] <genomeID> <anchors.bed> <file.bedgraph> " +
            "<out.stem>";

    public int resolution = 1000;

    public IntersectBedWithBedgraph(String[] args, CommandLineParser parser, String command) {
        // sieve <loops.bedpe> <output.bedpe> <file1.hic> <res1,res2,...>
        if (args.length != 5 && args.length != 6) {
            Main.printGeneralUsageAndExit(5, usage);
        }

        String genomeID = args[1];
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);

        String anchorsBed = args[2];
        Map<String, List<int[]>> anchors = BedFileParser.simpleParser(anchorsBed);

        resolution = parser.getResolutionOption(resolution);

        Map<String, float[]> bedgraphData = BedGraphParser.parse(handler, args[3], resolution);
        String outStem = args[4];

        Map<String, List<IntervalEntry>> result = intersect(anchors, bedgraphData, resolution);
        bedgraphData.clear();
        anchors.clear();

        try {
            exportToBedgraph(outStem + ".bedgraph", result);
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("intersect-bed-bedgraph complete");
    }

    public static void exportToBedgraph(String path, Map<String, List<IntervalEntry>> result) throws IOException {
        File outputFileName = new File(path);
        outputFileName.createNewFile();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
        // write for every chromosome
        for (String chromosome : result.keySet()) {
            for (IntervalEntry interval : result.get(chromosome)) {
                bw.write(chromosome + " " + interval.start + " " + interval.end + " " + interval.value);
                bw.newLine();
            }
        }
        bw.close();
    }

    private Map<String, List<IntervalEntry>> intersect(Map<String, List<int[]>> anchors, Map<String, float[]> bedgraphData, int resolution) {

        Map<String, List<IntervalEntry>> result = new HashMap<>();
        for (String chr : anchors.keySet()) {
            result.put(chr, new LinkedList<>());
        }

        for (String chr : anchors.keySet()) {
            List<int[]> anchorList = anchors.get(chr);
            float[] data = bedgraphData.get(chr);
            for (int[] anchor : anchorList) {
                long start = anchor[0];
                long end = anchor[1];
                float val = getValAtMidpoint(data, start, end, resolution);
                result.get(chr).add(new IntervalEntry(start, end, val));
            }
        }

        return result;
    }

    private float getValAtMidpoint(float[] data, long start, long end, int resolution) {
        if (data != null && data.length > 0) {
            return data[(int) (((start + end) / 2) / resolution)];
        }
        return 0;
    }
}
