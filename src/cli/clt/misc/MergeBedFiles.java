package cli.clt.misc;

import cli.clt.CommandLineParser;
import cli.utils.general.BedFileParser;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class MergeBedFiles {

    public static String usage = "merge-bed-files <genomeID> <out.bed> <file1.bed> <file2.bed> ...";
    private final int resolution = 50;

    public MergeBedFiles(String[] args, CommandLineParser parser, String command) {
        // merge-bed-files <genomeID> <out.bed> <file1.bed> <file2.bed> ...
        if (args.length < 5) {
            System.err.println(usage);
            System.exit(58);
        }

        String genomeID = args[1];
        String outPath = args[2];

        String[] bedFiles = new String[args.length - 3];
        System.arraycopy(args, 3, bedFiles, 0, bedFiles.length);

        mergeBedFiles(genomeID, outPath, bedFiles);

        System.out.println("merge-bed-files complete");
    }

    private void mergeBedFiles(String genomeID, String outPath, String[] bedFiles) {
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Map<String, int[]> globalResults = createInitialMapping(handler);
        populateGlobalResults(globalResults, bedFiles);
        threshold(globalResults, 3);
        Map<String, List<int[]>> listOfBedFileRegions = extractBounds(globalResults);
        exportToBed(outPath, listOfBedFileRegions);
    }

    private void exportToBed(String outPath, Map<String, List<int[]>> listOfBedFileRegions) {
        File file = new File(outPath);
        try (PrintWriter writer = new PrintWriter(file.getPath())) {

            for (Map.Entry<String, List<int[]>> entry : listOfBedFileRegions.entrySet()) {
                String chr = entry.getKey();
                List<int[]> regions = entry.getValue();
                for (int[] region : regions) {
                    writer.println(chr + "\t" + region[0] * resolution + "\t" + (region[1] * resolution + resolution));
                }
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error exporting bed file (" + file.getPath() + "): " + e.getMessage());
        }
    }

    private Map<String, List<int[]>> extractBounds(Map<String, int[]> globalResults) {
        Map<String, List<int[]>> listOfBedFileRegions = new HashMap<>();
        for (Map.Entry<String, int[]> entry : globalResults.entrySet()) {
            List<int[]> regions = new LinkedList<>();
            int[] values = entry.getValue();
            int start = -1;
            int end = -1;
            for (int i = 0; i < values.length; i++) {
                if (values[i] > 0) {
                    if (start == -1) {
                        start = i;
                    }
                    end = i;
                } else {
                    if (start != -1) {
                        regions.add(new int[]{start, end});
                        start = -1;
                        end = -1;
                    }
                }
            }
            if (start != -1) {
                regions.add(new int[]{start, end});
            }
            listOfBedFileRegions.put(entry.getKey(), regions);
        }
        return listOfBedFileRegions;
    }

    private void threshold(Map<String, int[]> globalResults, int minValNeeded) {
        for (Map.Entry<String, int[]> entry : globalResults.entrySet()) {
            int[] values = entry.getValue();
            for (int i = 0; i < values.length; i++) {
                if (values[i] < minValNeeded) {
                    values[i] = 0;
                }
            }
        }
    }

    private void populateGlobalResults(Map<String, int[]> globalResults, String[] bedFiles) {
        for (String bedFile : bedFiles) {
            Map<String, List<int[]>> anchors = BedFileParser.simpleParser(bedFile);
            for (Map.Entry<String, List<int[]>> entry : anchors.entrySet()) {
                int[] global = globalResults.get(entry.getKey());
                for (int[] interval : entry.getValue()) {
                    for (int i = interval[0] / resolution; i < (interval[1] / resolution) + 1; i++) {
                        global[i]++;
                    }
                }
            }
        }
    }

    private Map<String, int[]> createInitialMapping(ChromosomeHandler handler) {
        Map<String, int[]> initialMap = new HashMap<>();
        for (Chromosome chrom : handler.getChromosomeArrayWithoutAllByAll()) {
            int chrLength = (int) (chrom.getLength() / resolution) + 1;
            initialMap.put(chrom.getName(), new int[chrLength]);
        }
        return initialMap;
    }
}
