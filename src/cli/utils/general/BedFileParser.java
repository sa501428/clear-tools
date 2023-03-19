package cli.utils.general;

import cli.utils.flags.Anchor;
import cli.utils.flags.AnchorWithScore;
import javastraw.StrawGlobals;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

public class BedFileParser {

    public static GenomeWide1DList<Anchor> loadFromBEDFile(ChromosomeHandler handler,
                                                           String bedFilePath, float cutoff, boolean isPercentile,
                                                           boolean saveScore) {
        List<Anchor> anchors = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(bedFilePath)),
                    StrawGlobals.bufferSize);
            anchors.addAll(parseBEDFile(br, handler, cutoff, isPercentile));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWide1DList<>(handler, anchors);
    }

    private static List<Anchor> parseBEDFile(BufferedReader bufferedReader,
                                             ChromosomeHandler handler,
                                             float scoreCutoff, boolean isPercentile) throws IOException {

        boolean useCutoff = scoreCutoff > 0;

        Set<Anchor> anchors = new HashSet<>();
        String nextLine;
        DescriptiveStatistics statistics = new DescriptiveStatistics();

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            String chr1Name;
            int start1, end1;

            if (tokens[0].startsWith("chr") && tokens.length > 2) {
                if (tokens[0].equals("chrom") && tokens[1].equals("x1")) continue;

                // valid lines should have at least 3 columns
                chr1Name = tokens[0];
                start1 = Integer.parseInt(tokens[1]);
                end1 = Integer.parseInt(tokens[2]);

                Chromosome chr = handler.getChromosomeFromName(chr1Name);
                if (chr == null) {
                    if (errorCount < 10) {
                        System.err.println("Skipping line: " + nextLine);
                    } else if (errorCount == 10) {
                        System.err.println("Maximum error count exceeded.  Further errors will not be logged");
                    }
                    errorCount++;
                    continue;
                }

                if (tokens.length < 4) {
                    anchors.add(new Anchor(chr.getName(), start1, end1, chr.getIndex()));
                } else {
                    String name = tokens[3];
                    float score = Float.parseFloat(tokens[4]);

                    if (isPercentile) {
                        statistics.addValue(score);
                        anchors.add(new AnchorWithScore(chr.getName(), start1, end1, score, chr.getIndex(),
                                name));
                    } else if (useCutoff) {
                        if (score > scoreCutoff) {
                            anchors.add(new Anchor(chr.getName(), start1, end1, chr.getIndex()));
                        }
                    } else {
                        anchors.add(new AnchorWithScore(chr.getName(), start1, end1, score, chr.getIndex(),
                                name));
                    }
                }
            }
        }

        bufferedReader.close();
        if (anchors.size() < 1) System.err.println("BED File empty - file may have problems or error was encountered");
        if (useCutoff && isPercentile) {
            List<Anchor> anchorsThatPass = new ArrayList<>();
            double cutoff = statistics.getPercentile(scoreCutoff);
            for (Anchor anchor : anchors) {
                AnchorWithScore anchor1 = (AnchorWithScore) anchor;
                if (anchor1.getScore() >= cutoff) {
                    anchorsThatPass.add(anchor);
                }
            }
            return anchorsThatPass;
        }
        return new ArrayList<>(anchors);
    }
}
