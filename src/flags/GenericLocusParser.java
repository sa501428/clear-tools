package flags;

import flags.apa.Anchor;
import javastraw.StrawGlobals;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

public class GenericLocusParser {

    public static GenomeWide1DList<Anchor> loadFromBEDFile(ChromosomeHandler handler,
                                                           String bedFilePath, float cutoff) {
        List<Anchor> anchors = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(bedFilePath)),
                    StrawGlobals.bufferSize);
            anchors.addAll(parseBEDFile(br, handler, cutoff));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWide1DList<>(handler, anchors);
    }

    private static List<Anchor> parseBEDFile(BufferedReader bufferedReader,
                                             ChromosomeHandler handler,
                                             float scoreCutoff) throws IOException {
        Set<Anchor> anchors = new HashSet<>();
        String nextLine;

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            String chr1Name;
            int start1, end1;

            if (tokens[0].startsWith("chr") && tokens.length > 4) {
                // valid line
                chr1Name = tokens[0];
                start1 = Integer.parseInt(tokens[1]);
                end1 = Integer.parseInt(tokens[2]);
                float score = Float.parseFloat(tokens[4]);

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

                if (score > scoreCutoff) {
                    anchors.add(new Anchor(chr.getName(), chr.getIndex(), start1, end1));
                }
            }
        }
        if (anchors.size() < 1) System.err.println("BED File empty - file may have problems or error was encountered");
        bufferedReader.close();
        return new ArrayList<>(anchors);
    }
}
