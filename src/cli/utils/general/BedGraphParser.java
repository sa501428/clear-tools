package cli.utils.general;

import javastraw.StrawGlobals;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class BedGraphParser {
    public static Map<String, float[]> parse(ChromosomeHandler handler, String path, int resolution) {
        Map<String, float[]> map = null;
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(path)),
                    StrawGlobals.bufferSize);
            map = parseBedGraphFile(br, handler, resolution);
            br.close();
        } catch (IOException ec) {
            ec.printStackTrace();
        }
        return map;
    }

    private static Map<String, float[]> parseBedGraphFile(BufferedReader br, ChromosomeHandler handler, int resolution) throws IOException {
        Map<String, float[]> map = new HashMap<>();
        String nextLine;
        while ((nextLine = br.readLine()) != null) {

            String[] tokens = nextLine.split("\\s+");
            String chrName;
            int start, end;
            if (tokens[0].startsWith("chr") && tokens.length > 3) {
                if (tokens[0].equals("chrom") && tokens[1].equals("x1")) continue;
                chrName = tokens[0];
                if (!map.containsKey(chrName)) {
                    Chromosome chromosome = handler.getChromosomeFromName(chrName);
                    float[] values = new float[(int) (chromosome.getLength() / resolution) + 1];
                    map.put(chrName, values);
                }

                start = Integer.parseInt(tokens[1]) / resolution;
                end = Integer.parseInt(tokens[2]) / resolution;

                float[] values = map.get(chrName);
                float value = Float.parseFloat(tokens[3]);
                for (int i = start; i < end; i++) {
                    values[i] = value;
                }
            }
        }
        return map;
    }
}
