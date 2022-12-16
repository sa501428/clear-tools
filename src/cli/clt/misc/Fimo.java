package cli.clt.misc;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class Fimo {

    // fimo <file.tsv>

    private final static String DATABASE = "_HUMAN.H11MO";

    public static void run(String[] args, String command) {

        Map<String, Map<String, List<String[]>>> data = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(args[1]))) {
            for (String line; (line = br.readLine()) != null; ) {
                try {
                    if (line.contains(DATABASE)) {
                        String[] row = line.split("\t");
                        String motif = cleanMotifName(row[0]);
                        if (!data.containsKey(motif)) {
                            data.put(motif, getNewMapForMotif());
                        }
                        data.get(motif).get(row[5]).add(new String[]{row[2], row[3], row[4]});
                    }
                } catch (Exception e) {
                    // ignore problematic lines like header/footer
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        for (String motif : data.keySet()) {
            for (String strand : data.get(motif).keySet()) {

                try {

                    File fout = new File(motif + getStrandName(strand) + ".bed");
                    FileOutputStream fos = new FileOutputStream(fout);
                    BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));

                    for (String[] entry : data.get(motif).get(strand)) {
                        bw.write(entry[0] + "\t" + entry[1] + "\t" + entry[2]);
                        bw.newLine();
                    }

                    bw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private static Map<String, List<String[]>> getNewMapForMotif() {
        Map<String, List<String[]>> data = new HashMap<>();
        data.put("+", new LinkedList<>());
        data.put("-", new LinkedList<>());
        return data;
    }

    private static String getStrandName(String strand) {
        if (strand.equals("+")) {
            return "_pos";
        } else if (strand.equals("-")) {
            return "_neg";
        } else {
            return null;
        }
    }

    private static String cleanMotifName(String param) {
        return param.split(DATABASE)[0];
    }
}
