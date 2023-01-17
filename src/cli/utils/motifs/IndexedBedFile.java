package cli.utils.motifs;

import cli.utils.general.IGVTools;
import javastraw.feature2D.Feature2D;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.track.TribbleFeatureSource;

import java.io.IOException;
import java.util.*;

public class IndexedBedFile {

    private static final int MIDPT = 0;
    private static final int START = 1;
    private static final int END = 2;

    public static Map<Integer, Map<Integer, List<int[]>>> index(String inputBedFile, ChromosomeHandler handler, int binSize) throws IOException, TribbleIndexNotFoundException {
        TribbleFeatureSource bedFile = IGVTools.loadBed(inputBedFile, handler);
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        Map<Integer, Map<Integer, List<int[]>>> bedMap = new HashMap<>();
        for (Chromosome chromosome : chromosomes) {
            Map<Integer, List<int[]>> binToMotif = new HashMap<>();
            Iterator<?> iterator = bedFile.getFeatures(chromosome.getName(), 0, (int) chromosome.getLength());
            while (iterator.hasNext()) {
                IGVFeature feature = (IGVFeature) iterator.next();
                int start = feature.getStart();
                int end = feature.getEnd();
                int mid = (start + end) / 2;
                binToMotif.computeIfAbsent(mid / binSize,
                        k -> new LinkedList<>()).add(new int[]{mid, start, end});
            }
            bedMap.put(chromosome.getIndex(), binToMotif);
        }
        return bedMap;
    }

    public static int[] getUniqueMotif(String localPos, Map<Integer, List<int[]>> binnedMotifs,
                                       int binSize, int window) {
        try {
            long x = Long.parseLong(localPos);
            int bin = (int) (x / binSize);

            List<int[]> motifs = new ArrayList<>();
            for (int i = bin - 1; i < bin + 2; i++) {
                if (binnedMotifs.containsKey(i)) {
                    for (int[] motif : binnedMotifs.get(i)) {
                        if (Math.abs(motif[MIDPT] - x) < window) {
                            motifs.add(motif);
                        }
                    }
                }
            }
            if (motifs.size() == 1) {
                return motifs.get(0);
            }
        } catch (Exception e) {
            return null;
        }
        return null;
    }

    public static void setMotifAttributes(Feature2D loop, int[] motif, boolean isUpStream) {
        String startKey, endKey, midKey;
        if (isUpStream) {
            startKey = "motif_start_1";
            endKey = "motif_end_1";
            midKey = "motif_mid_1";
        } else {
            startKey = "motif_start_2";
            endKey = "motif_end_2";
            midKey = "motif_mid_2";
        }
        loop.setAttribute(startKey, "" + motif[START]);
        loop.setAttribute(endKey, "" + motif[END]);
        loop.setAttribute(midKey, "" + motif[MIDPT]);
    }
}
