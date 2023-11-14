package cli.utils.anchors;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.ParallelizationTools;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class AnchorFrequencyCounter {

    public static Map<String, int[]> calculateHiResAnchorCounts(List<Feature2DList> loopLists, ChromosomeHandler handler,
                                                                int highResolution) {
        Map<String, int[]> globalChromToHiResAnchorCounts = new HashMap<>();

        AtomicInteger globalIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            Map<String, int[]> localChromToHiResAnchorCounts = generateNewMap(handler, highResolution);
            int index = globalIndex.getAndIncrement();
            while (index < loopLists.size()) {
                Feature2DList loopList = loopLists.get(index);
                loopList.processLists((s, list) -> {
                    for (Feature2D loop : list) {
                        int x, y;
                        try {
                            x = (int) (Long.parseLong(loop.getAttribute("localX")) / highResolution);
                            y = (int) (Long.parseLong(loop.getAttribute("localY")) / highResolution);
                        } catch (Exception e) {
                            continue;
                        }
                        localChromToHiResAnchorCounts.get(loop.getChr1())[x]++;
                        localChromToHiResAnchorCounts.get(loop.getChr2())[y]++;
                    }
                });

                index = globalIndex.getAndIncrement();
            }

            synchronized (globalChromToHiResAnchorCounts) {
                for (String key : localChromToHiResAnchorCounts.keySet()) {
                    if (globalChromToHiResAnchorCounts.containsKey(key)) {
                        int[] localCounts = localChromToHiResAnchorCounts.get(key);
                        int[] globalCounts = globalChromToHiResAnchorCounts.get(key);
                        for (int k = 0; k < globalCounts.length; k++) {
                            globalCounts[k] += localCounts[k];
                        }
                    } else {
                        globalChromToHiResAnchorCounts.put(key, localChromToHiResAnchorCounts.get(key));
                    }
                }
            }
        });

        return globalChromToHiResAnchorCounts;
    }

    public static Map<String, int[]> generateNewMap(ChromosomeHandler handler, int highResolution) {
        Map<String, int[]> map = new HashMap<>();
        for (Chromosome chromosome : handler.getChromosomeArrayWithoutAllByAll()) {
            map.put(chromosome.getName(),
                    new int[(int) (chromosome.getLength() / highResolution) + 1]);
        }
        return map;
    }
}
