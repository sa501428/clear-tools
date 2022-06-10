package cli.utils;

import cli.utils.flags.RegionConfiguration;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.type.HiCZoom;

import java.util.List;
import java.util.Map;

public class HiCUtils {
    public static int populateChromosomePairs(Map<Integer, RegionConfiguration> chromosomePairs,
                                              Chromosome[] chromosomes, boolean includeInter) {
        int pairCounter = 0;
        for (int i = 0; i < chromosomes.length; i++) {
            if (includeInter) {
                for (int j = i; j < chromosomes.length; j++) {
                    RegionConfiguration config = new RegionConfiguration(chromosomes[i], chromosomes[j]);
                    chromosomePairs.put(pairCounter, config);
                    pairCounter++;
                }
            } else {
                RegionConfiguration config = new RegionConfiguration(chromosomes[i], chromosomes[i]);
                chromosomePairs.put(pairCounter, config);
                pairCounter++;
            }
        }
        return pairCounter;
    }

    public static HiCZoom getHighestResolution(List<HiCZoom> zooms) {
        HiCZoom highest = zooms.get(0);
        for (HiCZoom zoom : zooms) {
            if (zoom.getBinSize() < highest.getBinSize()) {
                highest = zoom;
            }
        }
        return highest;
    }
}
