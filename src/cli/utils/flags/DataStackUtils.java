package cli.utils.flags;

import javastraw.reader.basics.Chromosome;

import java.io.File;

public class DataStackUtils {

    public static FlagsDataStack[] initialize(Chromosome[] chromosomeArray, int matrixWidth, File outputDirectory,
                                              String prefix) {
        int maxNumIntraStacks = getMaxNumIntraStacks(chromosomeArray);
        FlagsDataStack[] stacks = new FlagsDataStack[maxNumIntraStacks];
        for (int z = 0; z < stacks.length; z++) {
            stacks[z] = new FlagsDataStack(matrixWidth, outputDirectory, prefix + "intra_" + z + "_");
        }
        return stacks;
    }

    private static int getMaxNumIntraStacks(Chromosome[] chromosomes) {
        int maxSizeInMB = (int) Math.ceil(getMaxSize(chromosomes) / 1000000.0);
        return (int) Math.ceil(Math.log(maxSizeInMB) / Math.log(2));
    }

    private static long getMaxSize(Chromosome[] chromosomes) {
        long maxSize = chromosomes[0].getLength();
        for (Chromosome chromosome : chromosomes) {
            if (chromosome.getLength() > maxSize) {
                maxSize = chromosome.getLength();
            }
        }
        return maxSize;
    }
}
