package cli.clt;

import cli.utils.seer.CumulativeDistributionFunction;
import cli.utils.seer.SeerUtils;
import cli.utils.sift.SimpleLocation;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import javastraw.tools.UNIXTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class Seer {
    /* takes in one file currently (for ease of testing: can change later to a list of files and easily iterate over).
    resolution also input but currently set to 50 manually.
    */

    public static void generateNewReads(String filename, int lowResolution, int highResolution,
                                        String outFolderPath, String possibleNorm, long seed, long numberOfContacts) {
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, false);
        NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);

        HiCZoom lowestResolution = getLowestResolution(ds.getBpZooms());
        Map<Chromosome, Long> contactsPerChromosome = generateCountsForEachChromosome(ds, lowestResolution);
        Map<Chromosome, Long> countsToGeneratePerChr = generateCountsToMake(numberOfContacts, contactsPerChromosome);

        contactsPerChromosome.clear();
        Random rand = new Random(seed);
        Random rand2 = new Random(seed * 3 + 5);

        AtomicInteger index0 = new AtomicInteger(0);
        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        ParallelizationTools.launchParallelizedCode(() -> {
            int indexForJob = index0.getAndIncrement();
            while (indexForJob < chromosomes.length) {
                Chromosome chromosome = chromosomes[indexForJob];
                Matrix matrix = ds.getMatrix(chromosome, chromosome);
                if (matrix != null) {
                    MatrixZoomData zdHigh = matrix.getZoomData(new HiCZoom(highResolution));
                    if (zdHigh != null) {
                        double[] rowSums = SeerUtils.convertToCDF(SeerUtils.getRowSumsForZD(chromosome, highResolution,
                                zdHigh.getDirectIterator()));
                        matrix.clearCacheForZoom(new HiCZoom(highResolution));

                        MatrixZoomData zdLow = matrix.getZoomData(new HiCZoom(lowResolution));
                        if (zdLow != null) {

                            CumulativeDistributionFunction cdf = new CumulativeDistributionFunction(zdLow.
                                    getNormalizedIterator(norm), 10000000, lowResolution);
                            String name = chromosome.getName();
                            matrix.clearCacheForZoom(new HiCZoom(lowResolution));

                            try {
                                File outputFileName = new File(outFolderPath, name + ".generated.contacts.mnd.txt");
                                outputFileName.createNewFile();
                                FileWriter fw = new FileWriter(outputFileName);
                                BufferedWriter bw = new BufferedWriter(fw);

                                // generate points at random
                                long numPointsToGenerate = countsToGeneratePerChr.get(chromosome);
                                for (long i = 0; i < numPointsToGenerate; i++) {
                                    SimpleLocation position = cdf.createRandomPoint(rand);
                                    position = SeerUtils.updateToHigherResPosition(position, rowSums, lowResolution,
                                            highResolution, rand2);
                                    bw.write(name + " " + position.getBinX() + " " + name + " " + position.getBinY());
                                    bw.newLine();
                                }
                                bw.close();
                            } catch (Exception e) {
                                e.printStackTrace();
                                System.exit(9);
                            }
                        }
                    }
                }
                indexForJob = index0.getAndIncrement();
            }
        });

        System.out.println("Seer completed");
    }

    private static Map<Chromosome, Long> generateCountsForEachChromosome(Dataset ds, HiCZoom lowestResolution) {
        Map<Chromosome, Long> results = new HashMap<>();
        for (Chromosome chromosome : ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(lowestResolution);
            if (zd == null) continue;

            long total = 0;
            for (Iterator<ContactRecord> it = zd.getDirectIterator(); it.hasNext(); ) {
                ContactRecord record = it.next();
                if (record.getCounts() > 0) {
                    total += record.getCounts();
                }
            }
            matrix.clearCacheForZoom(lowestResolution);
            results.put(chromosome, total);
        }
        return results;
    }

    private static HiCZoom getLowestResolution(List<HiCZoom> zooms) {
        HiCZoom zoom = zooms.get(0);
        for (HiCZoom z : zooms) {
            if (z.getBinSize() > zoom.getBinSize()) {
                zoom = z;
            }
        }
        return zoom;
    }

    private static Map<Chromosome, Long> generateCountsToMake(double totalNumberOfContacts, Map<Chromosome, Long> countMap) {
        double total = 0;
        for (Long value : countMap.values()) {
            total += value;
        }
        Map<Chromosome, Long> countsToMake = new HashMap<>();
        for (Chromosome chromosome : countMap.keySet()) {
            Long numToMake = Math.round(totalNumberOfContacts * (countMap.get(chromosome) / total));
            countsToMake.put(chromosome, numToMake);
        }
        return countsToMake;
    }

    public static void run(String[] args, CommandLineParser parser) {
        // check length of arguments equal to 4
        if (args.length != 4) {
            printUsageAndExit();
        }
        int highResolution = parser.getResolutionOption(50);
        int lowResolution = parser.getLowResolutionOption(5000);
        String possibleNorm = parser.getNormalizationStringOption();
        long seed = parser.getSeedOption(0);
        UNIXTools.makeDir(args[2]);
        long numContactsToGenerate = Long.parseLong(args[3]);
        generateNewReads(args[1], lowResolution, highResolution, args[2], possibleNorm, seed, numContactsToGenerate);
    }

    private static void printUsageAndExit() {
        System.out.println("seer [-r 50 (high res)] [--low-res 5000 (low res)] [-k SCALE (normalization)] " +
                "[--seed 0 (seed for random number generator)] " +
                "<input.hic> <output_folder> <number of contacts to generate>");
        System.exit(19);
    }
}