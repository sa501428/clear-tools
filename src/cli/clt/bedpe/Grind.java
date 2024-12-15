package cli.clt.bedpe;

import cli.clt.CommandLineParser;
import cli.utils.apa.MultiAPAManager;
import cli.utils.data.SparseContactMatrixWithMasking;
import io.jhdf.HdfFile;
import io.jhdf.WritableHdfFile;
import io.jhdf.api.WritableGroup;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class Grind {

    public static String usage = "grind [-k NORM] [-r resolution] [--window half-width] [--npy] <hic file> <bedpe> <directory> <stem>\n" +
            "\t\tsplit up the bedpe into multiple lists; number < 2 splits by chromosome";
    private final boolean useObservedOverExpected = false;
    private final boolean useNpy;
    private final Dataset ds;
    private final File outputDirectory;
    private final HiCZoom zoom;
    private final Feature2DList loopList;
    private final ChromosomeHandler handler;
    private int matrixHalfWidth = 10;
    private Integer resolution = 1000;
    private NormalizationType norm;
    private String stemName;


    // Lists to store loop information
    private final List<String> chr1List = new ArrayList<>();
    private final List<Long> start1List = new ArrayList<>();
    private final List<Long> end1List = new ArrayList<>();
    private final List<String> chr2List = new ArrayList<>();
    private final List<Long> start2List = new ArrayList<>();
    private final List<Long> end2List = new ArrayList<>();
    private final List<Integer> chunkIndexList = new ArrayList<>();


    // Lock object for synchronizing access to loop info lists
    private final Object loopInfoLock = new Object();

    public Grind(String[] args, CommandLineParser parser) {
        if (args.length != 5) {
            printUsageAndExit();
        }

        resolution = parser.getResolutionOption(5000);
        boolean useBI = resolution >= 50;

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, useBI);
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        stemName = sanitizeStemName(args[4]);

        String possibleNorm = parser.getNormalizationStringOption();

        try {
            norm = ds.getNormalizationHandler().getNormTypeFromString(possibleNorm);
        } catch (Exception e) {
            norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{possibleNorm, "SCALE", "KR", "NONE"});
        }

        System.out.println("Using normalization: " + norm.getLabel());


        matrixHalfWidth = parser.getWindowSizeOption(10);

        useNpy = parser.getNpyOption();

        zoom = new HiCZoom(resolution);
        handler = ds.getChromosomeHandler();


        loopList = Feature2DParser.loadFeatures(args[2], handler, false, null, false);
        if (loopList.getNumTotalFeatures() < 1) {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }
        System.out.println("Using stem name: " + stemName);
    }


    public void run() {
        if (useNpy){
            System.out.println("Using Numpy file type");
             LoopDumper.dump(ds, loopList, outputDirectory, handler, norm,
                    useObservedOverExpected, resolution, matrixHalfWidth);
        }
        else{
        System.out.println("Using HDF5 file type");
        String resolutionGroupName = String.valueOf(resolution);

        File hdf5File = new File(outputDirectory, stemName+"_loops_at_"+resolutionGroupName+"_bp.hdf5");

        try (WritableHdfFile writableHdfFile = HdfFile.write(Paths.get(hdf5File.getAbsolutePath()))) {

            WritableGroup resolutionGroup = writableHdfFile.putGroup(resolutionGroupName);

            WritableGroup chunksGroup = resolutionGroup.putGroup("chunks");
            WritableGroup loopInfoGroup = resolutionGroup.putGroup("loop_info");

            AtomicInteger chunkIndex = new AtomicInteger(0);

            List<float[][]> currentChunk = new ArrayList<>(100);

            loopList.processLists((chr, feature2DList) -> {

                System.out.println("Currently on: " + chr);

                Chromosome chrom = handler.getChromosomeFromName(feature2DList.get(0).getChr1());

                Matrix matrix = ds.getMatrix(chrom, chrom);
                if (matrix == null) return;

                HiCZoom zoom = ds.getZoomForBPResolution(resolution);
                final MatrixZoomData zd = matrix.getZoomData(zoom);

                if (zd == null) return;

                int matrixWidth = 2 * matrixHalfWidth + 1;

                try {
                    SparseContactMatrixWithMasking scm = new SparseContactMatrixWithMasking(zd,
                            feature2DList, resolution, matrixHalfWidth, matrixWidth, norm, true); // intra here

                    AtomicInteger index = new AtomicInteger(0);

                    ParallelizationTools.launchParallelizedCode(10, () -> {
                        int currIndex = index.getAndIncrement();

                        while (currIndex < feature2DList.size()) {

                            Feature2D loop = feature2DList.get(currIndex);
                            float[][] output = new float[matrixWidth][matrixWidth];
                            MultiAPAManager.addToMatrix(output, scm, loop, matrixHalfWidth, resolution, matrixWidth);

                            String serializedLoop = loop.toString();

                            String chr1 = loop.getChr1();
                            long start1 = loop.getStart1();
                            long end1 = loop.getEnd1();
                            String chr2 = loop.getChr2();
                            long start2 = loop.getStart2();
                            long end2 = loop.getEnd2();

                            int assignedChunkIndex = -1;

                            synchronized (loopInfoLock) {
                                currentChunk.add(output);
                                if (currentChunk.size() == 100) {
                                    float[][][] chunkArray = new float[100][matrixWidth][matrixWidth];
                                    for (int i = 0; i < 100; i++) {
                                        chunkArray[i] = currentChunk.get(i);
                                    }

                                    assignedChunkIndex = chunkIndex.getAndIncrement();

                                    String datasetName = "chunk_" + assignedChunkIndex;
                                    chunksGroup.putDataset(datasetName, chunkArray);
                                    currentChunk.clear();
                                }
                            }

                            synchronized (loopInfoLock) {
                                chr1List.add(chr1);
                                start1List.add(start1);
                                end1List.add(end1);
                                chr2List.add(chr2);
                                start2List.add(start2);
                                end2List.add(end2);
                                assignedChunkIndex = chunkIndex.get() ; // Last assigned chunk
                                chunkIndexList.add(assignedChunkIndex);
                            }

                            currIndex = index.getAndIncrement();
                        }
                    });

                } catch (Exception ex) {
                    System.err.println("Error processing: " + chr);
                    ex.printStackTrace();
                }

            });

            synchronized (loopInfoLock) {
                if (!currentChunk.isEmpty()) {
                    int remaining = currentChunk.size();
                    float[][][] chunkArray = new float[remaining][2 * matrixHalfWidth + 1][2 * matrixHalfWidth + 1];
                    for (int i = 0; i < remaining; i++) {
                        chunkArray[i] = currentChunk.get(i);
                    }
                    String datasetName = "chunk_" + chunkIndex.getAndIncrement();
                    chunksGroup.putDataset(datasetName, chunkArray);
                    currentChunk.clear();
                }
            }

            loopInfoGroup.putDataset("chr1", chr1List.toArray(new String[0]));
            loopInfoGroup.putDataset("start1", start1List.stream().mapToLong(Long::intValue).toArray());
            loopInfoGroup.putDataset("end1", end1List.stream().mapToLong(Long::intValue).toArray());
            loopInfoGroup.putDataset("chr2", chr2List.toArray(new String[0]));
            loopInfoGroup.putDataset("start2", start2List.stream().mapToLong(Long::intValue).toArray());
            loopInfoGroup.putDataset("end2", end2List.stream().mapToLong(Long::intValue).toArray());
            loopInfoGroup.putDataset("chunk_index", chunkIndexList.stream().mapToInt(Integer::intValue).toArray());

            System.out.println("All loops have been successfully written to " + hdf5File.getAbsolutePath());

        } catch (Exception e) {
            System.err.println("Error writing to HDF5 file.");
            e.printStackTrace();
            System.exit(4);
        }
    }
}

    private String sanitizeStemName(String stemName) {
        return stemName.replaceAll("[\\\\/:*?\"<>|]", "_");
    }

    private void printUsageAndExit() {
        System.out.println(usage);
        System.exit(19);
    }
}
