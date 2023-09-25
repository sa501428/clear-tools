package cli.clt.bedpe;

import cli.utils.apa.MultiAPAManager;
import cli.utils.data.SparseContactMatrixWithMasking;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;

import java.io.File;
import java.util.concurrent.atomic.AtomicInteger;

public class LoopDumper {
    public static void dump(Dataset ds, Feature2DList loopList, File outputDirectory,
                            ChromosomeHandler handler, NormalizationType norm,
                            boolean useObservedOverExpected, int resolution,
                            int window) {


        if (!outputDirectory.isDirectory()) {
            outputDirectory.mkdir();
        }

        File folder = new File(outputDirectory, String.valueOf(resolution));
        if (!folder.exists()) {
            folder.mkdir();
        }

        try {

            loopList.processLists((chr, feature2DList) -> {

                System.out.println("Currently on: " + chr);

                Chromosome chrom = handler.getChromosomeFromName(feature2DList.get(0).getChr1());

                Matrix matrix = ds.getMatrix(chrom, chrom);
                if (matrix == null) return;

                HiCZoom zoom = ds.getZoomForBPResolution(resolution);
                final MatrixZoomData zd = matrix.getZoomData(zoom);

                if (zd == null) return;

                int matrixWidth = 2 * window + 1;

                SparseContactMatrixWithMasking scm = new SparseContactMatrixWithMasking(zd,
                        feature2DList, resolution, window, matrixWidth, norm);


                AtomicInteger index = new AtomicInteger(0);

                ParallelizationTools.launchParallelizedCode(10, new Runnable() {
                    @Override
                    public void run() {
                        int currIndex = index.getAndIncrement();

                        while (currIndex < feature2DList.size()) {

                            Feature2D loop = feature2DList.get(currIndex);
                            float[][] output = new float[matrixWidth][matrixWidth];
                            MultiAPAManager.addToMatrix(output, scm, loop, window, resolution, matrixWidth);
                            String name = loop.getChr1() + "_" + loop.getStart1() + "_" + loop.getEnd1() + "_" +
                                    loop.getChr2() + "_" + loop.getStart2() + "_" + loop.getEnd2();
                            String path = new File(folder, name + ".npy").getAbsolutePath();
                            MatrixTools.saveMatrixTextNumpy(path, output);
                            currIndex = index.getAndIncrement();
                        }
                    }
                });

            });
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
