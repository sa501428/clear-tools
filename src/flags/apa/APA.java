/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package flags.apa;

import javastraw.StrawGlobals;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.feature2D.Feature2D;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class APA {

    private final Dataset ds;
    private final NormalizationType norm;
    private final File outputDirectory;
    private final int window = 20;
    private final int resolution = 10000;
    private final Object key = new Object();
    private final GenomeWide1DList<Anchor> anchors;


    public APA(Dataset ds, String outfolder, NormalizationType norm, GenomeWide1DList<Anchor> anchors) {
        this.ds = ds;
        this.norm = norm;
        this.outputDirectory = HiCFileTools.createValidDirectory(outfolder);
        this.anchors = anchors;
    }

    public void run() {
        int matrixWidth = 2 * window + 1;

        System.out.println("Processing APA for resolution " + resolution);
        HiCZoom zoom = new HiCZoom(HiCZoom.HiCUnit.BP, resolution);

        ChromosomeHandler handler = ds.getChromosomeHandler();

        if (anchors.size() < 2) {
            System.err.println("Loop list is empty or incorrect path provided.");
            System.exit(3);
        }


        APADataStack interDataStack = new APADataStack(matrixWidth, outputDirectory, "inter_");
        APADataStack[] intraDataStacks = DataStackUtils.initialize(handler.getChromosomeArrayWithoutAllByAll(),
                matrixWidth, outputDirectory, resolution);

        double maxProgressStatus = handler.size();
        final AtomicInteger currentProgressStatus = new AtomicInteger(0);

        Map<Integer, Chromosome[]> chromosomePairs = new ConcurrentHashMap<>();
        int pairCounter = 0;
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                Chromosome[] chromosomePair = {chromosomes[i], chromosomes[j]};
                chromosomePairs.put(pairCounter, chromosomePair);
                pairCounter++;
            }
        }
        final int chromosomePairCounter = pairCounter;
        final AtomicInteger chromosomePair = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(() -> {

            int threadPair = chromosomePair.getAndIncrement();
            while (threadPair < chromosomePairCounter) {
                Chromosome chr1 = chromosomePairs.get(threadPair)[0];
                Chromosome chr2 = chromosomePairs.get(threadPair)[1];

                MatrixZoomData zd;
                synchronized (key) {
                    zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, zoom);
                }

                if (zd == null) {
                    threadPair = chromosomePair.getAndIncrement();
                    continue;
                }

                if (StrawGlobals.printVerboseComments) {
                    System.out.println("CHR " + chr1.getName() + " " + chr1.getIndex() + " CHR " + chr2.getName() + " " + chr2.getIndex());
                }

                for (int distBin = 0; distBin < intraDataStacks.length; distBin++) {
                    // inter only done once
                    if (chr1.getIndex() != chr2.getIndex() && distBin > 0) continue;

                    APADataStack accumDataStack;
                    if (chr1.getIndex() == chr2.getIndex()) {
                        accumDataStack = intraDataStacks[distBin];
                    } else {
                        accumDataStack = interDataStack;
                    }

                    long minDist = (long) (Math.pow(2, distBin - 1) * 1000000L);
                    long maxDist = (long) (Math.pow(2, distBin) * 1000000L);
                    if (distBin == 0) {
                        minDist = 200000;
                    }

                    List<Feature2D> loops = LoopGenerator.generate(anchors, chr1, chr2, minDist, maxDist);
                    if (loops.size() < 1) {
                        if (StrawGlobals.printVerboseComments) {
                            System.out.println("CHR " + chr1.getName() + " CHR " + chr2.getName() + " - no loops, check loop filtering constraints");
                        }
                        threadPair = chromosomePair.getAndIncrement();
                        continue;
                    }

                    APADataStack temp = new APADataStack(matrixWidth, outputDirectory, "temp");
                    for (Feature2D loop : loops) {
                        try {
                            RealMatrix newData;
                            synchronized (key) {
                                newData = APAUtils.extractLocalizedData(zd, loop, matrixWidth, resolution, window, norm);
                            }
                            temp.addData(newData);
                        } catch (Exception e) {
                            System.err.println(e.getMessage());
                            System.err.println("Unable to find data for loop: " + loop);
                        }
                    }
                    synchronized (key) {
                        accumDataStack.addData(temp.getData());
                    }
                }

                if (chr2.getIndex() == chr1.getIndex()) {
                    System.out.print(((int) Math.floor((100.0 * currentProgressStatus.incrementAndGet()) / maxProgressStatus)) + "% ");
                }

                threadPair = chromosomePair.getAndIncrement();
            }
        });

        System.out.println("Exporting APA results...");
        for (APADataStack dataStack : intraDataStacks) {
            dataStack.exportData();
        }
        interDataStack.exportData();
        System.out.println("APA complete");
        //if no data return null
    }

}