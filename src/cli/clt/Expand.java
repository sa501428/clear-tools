package cli.clt;

import cli.Main;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class Expand {
    public static String usage = "expand[-exact][-clean] <genomeID> <anchor_size> <input.bedpe> <output.bedpe>\n" +
            "\t\tdefault behavior expands the size of anchors for each loop to be anchor_size\n" +
            "unless the loop's width is already bigger than anchor_size; expansion is in both \n" +
            "directions, so the midpoint is preserved\n" +
            "\t\texact will explicitly set the size to anchor_size even if the anchor was bigger\n" +
            "\t\tclean avoids saving old attributes";

    public static void run(String[] args, String command) {
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(5, usage);
        }
        String genomeID = args[1];
        int newSize = Integer.parseInt(args[2]);
        String inFile = args[3];
        String outFile = args[4];
        setAnchorSize(inFile, genomeID, outFile, newSize,
                command.contains("clean"), command.contains("exact"));
        System.out.println("anchor expansion complete");
    }

    private static void setAnchorSize(String inputBedpe, String genomeID, String outFile,
                                      int newAnchorSize, boolean noAttributes, boolean setExactSize) {
        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(genomeID);
        Feature2DList loopList = Feature2DParser.loadFeatures(inputBedpe, handler, !noAttributes,
                null, false);

        loopList.filterLists((s, list) -> {
            List<Feature2D> list2 = new ArrayList<>(list.size());
            for (Feature2D feature : list) {
                list2.add(getUpdatedLoop(feature, newAnchorSize, setExactSize));
            }
            return list2;
        });

        loopList.exportFeatureList(new File(outFile), false, Feature2DList.ListFormat.NA);
    }

    private static Feature2D getUpdatedLoop(Feature2D feature, int newAnchorSize, boolean setExactSize) {
        long mid1 = feature.getMidPt1();
        long mid2 = feature.getMidPt2();
        long start1 = mid1 - (newAnchorSize / 2);
        long start2 = mid2 - (newAnchorSize / 2);
        long end1 = start1 + newAnchorSize;
        long end2 = start2 + newAnchorSize;

        if (!setExactSize) {
            if (feature.getWidth1() > newAnchorSize) {
                start1 = feature.getStart1();
                end1 = feature.getEnd1();
            }
            if (feature.getWidth2() > newAnchorSize) {
                start2 = feature.getStart2();
                end2 = feature.getEnd2();
            }
        }
        return new Feature2D(feature.getFeatureType(),
                feature.getChr1(), start1, end1,
                feature.getChr2(), start2, end2,
                feature.getColor(), feature.getAttributes());
    }
}
