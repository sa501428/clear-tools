package cli.clt.flat.file;

import cli.Main;
import cli.clt.CommandLineParser;
import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.feature2D.Feature2DParser;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.UNIXTools;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

public class GetCommonVsVariableDiffsFromFlatFile {

    public static String usage = "get-common-vs-variable-from-flat-file <genomeID> <flat.file.bedpe> <outfolder> " +
            "<stem1,stem2,...> ...\n" +
            "creates differential loop lists from the flat file\n";
    private static int minNumForVariable = 3;

    public static void run(String[] args, String command, CommandLineParser parser) {
        if (args.length != 5) {
            Main.printGeneralUsageAndExit(58, usage);
        }

        String genomeID = args[1];
        String flatFile = args[2];
        String outFolder = args[3];

        String[] stems = args[4].split(",");
        minNumForVariable = stems.length / 3;

        File folder = UNIXTools.makeDir(new File(outFolder));

        getFromFlatFile(genomeID, flatFile, stems, folder);

        System.out.println("extraction complete");
    }

    private static void getFromFlatFile(String genomeID, String flatFile, String[] stems, File folder) {
        Feature2DList flat = Feature2DParser.loadFeatures(flatFile, ChromosomeTools.loadChromosomes(genomeID),
                true, null, false);

        Feature2DList alwaysPresent = new Feature2DList();
        Feature2DList mostlyPresent = new Feature2DList();
        Feature2DList variablePresent = new Feature2DList();
        Feature2DList neverPresent = new Feature2DList();

        flat.processLists((s, list) -> {

            List<Feature2D> always = new LinkedList<>();
            List<Feature2D> mostly = new LinkedList<>();
            List<Feature2D> variable = new LinkedList<>();
            List<Feature2D> never = new LinkedList<>();


            System.out.println("processing " + s);

            for (Feature2D loop : list) {
                int numPresent = 0;
                int numAbsent = 0;
                int numIndeterminate = 0;
                for (String stem : stems) {
                    int status = Integer.parseInt(loop.getAttribute(stem));
                    if (status == 1) {
                        numPresent++;
                    } else if (status == -1) {
                        numAbsent++;
                    } else {
                        numIndeterminate++;
                    }
                }

                if (numPresent == stems.length) {
                    always.add(loop);
                } else if (numAbsent == stems.length) {
                    never.add(loop);
                } else if (numPresent + numIndeterminate == stems.length) {
                    mostly.add(loop);
                } else if (numPresent >= minNumForVariable && numAbsent >= minNumForVariable) {
                    variable.add(loop);
                }
            }

            synchronized (alwaysPresent) {
                alwaysPresent.addByKey(s, always);
            }
            synchronized (mostlyPresent) {
                mostlyPresent.addByKey(s, mostly);
            }
            synchronized (variablePresent) {
                variablePresent.addByKey(s, variable);
            }
            synchronized (neverPresent) {
                neverPresent.addByKey(s, never);
            }
        });

        alwaysPresent.exportFeatureList(new File(folder, "always.bedpe"), false, Feature2DList.ListFormat.NA);
        mostlyPresent.exportFeatureList(new File(folder, "mostly.bedpe"), false, Feature2DList.ListFormat.NA);
        variablePresent.exportFeatureList(new File(folder, "variable.bedpe"), false, Feature2DList.ListFormat.NA);
        neverPresent.exportFeatureList(new File(folder, "never.bedpe"), false, Feature2DList.ListFormat.NA);
    }
}
