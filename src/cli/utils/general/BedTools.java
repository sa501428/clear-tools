package cli.utils.general;

import cli.utils.flags.Anchor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Set;

public class BedTools {
    public static void exportBedFile(File file, Set<Anchor> anchors) {
        try (PrintWriter writer = new PrintWriter(file.getPath())) {
            for (Anchor anchor : anchors) {
                writer.println(anchor.toString());
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error exporting bed file (" + file.getPath() + "): " + e.getMessage());
        }
    }
}
