package cli.utils.general;

import cli.utils.flags.Anchor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class BedTools {
    public static void exportBedFile(File file, Set<Anchor> anchorSet) {
        List<Anchor> anchors = new ArrayList<>(anchorSet);
        anchors.sort((o1, o2) -> {
            if (o1.getKey().equals(o2.getKey())) {
                return Long.compare(o1.getStart(), o2.getStart());
            }
            return o1.getKey().compareTo(o2.getKey());
        });

        try (PrintWriter writer = new PrintWriter(file.getPath())) {
            for (Anchor anchor : anchors) {
                writer.println(anchor.toString());
            }
        } catch (FileNotFoundException e) {
            System.err.println("Error exporting bed file (" + file.getPath() + "): " + e.getMessage());
        }
    }
}
