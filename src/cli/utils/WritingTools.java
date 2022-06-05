package cli.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

public class WritingTools {

    public static void writeToMND(int[][] matrix, int resolution,
                                  String xChrom, String yChrom,
                                  int xOrigin, int yOrigin, BufferedWriter bwMND) throws IOException {
        for(int i = 0; i < matrix.length; i++){
            for(int j = 0; j < matrix[i].length; j++){
                if (matrix[i][j] > 0) {
                    long gx = (long) (i + xOrigin) * resolution;
                    long gy = (long) (j + yOrigin) * resolution;
                    bwMND.write(xChrom + " " + gx + " " + yChrom + " " + gy + " " + matrix[i][j]);
                    bwMND.newLine();
                }
            }
        }
    }

    public static String buildCatScript(List<String> filePaths, String outFolder) {
        filePaths.sort(String::compareTo);
        String scriptPath = new File(outFolder, "cat_outputs.sh").getAbsolutePath();
        String resultPath = new File(outFolder, "enhance.short.mnd.txt").getAbsolutePath();

        try {
            PrintWriter finalOutput = new PrintWriter(scriptPath);
            StringBuilder catOutputLine = new StringBuilder();
            StringBuilder removeLine = new StringBuilder();

            catOutputLine.append("cat ").append(filePaths.get(0));
            removeLine.append("rm ").append(filePaths.get(0));
            for(int i = 1; i < filePaths.size(); i++){

                catOutputLine.append(" ").append(filePaths.get(i));
                removeLine.append(" ").append(filePaths.get(i));
            }

            catOutputLine.append(" ").append(" > ").append(resultPath).append("\n");
            finalOutput.println(catOutputLine);
            finalOutput.println(removeLine);
            finalOutput.close();

            System.out.println("Run the shell script at: " + scriptPath);
            //ShellCommandRunner.runShellFile("sh", scriptPath);

        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("Unable to write to catOutputs.sh");
            System.exit(70);
        }

        return resultPath;
    }

    public static String getResolutionsToBuild(int highestResolution) {
        StringBuilder resolutionsToBuild = new StringBuilder("2500000");
        int[] bpBinSizes = {1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 100, 50, 20, 10};
        for (int res : bpBinSizes) {
            if (res >= highestResolution) {
                resolutionsToBuild.append(",").append(res);
            }
        }
        return resolutionsToBuild.toString();
    }

    public static String cleanGenome(String genomeId) {
        if(genomeId.contains("hg38")) return "hg38";
        if(genomeId.contains("hg19")) return "hg19";
        return genomeId;
    }
}
