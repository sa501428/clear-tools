package cli.clt;

public class SplitOrJoin {
    public SplitOrJoin(String name, String[] args) {

        String genomeID = args[1];
        if (name.contains("split")) {
            // split <genomeID> <number> <input.bedpe>
            int n = Integer.parseInt(args[2]);


        } else {
            // join <genomeID> <output.bedpe> <input1.bedpe> <input2.bedpe> ...
            String outFile = args[2];
            String[] bedpeFiles = new String[args.length - 3];
            System.arraycopy(args, 3, bedpeFiles, 0, bedpeFiles.length);

        }


        System.out.println(name + " complete");
    }
}
