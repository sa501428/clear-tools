package flags;

public class Main {

    public static void main(String[] args) {
        // 0 filename
        // 1 bed file
        AVD avd = new AVD(args[0], args[1], Integer.parseInt(args[2]));
    }
}
