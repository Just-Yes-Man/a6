/** Uso: java Main p dof */
public class Main {

    public static void main(String[] args) {
        if (args.length != 2) {
            System.out.println("Uso: java Main <p> <dof>");
            return;
        }
        double p = Double.parseDouble(args[0]);
        int dof = Integer.parseInt(args[1]);
        System.out.printf("p=%s dof=%d x=%.5f%n", args[0], dof, SearchX.find(p, dof));
    }
}
