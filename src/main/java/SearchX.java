public class SearchX {

    public static final double ERROR = 1e-9;
    private static final double START_X = 1.0;
    private static final double START_D = 0.5;
    static int maxIterations = 10_000;

    public static double find(double p, int dof) {
        if (p <= 0.0 || p >= 0.5) {
            throw new IllegalArgumentException("p debe estar entre 0 y 0.5 (exclusivo)");
        }
        TDistribution t = new TDistribution(dof);
        double x = START_X;
        double d = START_D;
        double error = t.integrate(x) - p;
        int lastSign = Double.compare(error, 0.0);
        int iterations = 0;

        while (Math.abs(error) > ERROR && iterations++ < maxIterations) {
            x = (error < 0) ? x + d : Math.max(0.0, x - d);
            error = t.integrate(x) - p;
            int sign = Double.compare(error, 0.0);
            if (sign != lastSign) {
                d = d / 2.0;
                lastSign = sign;
            }
        }
        return x;
    }
}
