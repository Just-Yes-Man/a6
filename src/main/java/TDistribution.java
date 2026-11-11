public class TDistribution {

    private static final double RELATIVE_ACCURACY = 1e-12;
    private static final double ABSOLUTE_ACCURACY = 1e-15;
    private static final int MIN_ITERATIONS = 3;
    private static final int MAX_ITERATIONS = 40;

    private final double n;
    private final double factor;

    public TDistribution(int dof) {
        if (dof < 1) {
            throw new IllegalArgumentException("dof debe ser >= 1");
        }
        this.n = dof;
        double nPlus1Over2 = (n + 1) / 2;
        this.factor = Gamma.logGamma(nPlus1Over2)
                - 0.5 * (Math.log(Math.PI) + Math.log(n))
                - Gamma.logGamma(n / 2);
    }

    private double logDensity(double x) {
        double nPlus1Over2 = (n + 1) / 2;
        return factor - nPlus1Over2 * Math.log(1 + x * x / n);
    }

    public double density(double x) {
        return Math.exp(logDensity(x));
    }

    public double integrate(double x) {
        if (x == 0.0) {
            return 0.0;
        }
        double oldt = trapezoid(0.0, x, 0, 0.0);
        double olds = 0.0;
        double t = oldt;
        double s = 0.0;
        for (int iterations = 1; iterations <= MAX_ITERATIONS; iterations++) {
            t = trapezoid(0.0, x, iterations, oldt);
            s = (4 * t - oldt) / 3.0;
            if (iterations >= MIN_ITERATIONS) {
                double delta = Math.abs(s - olds);
                double rLimit = RELATIVE_ACCURACY * (Math.abs(olds) + Math.abs(s)) * 0.5;
                if (hasConverged(delta, rLimit)) {
                    return s;
                }
            }
            olds = s;
            oldt = t;
        }
        return s;
    }

    static boolean hasConverged(double delta, double rLimit) {
        return delta <= rLimit || delta <= ABSOLUTE_ACCURACY;
    }

    private double trapezoid(double min, double max, int stage, double previous) {
        if (!Double.isFinite(max)) {
            return Double.NaN;
        }
        if (stage == 0) {
            return 0.5 * (max - min) * (density(min) + density(max));
        }
        long np = 1L << (stage - 1);
        double sum = 0.0;
        double spacing = (max - min) / np;
        double x = min + 0.5 * spacing;
        for (long i = 0; i < np; i++) {
            sum += density(x);
            x += spacing;
        }
        return 0.5 * (previous + sum * spacing);
    }
}
