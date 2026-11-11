/**
 * Densidad de la distribución t de Student e integral por la regla de Simpson.
 *
 * CÓDIGO TOMADO / ADAPTADO DE PROGRAMAS EXISTENTES:
 *
 * 1) Densidad t (constructor y density): adaptado de la clase TDistribution
 *    de Apache Commons Math 3.6.1 (Apache License 2.0), métodos del
 *    constructor (cálculo de "factor") y logDensity/density:
 *    https://code.yawk.at/org.apache.commons/commons-math3/3.6.1/org/apache/commons/math3/distribution/TDistribution.java
 *    Fórmulas de esa clase (n = grados de libertad):
 *        factor     = logGamma((n + 1) / 2) - 0.5 * (log(pi) + log(n)) - logGamma(n / 2)
 *        logDensity = factor - ((n + 1) / 2) * log(1 + x * x / n)
 *        density    = exp(logDensity)
 *    Cambio hecho: FastMath se sustituyó por java.lang.Math.
 *    La forma equivalente de la fórmula también aparece en
 *    https://en.wikipedia.org/wiki/Student%27s_t-distribution
 *
 * 2) Simpson (integrate): adaptado de SimpsonIntegrator.doIntegrate de
 *    Apache Commons Math / Hipparchus (Apache License 2.0):
 *    https://hipparchus.org/hipparchus-core/jacoco/org.hipparchus.analysis.integration/SimpsonIntegrator.java.html
 *    De ahí se tomó:
 *      - Simpson como combinación de dos trapecios consecutivos:
 *            s = (4 * t - oldt) / 3.0
 *      - prueba de convergencia con delta = |s - olds|:
 *            rLimit = relativeAccuracy * (|olds| + |s|) * 0.5
 *            termina si delta <= rLimit  o  delta <= absoluteAccuracy
 *      - número mínimo de iteraciones antes de probar la convergencia.
 *    Los valores por defecto de la librería son relativeAccuracy = 1e-6,
 *    absoluteAccuracy = 1e-15 y mínimo de 3 iteraciones; aquí se usa
 *    relativeAccuracy = 1e-12 para obtener los 5 decimales que pide la tabla.
 *
 * 3) trapezoid (refinamiento de trapecios): es el esquema estándar que usa
 *    TrapezoidIntegrator de la misma librería: la etapa 0 es
 *    (b - a) / 2 * (f(a) + f(b)) y cada etapa n suma los puntos medios de
 *    los subintervalos nuevos y promedia con la etapa anterior.
 *    Esta parte se escribió siguiendo ese método, no es copia literal.
 */
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
        // Commons Math TDistribution: factor del constructor
        double nPlus1Over2 = (n + 1) / 2;
        this.factor = Gamma.logGamma(nPlus1Over2)
                - 0.5 * (Math.log(Math.PI) + Math.log(n))
                - Gamma.logGamma(n / 2);
    }

    /** Commons Math TDistribution.logDensity. */
    private double logDensity(double x) {
        double nPlus1Over2 = (n + 1) / 2;
        return factor - nPlus1Over2 * Math.log(1 + x * x / n);
    }

    /** Commons Math TDistribution.density: exp(logDensity). Siempre mayor que 0. */
    public double density(double x) {
        return Math.exp(logDensity(x));
    }

    /** Integral de 0 a x (adaptado de SimpsonIntegrator.doIntegrate). */
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

    /** Etapa n de la regla del trapecio (esquema de TrapezoidIntegrator). */
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
