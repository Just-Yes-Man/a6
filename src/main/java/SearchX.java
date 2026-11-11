/**
 * Busca x tal que la integral de la t de 0 a x sea p.
 *
 * Fuente del algoritmo de búsqueda: no es una fórmula de internet. Viene del
 * enunciado del programa 6 (6a.pdf, "Assignment Kit for Program 6", PSP for
 * Engineers Part II, Carnegie Mellon University / SEI), sección
 * "Search algorithm":
 *   - valor inicial x = 1.0 y d = 0.5
 *   - si la integral queda baja se suma d a x, si queda alta se resta d
 *   - cuando cambia el signo del error se divide d entre 2
 *   - se repite hasta que el error sea aceptable
 *
 * Las fórmulas matemáticas (densidad t, Gamma y Simpson) están en
 * TDistribution.java y Gamma.java, cada una con su fuente.
 */
public class SearchX {

    // Error aceptable. El PDF sugiere 0.00001, pero con ese valor el caso
    // p = 0.495, dof = 4 daba x = 4.60156 en vez de 4.60409, por lo que se
    // usa uno más estricto (decisión propia).
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
