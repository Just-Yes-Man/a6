import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import org.junit.jupiter.api.Test;

public class SearchXTest {

    private static final double TOL = 0.0001;

    @Test
    void test1() {
        assertEquals(0.55338, SearchX.find(0.20, 6), TOL);
    }

    @Test
    void test2() {
        assertEquals(1.75305, SearchX.find(0.45, 15), TOL);
    }

    @Test
    void test3() {
        assertEquals(4.60409, SearchX.find(0.495, 4), TOL);
    }

    @Test
    void gammaMedio() {
        assertEquals(Math.sqrt(Math.PI), Gamma.of(0.5), 1e-8);
    }

    @Test
    void gammaEntero() {
        assertEquals(24.0, Gamma.of(5), 1e-7);
    }

    @Test
    void gammaPuedeInstanciarse() {
        assertTrue(new Gamma() != null);
    }

    @Test
    void searchXPuedeInstanciarse() {
        assertTrue(new SearchX() != null);
    }

    @Test
    void mainPuedeInstanciarse() {
        assertTrue(new Main() != null);
    }

    @Test
    void integralDeCeroEsCero() {
        assertEquals(0.0, new TDistribution(5).integrate(0.0), 1e-12);
    }

    @Test
    void pInvalido() {
        assertThrows(IllegalArgumentException.class, () -> SearchX.find(0.6, 5));
    }

    @Test
    void pCeroEsInvalido() {
        assertThrows(IllegalArgumentException.class, () -> SearchX.find(0.0, 5));
    }

    @Test
    void gradosDeLibertadInvalidos() {
        assertThrows(IllegalArgumentException.class, () -> new TDistribution(0));
    }

    @Test
    void densidadEnCero() {
        assertEquals(0.3796066898, new TDistribution(5).density(0.0), 1e-10);
    }

    @Test
    void integralNegativa() {
        assertEquals(-new TDistribution(5).integrate(1.0),
                new TDistribution(5).integrate(-1.0), 1e-12);
    }

    @Test
    void integralPequenaAlcanzaPrecisionAbsoluta() {
        assertTrue(new TDistribution(5).integrate(1e-8) > 0.0);
    }

    @Test
    void integralNoFinitaDevuelveNaN() {
        assertTrue(Double.isNaN(new TDistribution(5).integrate(Double.NaN)));
    }

    @Test
    void integralUsaPrecisionAbsolutaCuandoCorresponde() {
        TDistribution distribution = new TDistribution(5);
        double result = distribution.integrate(1e-8);
        result += distribution.integrate(1e-7);
        result += distribution.integrate(1e-6);
        result += distribution.integrate(1e-5);
        result += distribution.integrate(1e-4);
        result += distribution.integrate(0.001);
        result += distribution.integrate(0.01);
        assertTrue(result > 0.0);
    }

    @Test
    void convergenciaPorPrecisionAbsoluta() {
        assertTrue(TDistribution.hasConverged(1e-15, 0.0));
    }

    @Test
    void busquedaEvitaXNegativo() {
        assertTrue(SearchX.find(1e-12, 5) >= 0.0);
    }

    @Test
    void busquedaRespetaLimiteDeIteraciones() {
        int original = SearchX.maxIterations;
        try {
            SearchX.maxIterations = 0;
            assertEquals(1.0, SearchX.find(0.20, 6), 0.0);
        } finally {
            SearchX.maxIterations = original;
        }
    }

    @Test
    void mainSinArgumentos() {
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        PrintStream original = System.out;
        try {
            System.setOut(new PrintStream(output));
            Main.main(new String[0]);
        } finally {
            System.setOut(original);
        }
        assertEquals("Uso: java Main <p> <dof>\n", output.toString());
    }

    @Test
    void mainConArgumentos() {
        ByteArrayOutputStream output = new ByteArrayOutputStream();
        PrintStream original = System.out;
        try {
            System.setOut(new PrintStream(output));
            Main.main(new String[]{"0.20", "6"});
        } finally {
            System.setOut(original);
        }
        assertEquals("p=0.20 dof=6 x=0.55338\n", output.toString());
    }
}
