
public class Gamma {

      private static final double[] COF = {
        76.18009172947146, -86.50532032941677, 24.01409824083091,
        -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
    };

       public static double logGamma(double x) {
        double y = x;
        double tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        double ser = 1.000000000190015;
        for (int j = 0; j < COF.length; j++) {
            ser += COF[j] / ++y;
        }
             return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    public static double of(double x) {
        return Math.exp(logGamma(x));
    }
}
