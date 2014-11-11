package se.dxtr;

import java.math.BigDecimal;
import java.math.BigInteger;

/**
 * Created by dexter on 11/11/14.
 */
public class QS {
    public static double getB(BigInteger n) {
        double c = 3;
        double tempN = n.doubleValue();
        double logs = Math.log(tempN) * Math.log(Math.log(tempN));
        double b = c * Math.exp(0.5 * Math.sqrt(logs));
        return Math.round(b);
    }

    public static BigInteger Q(BigInteger x, BigInteger n) {
        double parenthesis = Math.floor(Math.sqrt(n.doubleValue())) + x.doubleValue();
        double q = Math.pow(parenthesis, 2) - n.doubleValue();
        return new BigDecimal(q).toBigInteger();
    }


}
