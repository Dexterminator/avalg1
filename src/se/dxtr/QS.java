package se.dxtr;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

/**
 * Created by dexter on 11/11/14.
 */
public class QS {
    public static int getB(BigInteger n) {
        double c = 3;
        double tempN = n.doubleValue();
        double logs = Math.log(tempN) * Math.log(Math.log(tempN));
        double b = c * Math.exp(0.5 * Math.sqrt(logs));
        return (int) Math.round(b);
    }

    public static BigInteger Q(BigInteger x, BigInteger n) {
        double parenthesis = Math.floor(Math.sqrt(n.doubleValue())) + x.doubleValue();
        double q = Math.pow(parenthesis, 2) - n.doubleValue();
        return new BigDecimal(q).toBigInteger();
    }

    public static int legendre(BigInteger N, int p) {
        int result = 1;
        int a = N.mod(BigInteger.valueOf(p)).intValue();
        int power = (p - 1) / 2;

        while (power > 0) {
            if (power % 2 == 1) {
                result = (result * a) % p;
            }
            a = (a * a) % p;
            power = power / 2;
        }
        if (result - p == -1) {
            result = result - p;
        }
        return result;
    }

    public static ArrayList<Integer> factorBase(BigInteger n, int b) {
        ArrayList<Integer> factorBase = new ArrayList<Integer>();
        for (int i = 2; i < b; i++) {
            if (BigInteger.valueOf(i).isProbablePrime(10)) {
                if (legendre(n, i) == 1) {
                    factorBase.add(i);
                }
            }

        }
        return factorBase;
    }
}
