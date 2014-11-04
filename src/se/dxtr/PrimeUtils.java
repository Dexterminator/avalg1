package se.dxtr;

import java.math.BigInteger;

/**
 * Created by Ludde on 2014-11-04.
 */
public class PrimeUtils {

    public static BigInteger gcd(BigInteger a, BigInteger b){
        if(b.equals(BigInteger.ZERO)){
            return a;
        } else {
            return gcd(a, a.mod(b));
        }
    }

    public static BigInteger pollardRho(BigInteger n) {
        BigInteger d = new BigInteger("1");
        BigInteger x = new BigInteger("2");
        BigInteger y;
        while (d.equals(BigInteger.ONE)){
            x = x.pow(2).add(BigInteger.ONE).mod(n);
            y = x.pow(2).add(BigInteger.ONE).mod(n);
            d = gcd(x.subtract(y).abs(), n);
        }
        if (d.equals(n))
            return new BigInteger("-1");
        return d;
    }


}
