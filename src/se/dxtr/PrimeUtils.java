package se.dxtr;

import java.math.BigInteger;

/**
 * Created by Ludde on 2014-11-04.
 */
public class PrimeUtils {

    public static BigInteger gcd(BigInteger a, BigInteger b){
        if(b.equals("0")){
            return a;
        } else {
            return gcd(a, a.mod(b));
        }
    }


}
