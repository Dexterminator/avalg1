package se.dxtr;

import java.math.BigInteger;
import java.util.Random;

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

    public static boolean helpMillerRabin(BigInteger n, Random r){
        BigInteger random = new BigInteger(n.bitLength()-1, r).mod(n); // Pick random number

        // Make sure this number is coprime with n
        if(!gcd(n, random).equals(BigInteger.ONE))
            return false;

        BigInteger t = n.subtract(BigInteger.ONE);
        int s = 0;
        BigInteger two = new BigInteger("2");
        // Find the number s and t as described in lecture notes
        while (t.mod(two).equals(BigInteger.ZERO)) {
            t = t.divide(two);
            s++;
        }

        BigInteger firstVal = modPow(random, t, n);

        /* According to the lecture notes, if firstVal == 1 then n is probably a prime */

        if(firstVal.mod(n).equals(BigInteger.ONE))
            return true;

        /*
            Also according to the lecture notes, if any of the numbers
            firstVal^i mod n == -1 (which is equivalent with n-1),
            for 1 <= i <= s, then n is probably prime
         */
        BigInteger iterVal = firstVal;
        for( int i = 0; i < s; i++){
            if(iterVal.equals(n.subtract(BigInteger.ONE))){
                return true;
            }
            // The next value in the sequence is this number raised in two mod n
            iterVal = modPow(iterVal, two, n);
        }

        /* If none of the numbers in the sequence are equal to n-1,
         then the number is definitely composite
         */
        return false;

    }

    public static boolean millerRabin(BigInteger n, int k){
        Random r = new Random();
        for(int i = 0; i < k; i++){
            /* If we receive false from the help function,
             n is definitely composite
            */
            if(!helpMillerRabin(n, r))
                return false;
        }

        // If all of the Miller Rabin tests pass, n is probably prime
        return true;
    }

    /**
     * Own version of the modPow function
     * @param a
     * @param b
     * @param n
     * @return
     */
    public static BigInteger modPow(BigInteger a, BigInteger b, BigInteger n){
        BigInteger temp =  a;
        for( int i = 0; i < b.intValue(); i++){
            temp = temp.multiply(temp);
            temp = temp.mod(n);
        }

        return temp.mod(n);
    }


}
