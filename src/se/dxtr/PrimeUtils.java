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
            return gcd(b, a.mod(b));
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

    private static boolean helpMillerRabin(BigInteger n, Random r){
        BigInteger two = new BigInteger("2");
        BigInteger three = new BigInteger("3");
        // Base cases
        if(n.equals(BigInteger.ZERO) || n.equals(BigInteger.ONE)){
            return false;
        }

        if(n.equals(two) || n.equals(three)){
            return true;
        }

        if(n.mod(two).equals(BigInteger.ZERO)){
            return false;
        }

        BigInteger random = BigInteger.ZERO;


        do {
            random =  new BigInteger(n.bitLength()-1, r);
        } while(random.compareTo(BigInteger.ONE) <= 0);

        // Make sure this number is coprime with n
        if(!(gcd(n, random).equals(BigInteger.ONE))) {
            return false;
        }

        BigInteger t = n.subtract(BigInteger.ONE);
        int s = 0;
        // Find the number s and t as described in lecture notes
        while ((t.mod(two)).equals(BigInteger.ZERO)) {
            t = t.divide(two);
            s++;
        }
        BigInteger firstVal = modPow(random, t, n);


        /* According to the lecture notes, if firstVal == 1 then n is probably a prime */

        if(firstVal.equals(BigInteger.ONE))
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
     * Own version of the BigInteger.modPow function
     * @param a
     * @param b
     * @param n
     * @return
     */
    public static BigInteger modPow(BigInteger a, BigInteger b, BigInteger n){
        BigInteger temp =  BigInteger.ONE;
        for( BigInteger i = b; i.compareTo(BigInteger.ZERO) > 0; i = i.subtract(BigInteger.ONE)){
            temp = temp.multiply(temp);
            temp = temp.mod(n);
        }
        return temp;
    }


}
