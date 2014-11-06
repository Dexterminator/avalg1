package se.dxtr;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Random;

/**
 * Created by Ludde on 2014-11-04.
 */
public class PrimeUtils {
    private final static SecureRandom random = new SecureRandom();

    public static BigInteger gcd(BigInteger a, BigInteger b) {
        if(b.equals(BigInteger.ZERO)){
            return a;
        }
        return gcd(b, a.mod(b));
    }

    private static BigInteger g(BigInteger x, BigInteger n, BigInteger add){
        return x
                .multiply(x)
                .add(add)
                .mod(n);
    }

    public static BigInteger helpPollardRho(BigInteger n) {

        BigInteger d = BigInteger.ONE;
        BigInteger add = new BigInteger(n.bitLength(), random);
        BigInteger x = new BigInteger(n.bitLength(), random);
        BigInteger y = x;
        /*
        if(n.mod(new BigInteger("2")).equals(BigInteger.ZERO)){
            return new BigInteger("2");
        }
        */
        if(n.mod(new BigInteger("2")).equals(BigInteger.ZERO)){
            return new BigInteger("2");
        }
        while (d.equals(BigInteger.ONE)){
            x = g(x, n, add);
            y = g(g(y, n, add), n, add);
            d = gcd(y.subtract(x).abs(), n);
        }
        /*
        if (d.equals(n))
            return new BigInteger("-1");
            */
        return d;
    }

    public static ArrayList<BigInteger> pollardRho(BigInteger n){
        ArrayList<BigInteger> factors = new ArrayList<BigInteger>();
        BigInteger possiblePrime = helpPollardRho(n);
        while(true){
            if(millerRabin(possiblePrime, 20)){
                System.out.println(possiblePrime);
                factors.add(possiblePrime);
                n = n.divide(possiblePrime);
                if(n.equals(BigInteger.ONE)){
                    break;
                }
                possiblePrime = helpPollardRho(n);
            } else {
                possiblePrime = helpPollardRho(possiblePrime);
            }
        }

        return factors;
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
        return temp.mod(n);
    }


}
