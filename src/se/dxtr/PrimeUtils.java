package se.dxtr;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Random;

/**
 * Created by Ludde on 2014-11-04.
 */
public class PrimeUtils {
    private final static SecureRandom random = new SecureRandom();
    private final static SecureRandom sRand = new SecureRandom();

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
            if(possiblePrime.isProbablePrime(5)){
                factors.add(possiblePrime);
                n = n.divide(possiblePrime);
                if(n.equals(BigInteger.ONE)){
                    break;
                } else if(n.isProbablePrime(5)){
                    // In this case, the resulting number is a factor
                    factors.add(n);
                    break;
                } else {
//                    System.out.println(n);
                }
                possiblePrime = helpPollardRho(n);
            } else {
                possiblePrime = helpPollardRho(possiblePrime);
            }
        }

        return factors;
    }

    public static BigInteger random(BigInteger n){
        BigInteger r;
        do {
            r = new BigInteger(n.bitLength(), sRand);
        } while (r.compareTo(n) >= 0);
        return r;
    }

    public static BigInteger minFunc (BigInteger a, BigInteger b) {
        if (a.compareTo(b) <= 0) {
            return a;
        } else return b;
    }
    public static BigInteger helpBrent(BigInteger n) {
        BigInteger ys = BigInteger.ONE;
        BigInteger x = BigInteger.ONE;

        if (n.mod(BigInteger.valueOf(2)).equals(BigInteger.ZERO)) {
            return BigInteger.valueOf(2);
        }
        BigInteger y = random(n);
        BigInteger c = random(n);
        BigInteger m = random(n);

        BigInteger g = BigInteger.ONE;
        BigInteger r = BigInteger.ONE;
        BigInteger q = BigInteger.ONE;

        while (g.equals(BigInteger.ONE)) {
            x = y;
            for (BigInteger i = BigInteger.ZERO; i.compareTo(r) <= 0; i = i.add(BigInteger.ONE)) { // zero == one??
                y = y.multiply(y).mod(n).add(c).mod(n);
            }
            BigInteger k = BigInteger.ZERO;
            while (k.compareTo(r) < 0 && g.equals(BigInteger.ONE)) {
                ys = y;
                BigInteger minRK = minFunc(m, r.subtract(k));
                for (BigInteger i = BigInteger.ZERO; i.compareTo(minRK) < 0; i = i.add(BigInteger.ONE)) { //zero = one??
                    y = y.multiply(y).mod(n).add(c).mod(n);
                    q = q.multiply(x.subtract(y).abs()).mod(n);
                }
                g = gcd(q, n);
                k = k.add(m);
            }
            r = r.multiply(BigInteger.valueOf(2));
        }

        if (g.equals(n)) {
            while (true) {
                ys = ys.multiply(ys).mod(n).add(c).mod(n);
                g = gcd(x.subtract(ys).abs(), n);
                if (g.compareTo(BigInteger.ONE) > 0) {
                    break;
                }
            }
        }
        return g;
    }

    public static ArrayList<BigInteger> brent(BigInteger n){
        ArrayList<BigInteger> factors = new ArrayList<BigInteger>();
        if(n.isProbablePrime(10)){
            factors.add(n);
            factors.add(BigInteger.valueOf(1));
        }
        BigInteger possiblePrime = helpBrent(n);
        while(true){
            if(possiblePrime.isProbablePrime(5)){
                factors.add(possiblePrime);
                n = n.divide(possiblePrime);
                if(n.equals(BigInteger.ONE)){
                    break;
                } else if (n.isProbablePrime(5)){
                    factors.add(n);
                    break;
                }
                possiblePrime = helpBrent(n);
            }
            else {
                possiblePrime = helpBrent(possiblePrime);
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

        BigInteger firstVal = modPow2(random, t, n);


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
            iterVal = modPow2(iterVal, two, n);
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

    public static BigInteger modPow2(BigInteger a, BigInteger b, BigInteger n){
        BigInteger z = b;
        BigInteger y = a;
        BigInteger x = BigInteger.ONE;
        BigInteger two = new BigInteger("2");
        while(!z.equals(BigInteger.ZERO)){
            if(z.mod(two).equals(BigInteger.ZERO)){
                z = z.divide(two);
                y = y.pow(2).mod(n);
            } else {
                z = z.subtract(BigInteger.ONE);
                x = x.multiply(y).mod(n);
            }
        }
        return x;
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

    /**
     * A quick trial division algorithm
     * @param a
     * @return
     */
    public static void trialDivision(BigInteger a, ArrayList<BigInteger> retList, int fileNum) throws FileNotFoundException{
        FileInputStream file = null;

        file = new FileInputStream("src/se/dxtr/primes"+fileNum+".txt");
        Kattio io = new Kattio(file);
        while(io.hasMoreTokens()){
            if(a.equals(BigInteger.ONE))
                break;
            BigInteger temp = BigInteger.valueOf(io.getInt());
            while(a.mod(temp).equals(BigInteger.ZERO)){
                a = a.divide(temp);
                retList.add(temp);
            }

        }
        if(a.isProbablePrime(10)){
            retList.add(a);
            a = a.divide(a);
        }
        if(!a.equals(BigInteger.ONE)){
            trialDivision(a, retList, fileNum+1);
        }

    }


}
