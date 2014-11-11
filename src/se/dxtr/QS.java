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

    public static double legendre(BigInteger N, int p) {
        double result = 1;
        double a = N.mod(BigInteger.valueOf(p)).doubleValue();
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
        factorBase.add(2);
        return factorBase;
    }

    public static double[] tonelliShanks(BigInteger n, int p){
        if(p == 2){
            if(n.mod(BigInteger.valueOf(2)).equals(BigInteger.valueOf(1))){
                return new double[]{1};
            } else{
                return new double[]{0};
            }
        }
        int Q = p-1;
        int S = 0;
        BigInteger R;
        while (Q % 2 == 0) {
            S++;
            Q /= 2;
        }
        if (S == 1) {
            // Case when p mod 4 = 3
            R = PrimeUtils.modPow2(n, BigInteger.valueOf((p + 1) / 4), BigInteger.valueOf(p));
        } else {
            double res = 0;
            int Z = 2;
            for (Z = 2; ; Z++) {
                res = legendre(BigInteger.valueOf(Z), p);
                if (res == -1) {
                    break;
                }
            }
            /*
            while(res != -1){
                res = legendre(BigInteger.valueOf(Z), p);
                Z++;
            }
            */
            BigInteger C = PrimeUtils.modPow2(BigInteger.valueOf(Z), BigInteger.valueOf(Q), BigInteger.valueOf(p));
            R = PrimeUtils.modPow2(n, BigInteger.valueOf((Q + 1) / 2), BigInteger.valueOf(p));
            BigInteger t = PrimeUtils.modPow2(n, BigInteger.valueOf(Q), BigInteger.valueOf(p));
            int M = S;

            while (!t.equals(BigInteger.ONE)) {
                int tempval = 1;
                int i = 1;
                while (i < M) {
                    tempval *= 2;

                    if (PrimeUtils.modPow2(t, BigInteger.valueOf(tempval),
                            BigInteger.valueOf(p)).equals(BigInteger.ONE)) {
                        break;
                    }
                    i++;

                }
                if (i == M && i != 1) {
                    i--;
                }
                BigInteger b = PrimeUtils.modPow2(C, BigInteger.valueOf(2).pow(M - i - 1), BigInteger.valueOf(p));
                R = R.multiply(b).mod(BigInteger.valueOf(p));
                t = t.multiply(b).multiply(b).mod(BigInteger.valueOf(p));
                C = b.multiply(b).mod(BigInteger.valueOf(p));
                M = i;
            }
        }
        double[] ret = new double[]{R.doubleValue(), p-R.doubleValue()};
        return ret;
    }

    public static BigInteger[] getSieveArray (BigInteger n, int size){
        double tmpFirstX = Math.sqrt(n.doubleValue());
        BigInteger firstX = new BigDecimal(tmpFirstX).toBigInteger();
        firstX = BigInteger.valueOf(1);
        BigInteger[] sieveArray = new BigInteger[size];
        for (int x = 0; x < size; x++) {
            sieveArray[x] = Q(firstX.add(BigInteger.valueOf(x)), n);
        }
        return sieveArray;
    }

    public static BigInteger[] performSieving (BigInteger[] sieveArray, ArrayList<Integer> factorBase, BigInteger n) {
        for (Integer prime : factorBase) {
            double[] roots = tonelliShanks(n, prime);
            for (double root : roots) {
                double x = root - Math.ceil(Math.sqrt(n.doubleValue())) % prime;
                sieveDivision(sieveArray, prime, (int) x);
            }
        }
        return sieveArray;
    }

    public static void sieveDivision(BigInteger[] sieveArray, Integer prime, int x) {
        while (x < sieveArray.length) {
            sieveArray[x] = sieveArray[x].divide(BigInteger.valueOf(prime));
            x += prime;
        }
    }
}
