package se.dxtr;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;

/**
 * Created by Ludde on 2014-11-25.
 */
public class BrentThread extends Thread{
    private final SecureRandom sRand = new SecureRandom();
    private BigInteger N;
    private BigInteger start;
    private BigInteger stop;
    public ArrayList<BigInteger> ans;

    public BrentThread(BigInteger N){
        this.N = N;
    }

    @Override
    public void run(){
        ans = brent(N);
    }

    public BigInteger gcd(BigInteger a, BigInteger b){
        if(b.equals(BigInteger.ZERO)){
            return a;
        }
        return gcd(b, a.mod(b));
    }

    public BigInteger random(BigInteger n){
        BigInteger r;
        do {
            r = new BigInteger(n.bitLength(), new SecureRandom());
        } while (r.compareTo(n) >= 0);
        return r;
    }

    public BigInteger random(BigInteger min, BigInteger max, BigInteger n){
        BigInteger r = random(n.add(BigInteger.ONE));
        r = r.mod(max.subtract(min).add(BigInteger.ONE));
        r = r.add(min);
        return r;
    }

    public BigInteger minFunc (BigInteger a, BigInteger b) {
        if (a.compareTo(b) <= 0) {
            return a;
        } else return b;
    }

    public BigInteger helpBrent(BigInteger n) {
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
            for (BigInteger i = BigInteger.ZERO; i.compareTo(r) <= 0; i = i.add(BigInteger.ONE)) {
                y = y.multiply(y).mod(n).add(c).mod(n);
            }
            BigInteger k = BigInteger.ZERO;
            while (k.compareTo(r) < 0 && g.equals(BigInteger.ONE)) {
                ys = y;
                BigInteger minRK = minFunc(m, r.subtract(k));
                for (BigInteger i = BigInteger.ZERO; i.compareTo(minRK) < 0; i = i.add(BigInteger.ONE)) {
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

    public ArrayList<BigInteger> brent(BigInteger n){
        ArrayList<BigInteger> factors = new ArrayList<BigInteger>();
        if(PrimeUtils.millerRabin(n, 10)){
            factors.add(n);
            factors.add(BigInteger.valueOf(1));
        }
        BigInteger possiblePrime = helpBrent(n);
        while(true){
            if(PrimeUtils.millerRabin(possiblePrime, 5)){
                factors.add(possiblePrime);
                n = n.divide(possiblePrime);
                if(n.equals(BigInteger.ONE)){
                    break;
                } else if (PrimeUtils.millerRabin(n, 5)){
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

}
