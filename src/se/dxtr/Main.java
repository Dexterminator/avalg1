package se.dxtr;

import java.math.BigInteger;

public class Main {

    public static void main(String[] args) {
        //BigInteger dex = new BigInteger("9110123735000000000000000000000000000000000000000000000000000000000001");
        BigInteger test = new BigInteger("123143254234122131245431514353453");
        System.out.println(test);
        System.out.println(PrimeUtils.millerRabin(test, 1));
        //System.out.println(PrimeUtils.pollardRho(dex));
        //System.out.println(dex);
    }
}
