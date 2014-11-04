package se.dxtr;

import java.math.BigInteger;

public class Main {

    public static void main(String[] args) {
        BigInteger dex = new BigInteger("9110123735000000000000000000000000000000000000000000000000000000000001");
        System.out.println(PrimeUtils.pollardRho(dex));
        System.out.println(dex);
    }
}
