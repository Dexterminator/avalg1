package se.dxtr;

import java.math.BigInteger;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) {
        BigInteger dex = new BigInteger("9110123735000000000000000000000000000000000000000000000000000000000001");

        for(long i = 0; i < 13; i++){

            BigInteger number = dex.add((BigInteger.valueOf(i)));
            System.out.println("******** NEW NUMBER: " + number + " ********");
            ArrayList<BigInteger> factors = PrimeUtils.pollardRho(number);
            if(factors.size() == 1){
                System.out.println(number + " is prime!");
            } else {
                System.out.println("The factors for " + number + " are " + factors.toString());
            }
        }
        BigInteger number = new BigInteger("1231242341234");

        //PrimeUtils.pollardRho(new BigInteger("25"));
        //System.out.println(PrimeUtils.helpPollardRho(new BigInteger("4")));


        //System.out.println(factors.toString());
        //System.out.println(PrimeUtils.pollardRho(dex));
        //System.out.println(dex);
    }
}
