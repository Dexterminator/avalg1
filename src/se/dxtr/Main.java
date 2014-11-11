package se.dxtr;

import java.io.*;
import java.math.BigInteger;
import java.util.*;

public class Main {
    // Change this to the instance of Pollard-Rho running
    private static final int CORE_NUMBER = 1;
    // Change this to your name
    private static final String GROUP_MEMBER = "ludjan";

    public static void main(String[] args) {

        BigInteger personnummer = new BigInteger("9112232872");
        BigInteger longPersonnummer = new BigInteger("9112232872000000000000000000000000000000000000000000000000000000000001");
        BigInteger n = BigInteger.valueOf(15347);
        System.out.println("N: " + n + "\n");

        ArrayList<Integer> factorBase = QS.factorBase(n, QS.getB(n));
        BigInteger[] sieve = QS.getSieveArray(n, 100);
        BigInteger[] originalSieve = Arrays.copyOf(sieve, sieve.length);
        int[][] yFactors = new int[sieve.length][factorBase.size()];
        ArrayList<Integer> smoothIndices = QS.performSieving(sieve, factorBase, n, yFactors);

//        for (BigInteger bigInteger : sieve) {
//            System.out.println(bigInteger.toString());
//        }

        System.out.println("Indices of smooth y values: " + smoothIndices);
        for (Integer smoothIndex : smoothIndices) {
            System.out.println("Smooth number: " + originalSieve[smoothIndex]);
            System.out.println(Arrays.toString(yFactors[smoothIndex]));
            for (int i = 0; i < yFactors[smoothIndex].length; i++) {
                for (int primeCount = 0; primeCount < yFactors[smoothIndex][i]; primeCount++) {
                    System.out.println(factorBase.get(i));
                }
            }
            System.out.println();
        }
        ArrayList<BigInteger> pollardFactors = PrimeUtils.pollardRho(longPersonnummer);
        System.out.println("Pollar factors of: " + longPersonnummer);
        System.out.println(pollardFactors.toString());
//        int j = 0;
//        initialPrint(personnummer, j);
//        // Change the index for each core
//        for(long i = 50; i < 75; i++){
//
//            BigInteger number = longPersonnummer.add((BigInteger.valueOf(i)));
//            System.out.println("******** NEW NUMBER: " + number + " ********");
//            ArrayList<BigInteger> factors = PrimeUtils.pollardRho(number);
//            if(factors.size() == 1){
//                System.out.println(number + " is prime!");
//            } else {
//                System.out.println("The factors for " + number + " are " + factors.toString());
//            }
//            LinkedHashMap<BigInteger, Integer> matrix = new LinkedHashMap<BigInteger, Integer>();
//            for(BigInteger factor : factors){
//                if(matrix.containsKey(factor)){
//                    int num = matrix.get(factor) + 1;
//                    matrix.put(factor, num);
//                } else {
//                    matrix.put(factor, 1);
//                }
//            }
//
//
//            writeToFile(matrix);
//        }



        //PrimeUtils.pollardRho(new BigInteger("25"));
        //System.out.println(PrimeUtils.helpPollardRho(new BigInteger("4")));


        //System.out.println(factors.toString());
        //System.out.println(PrimeUtils.pollardRho(dex));
        //System.out.println(dex);
    }

    private static void initialPrint(BigInteger personnummer, int j){
        try {
            PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(GROUP_MEMBER + "_output_core_"+String.valueOf(CORE_NUMBER)+".in", true)));
            writer.println(personnummer + " " + j);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
    private static void writeToFile(LinkedHashMap<BigInteger, Integer> matrix){
        try {
            PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(GROUP_MEMBER + "_output_core_"+String.valueOf(CORE_NUMBER)+".in", true)));
            Iterator it = matrix.entrySet().iterator();
            StringBuilder sb = new StringBuilder();
            while(it.hasNext()){
                Map.Entry pairs = (Map.Entry)it.next();
                sb.append(pairs.getKey() + " " + pairs.getValue() + " ");
                it.remove();
            }
            writer.println(sb.toString());
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
