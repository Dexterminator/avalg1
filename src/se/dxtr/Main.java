package se.dxtr;

import java.io.*;
import java.math.BigInteger;
import java.util.*;

public class Main {
    // Change this to the instance of Pollard-Rho running
    private static final int CORE_NUMBER = 1;
    // Change this to your name
    private static final String GROUP_MEMBER = "ludjan";
    // Time limit per number
    //private final static long TIME_LIMIT = 1000; // One second
    private final static long TIME_LIMIT = 1000*60*5; // Five minutes
    private final static long SECOND_TIME_LIMIT = 1000*60*60*10; // Ten hours
    // Personnummer of this group member
    private final static BigInteger personnummer = new BigInteger("9112232872");
    // The long personnummer of this group member
    private final static BigInteger longPersonnummer = new BigInteger("9112232872000000000000000000000000000000000000000000000000000000000000");

    public static void main(String[] args) {
        ArrayList<BigInteger> restList = new ArrayList<BigInteger>();

        /*
        This fucking piece of shit is the failed attempt to implement the Quadratic sieve. May it never be touched again.
        BigInteger personnummer = new BigInteger("9112232872");
        //BigInteger longPersonnummer = new BigInteger("9112232872000000000000000000000000000000000000000000000000000000000001");
        BigInteger longPersonnummer = new BigInteger("911223287200000000000000000001");

        BigInteger n = new BigInteger("1534000007");
        n = longPersonnummer;
        System.out.println("N: " + n + "\n");

        ArrayList<Integer> factorBase = QS.factorBase(n, QS.getB(n));
        float[] sieve = QS.getSieveArray(n, 0, (int) Math.round(Math.pow(factorBase.size(), 2)));
        ArrayList<Integer> smoothIndices = QS.performSieving(sieve, factorBase, n);
        BitSet[] expMatrix = QS.getExpMatrix(factorBase, smoothIndices, n);
        ArrayList<Integer>[] subsets = QS.processMatrix(expMatrix, factorBase.size());
        BigInteger factor = QS.getNonTrivialFactor(subsets, smoothIndices, n);
        if (factor != null) {
            System.out.println(factor);
        }
        */
        // factor may be null here if we fucked up

//        for (Integer index : indices) {
//            int x = index + root;
//            BigInteger y = originalSieve[smoothIndices.get(index)];
//        }
        int j = 0;
        initialPrint(personnummer, j);
        // Change the index for each core
        int numFound = 0;
        for(long i = 1; i < 101; i++){
            BigInteger number = longPersonnummer.add((BigInteger.valueOf(i)));
            System.out.println("******** NEW NUMBER: " + number + " ********");
            boolean found = false;
            BrentThread thread = new BrentThread(number);
            ArrayList<BigInteger> factors = new ArrayList<BigInteger>();
            long startTime = System.currentTimeMillis();
            long totalTime = TIME_LIMIT;
            thread.start();
            while(!found && (System.currentTimeMillis() - startTime) < totalTime){
                try{
                    Thread.sleep(500L);
                } catch (InterruptedException e){
                    return;
                }
                if(thread.ans != null){
                    factors = thread.ans;
                    numFound ++;
                    found = true;
                }
                // Shut down thread
                if((System.currentTimeMillis() - startTime) > totalTime){
                    thread.interrupt();
                }
            }

            if(factors == null || factors.size() < 1){
                restList.add(number);
                continue; // Iterate to next number
            }
            if(factors != null && factors.size() == 1){
                System.out.println(number + " is prime!");
            }
            LinkedHashMap<BigInteger, Integer> matrix = new LinkedHashMap<BigInteger, Integer>();
            for(BigInteger factor : factors){
                if(matrix.containsKey(factor)){
                    int num = matrix.get(factor) + 1;
                    matrix.put(factor, num);
                } else {
                    matrix.put(factor, 1);
                }
            }


            writeToFile(number, matrix);
        }
        System.out.println(numFound);
        for(BigInteger number: restList){
            System.out.println("******** NEW NUMBER: " + number + " ********");
            boolean found = false;
            BrentThread thread = new BrentThread(number);
            ArrayList<BigInteger> factors = new ArrayList<BigInteger>();
            long startTime = System.currentTimeMillis();
            long totalTime = SECOND_TIME_LIMIT;
            thread.run();
            while(!found && (System.currentTimeMillis() - startTime) < totalTime){
                try{
                    Thread.sleep(500L);
                } catch (InterruptedException e){
                    return;
                }
                if(thread.ans != null){
                    factors = thread.ans;
                    found = true;
                }
            }

            if(factors == null || factors.size() < 1){
                restList.remove(number);
                continue; // Iterate to next number
            }
            if(factors != null && factors.size() == 1){
                System.out.println(number + " is prime!");
            }
            LinkedHashMap<BigInteger, Integer> matrix = new LinkedHashMap<BigInteger, Integer>();
            for(BigInteger factor : factors){
                if(matrix.containsKey(factor)){
                    int num = matrix.get(factor) + 1;
                    matrix.put(factor, num);
                } else {
                    matrix.put(factor, 1);
                }
            }


            writeToFile(number, matrix);
        }



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
    private static void writeToFile(BigInteger number, LinkedHashMap<BigInteger, Integer> matrix){
        try {
            PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(GROUP_MEMBER + "_output_core_"+String.valueOf(CORE_NUMBER)+".in", true)));
            Iterator it = matrix.entrySet().iterator();
            StringBuilder sb = new StringBuilder();
            sb.append(number + ": ");
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
