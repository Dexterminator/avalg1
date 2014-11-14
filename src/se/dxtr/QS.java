package se.dxtr;

import com.google.common.math.BigIntegerMath;
import sun.security.util.BigInt;

import java.io.FileInputStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by dexter on 11/11/14.
 */
public class QS {
    private static BigInteger nRoot;
    public static int getB(BigInteger n) {
        double c = 3;
        double tempN = Math.log(n.doubleValue());
        double b = c*Math.exp(0.5 * Math.sqrt(tempN * Math.log(tempN)));
        return (int) b;
    }

    public static BigInteger Q(BigInteger x, BigInteger n) {
        BigInteger parenthesis = nRoot.add(x);
        BigInteger q = parenthesis.multiply(parenthesis).subtract(n);
        return q;
        //double parenthesis = Math.floor(Math.sqrt(n.doubleValue())) + x.doubleValue();
        //System.out.println(Math.pow(parenthesis, 2));
        //double q = Math.pow(parenthesis, 2) - n.doubleValue();
        //return new BigDecimal(q).toBigInteger();
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
        factorBase.add(2);
        for (int i = 2; i < b; i++) {
            if (BigInteger.valueOf(i).isProbablePrime(10)) {
                if (legendre(n, i) == 1) {
                    factorBase.add(i);
                }
            }
        }
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

    public static float[] getSieveArray (BigInteger n, int firstVal, int size){
        if(size < 100)
            size = 100;
        //TODO: CHECK IF THIS NEEDS TO BE ROUNDED
        double tmpFirstX = Math.ceil(Math.sqrt(n.doubleValue())) + firstVal;
        BigInteger firstX = new BigDecimal(tmpFirstX).toBigInteger();
        firstX = BigInteger.valueOf(0);
        float[] sieveArray = new float[size];
        nRoot = BigIntegerMath.sqrt(n, RoundingMode.FLOOR);
        for (int x = 0; x < size; x++) {
            sieveArray[x] = (float) Math.log(Q(firstX.add(BigInteger.valueOf(x)), n).doubleValue());
        }
        return sieveArray;
    }

    public static ArrayList<Integer> performSieving(float[] sieveArray, ArrayList<Integer> factorBase, BigInteger n) {
        ArrayList<Integer> smoothIndices = new ArrayList<Integer>();
        //System.out.println(sieveArray.length);
        /*
        //System.out.println(Arrays.toString(sieveArray));
        for (Integer prime : factorBase) {
            double[] roots = tonelliShanks(n, prime);
            for (double root : roots) {
                double x = (BigInteger.valueOf((long )root).subtract(BigIntegerMath.sqrt(n, RoundingMode.CEILING)).mod(BigInteger.valueOf(prime))).doubleValue();
                //double x = (root - Math.ceil(Math.sqrt(n.doubleValue()))) % prime;
                if (x < 0)
                    x += prime;
                sieveDivision(sieveArray, prime, (int) x, smoothIndices, n);
            }
        }
        */
        boolean first = true;
        int lastIndexPlusOne;
        while(factorBase.size() > smoothIndices.size()){
            if(first) {
                lastIndexPlusOne = 0;
                first = !first;
            } else {
                lastIndexPlusOne = sieveArray.length;
            }
            sieveArray = getSieveArray(n, lastIndexPlusOne, (int) Math.pow(factorBase.size(), 2));
            for (Integer prime : factorBase) {
                double[] roots = tonelliShanks(n, prime);
                for (double root : roots) {
                    double x = (BigInteger.valueOf((long )root).subtract(nRoot).mod(BigInteger.valueOf(prime))).doubleValue();
                    //double x = (root - Math.ceil(Math.sqrt(n.doubleValue()))) % prime;
                    if (x < 0)
                        x += prime;
                    sieveDivision(sieveArray, prime, (int) x, smoothIndices, n);
                }
            }
        }
        System.out.println(smoothIndices);
        return smoothIndices;
    }

    public static void sieveDivision(float[] sieveArray, Integer prime, int x, ArrayList<Integer> smoothIndices, BigInteger n) {
        //System.out.println(Arrays.toString(sieveArray));
        while (x < sieveArray.length) {
            //sieveArray[x] -= Math.log(prime);
            /*
            while(sieveArray[x].mod(BigInteger.valueOf(prime)).equals(BigInteger.ZERO)){
                sieveArray[x] = sieveArray[x].divide(BigInteger.valueOf(prime));
            }
            */
            //TODO: use x to calculate Q(x) again, check the mod. Done?
            BigInteger temp = Q(BigInteger.valueOf(x), n);
            while(temp.mod(BigInteger.valueOf(prime)).equals(BigInteger.ZERO)){
                sieveArray[x] -= Math.log(prime);
                temp = temp.divide(BigInteger.valueOf(prime));
            }
            if (sieveArray[x] < 0.1) {
                System.out.println(sieveArray[x]);
                smoothIndices.add(x);
            } else {
                //System.out.println(sieveArray[x]);
            }
            x += prime;
        }
    }

    public static int[][] getExpMatrix(ArrayList<Integer> factorBase, ArrayList<Integer> smoothIndices, BigInteger n) {
        int[][] expMatrix = new int[smoothIndices.size()][factorBase.size() + smoothIndices.size()];
        for (int i = 0; i < smoothIndices.size(); i++) {
            int smoothIndex = smoothIndices.get(i);
            // TODO: Use smoothIndex to calculate Q(x) again, as well as probably use Trial Division instead. Done?
            ArrayList<BigInteger> factors = PrimeUtils.trialDivision(Q(BigInteger.valueOf(smoothIndex), n));
            //ArrayList<BigInteger> factors = PrimeUtils.pollardRho(BigInteger.valueOf(Math.round(Math.exp(originalSieve[smoothIndex]))));
            //System.out.println(Math.exp(originalSieve[smoothIndex]));
            //System.out.println(factors);
            //System.out.println(factorBase);
            for (BigInteger factor : factors) {
                expMatrix[i][factorBase.indexOf(factor.intValue())] ^= 1;
            }
        }

        // Append the identity matrix to the exponent matrix
        int j = factorBase.size();
        for (int i = 0; i < expMatrix.length; i++) {
            expMatrix[i][j] = 1;
            j++;
        }
        return expMatrix;
    }

    public static ArrayList<Integer>[] processMatrix(int[][] matrix, int baseSize) {
        printMatrix(matrix);
        // Perform gaussian elimination
        for (int col = 0; col < baseSize; col++) {
            gaussForColumn(matrix, col);
        }
        //printMatrix(matrix);

        // Find all zero rows
        ArrayList<Integer> zeroRows = new ArrayList<Integer>();
        for (int i = 0; i < matrix.length; i++) {
            if (isZeroRow(matrix[i], baseSize))
                zeroRows.add(i);
        }

        // Find indices corresponding to zero rows by traversing the identity matrix
        ArrayList<Integer>[] subsets = new ArrayList[zeroRows.size()];
        for (int i = 0; i < zeroRows.size(); i++) {
            ArrayList<Integer> indices = new ArrayList<Integer>();
            for (int col = baseSize; col < matrix[zeroRows.get(i)].length; col++) {
                if (matrix[zeroRows.get(i)][col] == 1) {
                    indices.add(col - baseSize);
                }
            }
            subsets[i] = indices;
        }

        System.out.println(zeroRows);
        /*
        for(ArrayList<Integer> array: subsets){
            System.out.println(array);
        }
        */
        return subsets;
    }

    private static boolean isZeroRow (int[] row, int baseSize) {
        for (int i = 0; i < baseSize; i++) {
            if (row[i] != 0)
                return false;
        }
        return true;
    }

    public static void gaussForColumn (int[][] matrix, int column) {
        for (int row = 0; row < matrix.length; row++) {
            if (matrix[row][column] == 1) {
                performSubtractions(matrix, row, column);
                return;
            }
        }
    }

    private static void performSubtractions(int[][] matrix, int rowToSubtract, int column) {
        for (int row = 0; row < matrix.length; row++) {
            if (row != rowToSubtract && matrix[row][column] == 1)
                subtractRow(matrix, rowToSubtract, row);
        }
    }

    private static void subtractRow(int[][] matrix, int row1, int row2) {
        for (int i = 0; i < matrix[row2].length; i++) {
            matrix[row2][i] ^= matrix[row1][i];
        }
    }

    private static void printMatrix(int[][] matrix) {
        for (int[] ints : matrix) {
            System.out.println(Arrays.toString(ints));
        }
    }

    public static BigInteger getNonTrivialFactor(ArrayList<Integer>[] subsets, ArrayList<Integer> smoothIndices, BigInteger n) {
        for (ArrayList<Integer> subset : subsets) {
            BigInteger factor = getFactor(subset, smoothIndices, n);
            if (!factor.equals(BigInteger.ONE) && !factor.equals(n))
                return factor;
        }
        return null;
    }

    /**
     * Here be dragons
     *
     * This code may summon Cthulhu if used improperly
     * TODO: Remove unnecessary Math.exp, use smoothIndicies instead.
     */
    private static BigInteger getFactor(ArrayList<Integer> subsetIndices, ArrayList<Integer> smoothIndices, BigInteger n) {
        //int root = (int) Math.ceil(Math.sqrt(n.doubleValue()));
        BigInteger temp = BigInteger.valueOf(smoothIndices.get(subsetIndices.get(0))).add(nRoot);
        BigInteger x = temp.multiply(temp);
        //double x = Math.pow(smoothIndices.get(subsetIndices.get(0)) + root, 2);
        BigInteger y = Q(BigInteger.valueOf(smoothIndices.get(subsetIndices.get(0))), n);
        //BigInteger y = BigInteger.valueOf((long) Math.round(Math.exp(originalSieve[smoothIndices.get(subsetIndices.get(0))])));
        for (int i = 1; i < subsetIndices.size(); i++) {
            BigInteger temp2 = BigInteger.valueOf(smoothIndices.get(subsetIndices.get(i))).add(nRoot);
            temp2 = temp2.multiply(temp2);
            x = x.multiply(temp2);
            y = y.multiply(Q(BigInteger.valueOf(smoothIndices.get(subsetIndices.get(i))), n));
            //y = y.multiply(BigInteger.valueOf((long) Math.round(Math.exp(originalSieve[smoothIndices.get(subsetIndices.get(i))]))));
        }
        // x congruent with y (mod n) at this point, sqrt to get the values to calculate gcd

        BigInteger a = BigIntegerMath.sqrt(x, RoundingMode.UNNECESSARY).mod(n);
        BigInteger b = BigIntegerMath.sqrt(y, RoundingMode.UNNECESSARY).mod(n);
        System.out.println(a);
        System.out.println(b);

        BigInteger factor = PrimeUtils.gcd(a.subtract(b), n);
        return factor;
    }
}
