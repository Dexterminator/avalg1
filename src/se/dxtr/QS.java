package se.dxtr;

import com.google.common.math.BigIntegerMath;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by dexter on 11/11/14.
 */
public class QS {
    public static int getB(BigInteger n) {
        double c = 3;
        double tempN = Math.log(n.doubleValue());
        double b = c*Math.exp(0.5*Math.sqrt(tempN*Math.log(tempN)));
        System.out.println(b);
        return (int) b;
    }

    public static BigInteger Q(BigInteger x, BigInteger n) {
        BigInteger parenthesis = BigIntegerMath.sqrt(n, RoundingMode.FLOOR).add(x);
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

    public static BigInteger[] getSieveArray (BigInteger n, int size){
        if(size < 100)
            size = 100;
        double tmpFirstX = Math.sqrt(n.doubleValue());
        BigInteger firstX = new BigDecimal(tmpFirstX).toBigInteger();
        firstX = BigInteger.valueOf(1);
        BigInteger[] sieveArray = new BigInteger[size];
        for (int x = 0; x < size; x++) {
            sieveArray[x] = Q(firstX.add(BigInteger.valueOf(x)), n);
        }
        return sieveArray;
    }

    public static ArrayList<Integer> performSieving(BigInteger[] sieveArray, ArrayList<Integer> factorBase, BigInteger n) {
        ArrayList<Integer> smoothIndices = new ArrayList<Integer>();
        for (Integer prime : factorBase) {
            double[] roots = tonelliShanks(n, prime);
            for (double root : roots) {
                double x = (root - Math.ceil(Math.sqrt(n.doubleValue()))) % prime;
                if (x < 0)
                    x += prime;
                sieveDivision(sieveArray, prime, (int) x, smoothIndices);
            }
        }
        return smoothIndices;
    }

    public static void sieveDivision(BigInteger[] sieveArray, Integer prime, int x, ArrayList<Integer> smoothIndices) {
        while (x < sieveArray.length) {
            sieveArray[x] = sieveArray[x].divide(BigInteger.valueOf(prime));
            while(sieveArray[x].mod(BigInteger.valueOf(prime)).equals(BigInteger.ZERO)){
                sieveArray[x] = sieveArray[x].divide(BigInteger.valueOf(prime));
            }
            if (sieveArray[x].equals(BigInteger.ONE)) {
                smoothIndices.add(x);
            }
            x += prime;
        }
    }

    public static int[][] getExpMatrix(ArrayList<Integer> factorBase, ArrayList<Integer> smoothIndices,
                                       BigInteger[] originalSieve) {
        // TODO: Make sure matrix is mod 2
        int[][] expMatrix = new int[smoothIndices.size()][factorBase.size() + smoothIndices.size()];
        for (int i = 0; i < smoothIndices.size(); i++) {
            int smoothIndex = smoothIndices.get(i);
            ArrayList<BigInteger> factors = PrimeUtils.pollardRho(originalSieve[smoothIndex]);
            for (BigInteger factor : factors) {
                expMatrix[i][factorBase.indexOf(factor.intValue())]++;
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

    public static ArrayList<Integer> processMatrix(int[][] matrix, int baseSize) {
        // Perform gaussian elimination
        for (int col = 0; col < baseSize; col++) {
            gaussForColumn(matrix, col);
        }
        printMatrix(matrix);

        // Find all zero rows
        ArrayList<Integer> zeroRows = new ArrayList<Integer>();
        for (int i = 0; i < matrix.length; i++) {
            if (isZeroRow(matrix[i], baseSize))
                zeroRows.add(i);
        }

        // Find indices corresponding to zero rows by traversing the identity matrix
        ArrayList<Integer> indices = new ArrayList<Integer>();
        for (Integer zeroRow : zeroRows) {
            for (int col = baseSize; col < matrix[zeroRow].length; col++) {
                if (matrix[zeroRow][col] == 1) {
                    indices.add(col - baseSize);
                }
            }
            System.out.println(indices);
        }
        System.out.println(zeroRows);
        return indices;
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
}
