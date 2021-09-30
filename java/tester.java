// import java.io.*;
import java.util.Arrays;
// import java.lang.Math;
import linalgexercise.*;

class tester {

  static RealSquare uniqueA = new RealSquare(new double[][]{
      {1, 2, 3, 4, 5},
      {3, 5, 11, -7, 9},
      {2, -11, 5, 2, 13}, 
      {5, 7, -3, 8, 4},
      {4, -1, 7, -5, 17}});
  static RealVector uniqueb = new RealVector(new double[]{-12, 57.5, 30, -54, 49.5});

  public static void basicOperationsTest() {
    RealMatrix A = new RealMatrix(new double[][]{{1, 2, 3}, {3, 4, 5}, {5, 6, 7}});
    System.out.println(A);
    RealVector v = new RealVector(new double[]{5, 12});
    System.out.println(v);
    System.out.println(v.modulus());
    System.out.println(v.normalize());
    RealVector v2 = new RealVector(new double[]{1, 2});
    System.out.println(RealUtil.dot(v, v2));
    RealMatrix B = new RealMatrix(new double[][]{{3, 4},{5, 6}});
    System.out.println(A.dot(B));
    System.out.println(A.dot(v));
    System.out.println(A.ref());
    //System.out.println(D.inverse());
    //System.out.println(RealUtil.add(A, D));
    //System.out.println(RealUtil.sub(A, D));
    RealVector b = new RealVector(new double[]{1, 2, 3});
    System.out.println(b.scalarMulti(2));
  }

  public static void rrefTest() {
    RealMatrix C = new RealMatrix(
        new double[][]{{1, 2, 2, 3}, {2, -2, -8, 4}, {1, 1, 0, 1}, {0, 2, 4, 1}});
    System.out.println(C.ref());
    System.out.println(new RealRref(C));
    System.out.println(C.rank());
  }

  public static void linSolveGaussianTest() {
    RealMatrix[] matricesA = new RealMatrix[3];
    RealVector[] vectorsb = new RealVector[3];
    matricesA[0] = new RealMatrix(new double[][]{{2, 1, -1}, {1, -2, -3}, {5, -1, 2}});
    vectorsb[0] = new RealVector(new double[]{5, 0, 1});
    matricesA[1] = new RealMatrix(new double[][]{{1, 1, -2}, {3, -3, 2}, {7, -2, -2}});
    vectorsb[1] = new RealVector(new double[]{-12, 22, 3});
    matricesA[2] = new RealMatrix(new double[][] {{6, 15, -21}, {8, 2, -11}, {-2, -5, 7}});
    vectorsb[2] = new RealVector(new double[]{4, 15, 17});
    for (int i = 0; i < 3; i++) {
      System.out.println("A = ");
      System.out.println(matricesA[i]);
      System.out.println("b = ");
      System.out.println(vectorsb[i]);
      RealVector solution = RealUtil.linSolveGaussian(matricesA[i], vectorsb[i]);
      if (solution != null){
        System.out.println(solution);
      }
    }
  }

  public static void lduDecompositionTest() {
    RealSquare E = new RealSquare(new double[][]{{2, 1, 0},{1, 2, 1},{0, 1, 2}});
    System.out.println(E.lduDecomposition()[0]);
    System.out.println(E.lduDecomposition()[1]);
    System.out.println(E.lduDecomposition()[2]);
    System.out.println(RealUtil.dot(E.luDecomposition()[0], E.luDecomposition()[1]));
  }

  public static void luSolveTest() {
    RealSquare A = new RealSquare
        (new double[][]{{1, 1, 1}, {1, 2, 3}, {1, 3, 6}});
    RealVector b = new RealVector(new double[]{5, 7, 11});
    // System.out.println(A.luDecomposition()[0]);
    // System.out.println(A.luDecomposition()[1]);
    System.out.println(RealUtil.luSolve(A, b));
    System.out.println(A.luInverse());
    }

  public static void plduTest() {
/*  RealMatrix A = new RealMatrix
        (new double[][]{{0, 1, 1}, {1, 2, 1}, {2, 7, 9}});*/
    RealSquare A = new RealSquare
        (new double[][]{{7, 0, 5, 1}, {-2, 2, 4, -5}, {-3, 1, 3, 6}, {1, -6, 2, -4}});
    RealPldu factors = new RealPldu(A);
    System.out.println(factors.getL());
    System.out.println(factors.getDU());
    System.out.println(factors.getP());
    System.out.println(factors.getD());
    System.out.println("LU: ");
    RealMatrix lu = RealUtil.dot(factors.getL(), factors.getDU());
    System.out.println(lu);
    System.out.println("PA: ");
    RealMatrix pa = RealUtil.dot(factors.getP(), A);
    System.out.println(pa);
    if (RealUtil.matrixCompare(lu, pa)) {
      System.out.println("Success!");
    }
    System.out.println(factors.restore());
    if (RealUtil.matrixCompare(A, factors.restore())) {
      System.out.println("Greater success!");
    }
    System.out.println(factors.inverse());
    System.out.println(A.dot(factors.inverse()));
    if (RealUtil.matrixCompare(A.dot(factors.inverse()), RealUtil.eye(A.dim()[0]))) {
      System.out.println("One more success!");
    }
  }

  public static void plduSolveTest() {
    RealSquare A = new RealSquare
        (new double[][]{{1, 1, 1}, {1, 2, 3}, {1, 3, 6}});
    RealVector b = new RealVector(new double[]{5, 7, 11});
    System.out.println("Start here: ");
    System.out.println(new RealPldu(A).solve(b));
    System.out.println(new RealPldu(A).inverse());
    System.out.println(A.dot(new RealPldu(A).inverse()));
    System.out.println(new RealPldu(uniqueA).solve(uniqueb));
  }

  public static void basisTest() {
    RealMatrix A = new RealMatrix(new double[][]{{1, 2, 3, 0, 0}, {4, 10, 0, 0, 1}});
    RealVector[] rowBasis = A.rowBasis();
    RealVector[] colBasis = A.colBasis();
    RealVector[] nullBasis = A.nullBasis();
    int i;
    System.out.println("Row Basis: ");
    for (i = 0; i < rowBasis.length; i++) {
      System.out.println(rowBasis[i]);
    }
    System.out.println("Column Basis: ");
    for (i = 0; i < colBasis.length; i++) {
      System.out.println(colBasis[i]);
    }
    System.out.println("Nullspace(Kernel) Basis: ");
    for (i = 0; i < nullBasis.length; i++) {
      System.out.println(nullBasis[i]);
    }
    System.out.println("Basis of row space and of nullspace are orthogonal: ");
    System.out.println(RealUtil.dot(new RealMatrix(rowBasis).transpose(), 
        new RealMatrix(nullBasis)));
    }

  public static void projectionTest() {
    RealMatrix basis = new RealMatrix(new double[][]
        {{1, 0, 1}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}});
    RealVector v = new RealVector(new double[]{1, 2, 1, 2});
    System.out.println(v.vectorProjectionRaw(basis));
    System.out.println(v.vectorProjection(basis));
    System.out.println(RealUtil.dot(basis.toProjectionMatrix(), v));
  }

  public static void fittingTest() {
    double[]x = new double[]{-2, -1, 0, 1, 2};
    double[]y = new double[]{18, 8.1, 3.8, 3.1, 6.1};
    System.out.println(RealUtil.fittedExpression(RealUtil.leastSquare(x, y)));
    System.out.println(RealUtil.fittedExpression(RealUtil.polyFit(x, y, 2)));
  }

  public static void qrTest() {
    RealMatrix A = new RealMatrix(new double[][]{{1, 2, 3}, {-1, 0, -3}, {0, -2, 3}});
    RealQr factors = new RealQr(A);
    System.out.println(factors.getQ());
    System.out.println(factors.getR());
    System.out.println(RealUtil.dot(factors.getQ(), factors.getR()));
    RealMatrix B = new RealMatrix(new double[][]{{3, 6, 5}, {4, 8, 12}});
    RealQr factorsB = new RealQr(B);
    System.out.println(factorsB.getQ());
    System.out.println(factorsB.getR());
    System.out.println(RealUtil.dot(factorsB.getQ(), factorsB.getR()));
    RealQr factorsBT = new RealQr(B.transpose());
    System.out.println(factorsBT.getQ());
    System.out.println(factorsBT.getR());
    System.out.println(RealUtil.dot(factorsBT.getQ(), factorsBT.getR()));
    System.out.println(factorsBT.getQFull());
    System.out.println(factorsBT.getQFull().isOrthonormal());
    System.out.println(factorsBT.getRFull());
    System.out.println(RealUtil.dot(factorsBT.getQFull(), factorsBT.getRFull()));
  }

  public static void qrSolveTest() {
    RealMatrix A = new RealMatrix
        (new double[][]{{1, 1, 1}, {1, 2, 3}, {1, -1, 1}, {2, -2, 2}});
    RealVector b = new RealVector(new double[]{6, 14, 2, 5});
    RealMatrix A2 = new RealMatrix(
      new double[][]{{1, 1, 1, 1}, {1, 2, 3, 4}, {2, 2, 2, 2}}
    );
    RealVector b2 = new RealVector(new double[]{10, 30, 20});
    System.out.println("QR start here: ");
    System.out.println(RealUtil.qrSolve(uniqueA, uniqueb));
    System.out.println(RealUtil.qrSolve(A, b));
    System.out.println(RealUtil.qrSolve(A2, b2));
  }

  public static void eigenTest() {
    RealSquare A = new RealSquare(new double[][]{{1, 0, 1}, {0, 1, -1}, {1, 1, 2}});
    RealSquare B = new RealSquare(new double[][]{{16, 2, 17}, {0, -43, 14}, {0, 0, 16}});
    RealSquare C = new RealSquare(new double[][]{{3, 2, 4}, {2, 0, 2}, {4, 2, 3}});
    System.out.println(Arrays.toString(A.eigenValue()));
    System.out.println(Arrays.toString(B.eigenValue()));
    System.out.println(Arrays.toString(C.eigenValue()));
    RealVector[] eigenVectorsA = A.eigenVector(1);
    for (int i = 0; i < eigenVectorsA.length; i++) {
      System.out.println(eigenVectorsA[i]);
    }
    RealVector[] eigenVectorsC = C.eigenVector(C.eigenValue()[0]);
    for (int i = 0; i < eigenVectorsC.length; i++) {
      System.out.println(eigenVectorsC[i]);
    }
  }

  public static void subMatrixTest() {
    RealMatrix A = new RealMatrix(new double[][]{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}});
    System.out.println(A.subMatrix(1, 2));
    System.out.println(A.subMatrix(0, 0));
    System.out.println(A.subMatrix(2, 3));
  }

  public static void houseHolderTest() {
    RealMatrix A =  new RealMatrix(new double[][]{{12, -51, 4}, {6, 167, -68}, {-4, 24, -41}});
    RealQr factors = new RealQr(A, "HH");
    System.out.println(factors.getQ());
    System.out.println(factors.getR());
    System.out.println(RealUtil.dot(factors.getQ(), factors.getR()));
    if (RealUtil.matrixCompare(A, RealUtil.dot(factors.getQ(), factors.getR()))) {
      System.out.println("Success!");
    }
    if (factors.getQ().isOrthonormal()) {
      System.out.println("Greater Success!");
    }
  }

  public static void cramerTest() {
    RealSquare A = new RealSquare
        (new double[][]{{1, 1, 1}, {1, 2, 3}, {1, 3, 6}});
    RealVector b = new RealVector(new double[]{5, 7, 11});
    System.out.println(RealUtil.cramerSolve(A, b));
    System.out.println(A.inversePainful());
  }

  public static void detTest() {
    // RealMatrix D = new RealMatrix(new double[][]{{2, -1, 0},{-1, 2, -1},{0, -1, 2}});
    RealSquare D = new RealSquare(new double[][]
        {{2, 1, 5, 2}, {4, 0, -3, -1}, {6, -2, 0, 8}, {1, 0, -1, 3}});
    System.out.println(D.det());
    System.out.println(D.detPainful());
  }

  public static void main(String[]args)
  {
/*  */
    // luSolveTest();
    // System.out.println("To compare...");
    // plduTest();
    // basisTest();
    // linSolveGaussianTest();
    // plduSolveTest();
    // rrefTest();
    // projectionTest();
    // fittingTest();
    qrTest();
    // qrSolveTest();
    // eigenTest();
    // subMatrixTest();
    // houseHolderTest();
    // detTest();
    // cramerTest();
    // System.out.println(uniqueA.inverseGaussian());
  }
}