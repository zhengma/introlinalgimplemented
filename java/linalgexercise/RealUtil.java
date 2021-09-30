package linalgexercise;

// import java.io.*;
// import java.util.Arrays;
import java.lang.Math;

public class RealUtil
{
  public final static double TOL = 1e-8;

  public RealUtil(){}

  /**
   * check whether two integers equal to each other.  Used for all the dimension-agreement
   *    checking in the library.
   * @param d1 {@code int} one dimension
   * @param d2 {@code int} the other dimension
   * @return {@code true} if {@code d1} and {@code d2} agree.  Otherwise print the two dimensions that
   *         disagree and return {@code false}.
   */
  public static boolean dimCheck(int d1, int d2) {
    if (d1 != d2) {
      System.out.println("Dimensions that should've matched: " + d1 + " and " + d2);
    }
    return (d1 == d2);
  }

  /**
   * Add two vectors
   * @param u {@code RealVector} n-dimensional vector
   * @param v {@code RealVector} n-dimensional vector
   * @return {@code RealVector} n-dimensional vector u + v
   */
  public static RealVector add(RealVector u, RealVector v) {
    if (!dimCheck(u.dim, v.dim)) {
      return null;
    }
    RealVector sum = new RealVector(u.dim);
    for (int i = 0; i < u.dim; i++) {
      sum.vec[i] = u.vec[i] + v.vec[i];
    }
    return sum;
  }

  /**
   * Return the sum of two matrices
   * @param A {@code RealMatrix} m × n matrix
   * @param B {@code RealMatrix} m × n matrix
   * @return {@code RealMatrix} m × n matrix A + B
   */
  public static RealMatrix add(RealMatrix A, RealMatrix B) {
    if (!dimCheck(A.nRow, B.nRow)) {
      return null;
    }
    if (!dimCheck(A.nCol, B.nCol)) {
      return null;
    }
    RealMatrix sum = new RealMatrix(A.nRow, A.nCol);
    // You can use row vectors instead, or even two layers of nested for-loop
    for (int col = 0; col < A.nCol; col++) { 
      sum.setCol(col, add(A.getCol(col), B.getCol(col)));
    }
    return sum;
  }

  /**
   * Return the difference between two vectors.
   * @param u {@code RealVector} n-dimensional vector
   * @param v {@code RealVector} n-dimensional vector
   * @return {@code RealVector} n-dimensional vector u - v
   */
  public static RealVector sub(RealVector u, RealVector v) {
    int n = u.dim;
    if (!dimCheck(n, v.dim)) {
      return null;
    }
    RealVector diff = new RealVector(n);
    for (int i = 0; i < n; i++) {
      diff.vec[i] = u.vec[i] - v.vec[i];
    }
    return diff;
  }

  /**
   * Return the difference of two matrices
   * @param A {@code RealMatrix} m × n matrix
   * @param B {@code RealMatrix} m × n matrix
   * @return {@code RealMatrix} m × n matrix A - B
   */
  public static RealMatrix sub(RealMatrix A, RealMatrix B) {
    if (!dimCheck(A.nRow, B.nRow)) {
      return null;
    }
    if (!dimCheck(A.nCol, B.nCol)) {
      return null;
    }
    
    RealMatrix diff = new RealMatrix(A.nRow, A.nCol);
    for (int col = 0; col < A.nCol; col++) {
      diff.setCol(col, sub(A.getCol(col), B.getCol(col)));
    }
    return diff;
  }

  /**
   * Returns the dot product of two vectors
   * @param u {@code RealVector} n-dimensional vector
   * @param v {@code RealVector} n-dimensional vector
   * @return {@code double} u·v
   */
  public static double dot(RealVector u, RealVector v)
  {
    if (!dimCheck(u.dim, v.dim)) {
      return 0.0;
    }
    double sum = 0.0;
    for (int i = 0; i < u.dim; i++) {
      sum += u.vec[i] * v.vec[i];
    }
    return sum;
  }

  /**
   * Returns the product of two matrices.
   * @param A {@code RealMatrix} m × k matrix
   * @param B {@code RealMatrix} k × n matrix
   * @return {@code RealMatrix} m × n matrix A·B
   */
  public static RealMatrix dot(RealMatrix A, RealMatrix B) {
    if (!dimCheck(A.nCol, B.nRow)) {
      return null;
    }

    RealMatrix product = new RealMatrix(A.nRow, B.nCol);
    for (int row = 0; row < A.nRow; row++) {
      for (int col = 0; col < B.nCol; col++) {
        product.mat[row][col] = dot(A.getRow(row), B.getCol(col));
      }
    }          
    return product;
  }

  /**
   * Returns the product Ab
   * @param A {@code RealMatrix} m × n matrix 
   * @param b {@code RealVector} n-dimensional vector
   * @return {@code RealVector} m-dimensional vector A·b
   */
  public static RealVector dot(RealMatrix A, RealVector b) {
    if (!dimCheck(A.nCol, b.dim)) {
      return null;
    }

    RealVector product = new RealVector(A.nRow);
    for(int i = 0; i < A.nRow; i++)
      product.vec[i] = dot(A.getRow(i), b);
    return product;
    // An even easier way is just "return dot(A, b.toMatrix(true))".
  }

  /**
   * Returns the product of a^T B
   * @param a {@code RealVector} m-dimensional vector
   * @param B {@code RealMatrix} m × n matrix
   * @return {@code RealVector} a (as a row vector) times B, a n-dimensional vector
   */
  public static RealVector outer(RealVector a, RealMatrix B) {
    if (!dimCheck(a.dim, B.nRow)) {
      return null;
    }
    RealVector product = new RealVector(B.nCol);
    for(int i = 0; i < B.nCol; i++) {
      product.vec[i] = dot(a, B.getCol(i));
    }
    return product;
  }

  /**
   * Returns the outer product of two vectors
   * @param u {@code RealVector} m-dimensional vector
   * @param v {@code RealVector} n-dimensional vector
   * @return {@code RealMatrix} u (as a column matrix) times v (as a row matrix),
   *         an m × n matrix.
   */
  public static RealMatrix outer(RealVector u, RealVector v) {
    RealMatrix product = new RealMatrix(u.dim, v.dim);
    for (int row = 0; row < u.dim; row++) {
      for (int col = 0; col < v.dim; col++) {
        product.mat[row][col] = u.vec[row] * v.vec[col];
      }
    }
    return product;
  }

/**
   * Return the cross product between two 3d vectors
   * @param u {@code RealVector} the first vector
   * @param v {@code RealVector} the second vector
   * @return {@code RealVector} u × v
   */
  public static RealVector cross(RealVector u, RealVector v) {
    if (!dimCheck(u.dim, 3) || !dimCheck(v.dim, 3)) {
      return null; // Conventional cross product only exists between 3d vectors 
    }
    double xComponent = u.vec[1] * v.vec[2] - u.vec[2] * v.vec[1];
    double yComponent = u.vec[2] * v.vec[0] - u.vec[0] * v.vec[2];
    double zComponent = u.vec[0] * v.vec[1] - u.vec[1] * v.vec[0]; // This is actually less code
    return new RealVector(new double[]{xComponent, yComponent, zComponent});
  }

  /**
   * Determine whether two vectors are the same within floating-point error
   * @param v1 {@code RealVector} one vector
   * @param v2 {@code RealVector} the other vector
   * @return {@code true} if the root-mean-square difference between the components of
   *         the two vectors is within TOL
   */
  public static boolean vectorCompare(RealVector v1, RealVector v2) {
    if (!RealUtil.dimCheck(v1.dim, v2.dim)) {
      return false;
    }
    RealVector diff = RealUtil.sub(v1, v2);
    return (diff.modulus() / v1.dim < RealUtil.TOL);
  }

  /**
   * Determine whether two matrices are the same withing floating-point error
   * @param m1 {@code RealMatrix} one matrix
   * @param m2 {@code RealMatrix} the other matrix
   * @return {@code true} if the two matrices are of the same dimensions,
   *         and root-mean-square difference between the corresponding
   *         entries of the two marices is within {@code TOL}.
   */
  public static boolean matrixCompare(RealMatrix m1, RealMatrix m2) {
    double totalMod = 0.0;
    
    if ((!RealUtil.dimCheck(m1.nRow, m2.nRow)) || 
        (!RealUtil.dimCheck(m1.nCol, m2.nCol))) {
      return false;
    }
    RealMatrix diff = RealUtil.sub(m1, m2);
    for (int i = 0; i < diff.nRow; i++) {
      totalMod += Math.pow(diff.getRow(i).modulus(), 2);
    }
    return (Math.sqrt(totalMod) / diff.nRow / diff.nCol < RealUtil.TOL);
  }

  /**
   * check if two vectors are orthogonal to each other.
   * @param u {@code RealVector} one vector
   * @param v {@code RealVector} the other vector
   * @return {@code true} iff u and v are orthogonal to each other within floating-point error.
   */
  public static boolean isOrthogonal(RealVector u, RealVector v) {
    return Math.abs(dot(u, v)) < TOL;
  }

  /**
   * Find the linear combination of a set of vectors
   * @param vectors {@code LinAlgMatrix} has n columns: v_0, v_1, ..., v_{n-1}
   * @param coefficients {@code double[n]} the coefficients c_0, ... c_{n-1}
   * @return {@code RealVector} c_0 v_0 + c_1 v_1 + ... + c_{n-1} v_{n-1}
   */
  public static RealVector linComb(RealMatrix vectors, double[] coefficients) {
    return dot(vectors, new RealVector(coefficients)); // Yes, it's this simple.
  }

  /**
   * Find the linear combination of a set of vectors
   * @param vectors {@code RealVectors[n]} vectors v_0, v_1, ..., v_{n-1}
   * @param coefficients {@code double[n]} the coefficients c_0, ... c_{n-1}
   * @return {@code RealVector} c_0 v_0 + c_1 v_1 + ... + c_{n-1} v_{n-1}
   */
  public static RealVector linComb(RealVector[] vectors, double[] coefficients) {
    RealMatrix vectorsMatrix = new RealMatrix(vectors);
    return linComb(vectorsMatrix, coefficients);
  }

  /**
   * Calculate the angle between two vectors.
   * @param u {@code RealVector} one of the vectors
   * @param v {@code RealVector} the other vector
   * @return {@code double} the angle between {@code u} and {@code v} in radian. between 0 and π
   */
  public static double angle(RealVector u, RealVector v) {
    return Math.acos(dot(u, v)/u.modulus()/v.modulus());
  }

  /**
   * Construct an identity matrix (similar to {@code Sympy.matrices.eye()})
   * @param n {@code int} the dimension of the matrix
   * @return {@code RealSquare} an n x n dimension matrix.
   */
  public static RealSquare eye(int n) {
    RealSquare idMatrix = new RealSquare(n);
    for (int i = 0; i < n; i++) {
      idMatrix.mat[i][i] = 1.0;
    }
    return idMatrix;
  }

  /**
   * Construct a diagonal matrix D with specified diagonal elements
   * @param diagonal {@code double[]} elements on the diagonal in top-left to bottom-right order
   * @return {@code RealSquare} that contains the diagonal matrix
   */
  public static RealSquare diag(double[] diagonal) {
    RealSquare matrixD = new RealSquare(diagonal.length);
    for (int i = 0; i < diagonal.length; i++) {
      matrixD.mat[i][i] = diagonal[i];
    }
    return matrixD;
  }

  /**
   * Construct a permutation matrix P based on the reordering specified.
   * @param order {@code int[n]} which contains a permutation of integers
   *              from {@code 0 } to {@code (n-1)}
   * @return a {@code RealSquare} P, where for all matrices A of appropriate dimension,
   *         the {@code i}th row of PA is the {@code p[i]}th row of A. 
   */
  public static RealSquare permutationMatrix(int[] order) {
    RealSquare matrixP = new RealSquare(order.length);
    for (int i = 0; i < order.length; i++) {
      matrixP.mat[i][order[i]] = 1.0;
    }
    return matrixP;
  }

  /**
   * Given a linear system Ax = b, print one of {@code caseNone}, {@code caseUnique}
   *    and {@code caseInf} as appropriate.
   * <p> If {@code caseUnique} is printed, follow it by printing the solution vector.
   * <p> If {@code caseInf} is printed, follow it by printing the general solution.
   *     You may use {@code caseInfParticular} and {@code caseInfNullSpace} to express the
   *     general solution.
   * <p >The task shall be accomplished using Gauss-Jordan elimination. 
   * @param A {@code RealMatrix} the coefficient matrix
   * @param b {@code RealVector} the right-hand vector
   * @return {@code RealVector} the solution vector <i>only when</i> the solution is unique.
   *         Otherwise return {@code null}
   */
  public static RealVector linSolveGaussian(RealMatrix A, RealVector b)
  {
    String caseNone = "The system Ax = b has no solution.";
    String caseUnique = "The system Ax = b has a unique solution: ";
    String caseInf = "The system Ax = b has infinitely many solutions.";
    String caseInfParticular = "But all solution vectors can be expressed as this: ";
    String caseInfNullSpace = "Plus any linear combination of the following vectors: ";

    if (!dimCheck(A.nRow, b.dim)) {
      return null;
    }
    RealMatrix augmented = new RealMatrix(A.nRow, A.nCol + 1);
    for (int j = 0; j < A.nCol; j++) {
      augmented.setCol(j, A.getCol(j));
    }
    augmented.setCol(A.nCol, b); // The augmented matrix is A|b
    RealRref augmentedReduced = new RealRref(augmented); // Find the RREF of A|b

    // Conclude that "the system has no solution" if in the RREF of A|b,
    // the pivot column of a row is the last column,
    // that is: the only non-zero entry of that row is the last one.
    if (augmentedReduced.pivots[A.nCol]) {
      System.out.println(caseNone);
      return null;
    } else if (augmentedReduced.rank == A.nCol) { // A|b have the same rank as A
      System.out.println(caseUnique);
      return augmentedReduced.getCol(A.nCol);
    } else {
      System.out.println(caseInf);
      // When the system has infinitely many solutions, the general solution is any particular
      // solution plus any vector int he nullspace of A.
      System.out.println(caseInfParticular);
      // One particular solution is constructed as below
      RealVector particular = new RealVector(A.nCol); // New a solution vector of the right size.
      int row = 0; // Start from the first row
      for (int col = 0; col < A.nCol; col++) {
        if (augmentedReduced.pivots[col]) {  // If it's a pivot column
          // Put the last element of the row in the corresponding entry of the solution
          particular.vec[col] = augmentedReduced.mat[row][A.nCol];
          row++;  // Move on to the next row
        }
      }
      System.out.println(particular);
      System.out.println(caseInfNullSpace);
      RealVector[] nullSpace = A.nullBasis();
      for (int i = 0; i < nullSpace.length; i++) {
        System.out.println(nullSpace[i]);
      }
      return null;
    }
  }

  /**
   * Solve a consistent and independent linear system Ax = b
   * @param A {@code RealSquare} the coefficient matrix.  Has to be non-singular,
   *          square, and its LU decompostion must not involve swapping of rows
   * @param b {@code RealVector} a vector of matching dimension
   * @return {@code RealVector} the solution vector x for Ax = b
   * @deprecated This is for familiarizing with LU decomposition only.
   * Use {@link RealPldu#solve(RealVector)} instead after successfully implementing the latter.
   */
  @Deprecated
  public static RealVector luSolve(RealSquare A, RealVector b) {
    if (!dimCheck(A.nRow, b.dim)) {
      return null;
    }

    int n = A.nRow; // Just use A.nRow instead is also okay
    // In Java, there isn't a truly convenient way to return more than one thing at once
    // You can think of your own workaround, such as defining a small class just for this.
    RealSquare[] lu = A.luDecomposition();
    RealSquare matrixL = lu[0];
    RealSquare matrixU = lu[1];

    RealVector c = b.copy();
    // This can serve the starting point of your implementation of RealPldu.forwardSub()
    for (int col = 0; col < n; col++) {
      for (int row = col + 1; row < n; row++) {
        c.vec[row] -= matrixL.mat[row][col] * c.vec[col];
      }
    }

    RealVector x = c.copy();
    // This can serve the starting point of your implementation of RealPldu.backSub()
    for (int col = n - 1; col > -1; col--) {
      x.vec[col] /= matrixU.mat[col][col];
      for (int row = col - 1; row > -1; row--) {
        x.vec[row] -= matrixU.mat[row][col] * x.vec[col];
      }
    }
    return x;
  }

  /**
   * Find the least-square line of best fit for a set of data points.
   * @param x {@code double[]} the x coordinates of the points
   * @param y {@code double[]} the matching y coordinates of the points. {@code x.length}
   *          and {@code y.length} have to match
   * @return {@code double[2]} array containing the intercept and the slope of the best fit line,
   *         in that order.
   * @see {@link RealMatrix#toProjectionMatrix()}, {@link RealVector#vectorProjection(RealMatrix)}
   */
  public static double[] leastSquare(double[] x, double[] y) {
    if(!dimCheck(x.length, y.length)) {
      return null;
    }
    RealMatrix basis = new RealMatrix(x.length, 2);
    for (int row = 0; row < x.length; row++) {
      basis.mat[row][0] = 1.0; // The first columns is filled with 1.0
      basis.mat[row][1] = x[row]; // The second column is filled with elements of x[]
    }
    RealSquare A = RealUtil.dot(basis.transpose(), basis).toSquare();
    RealVector b = RealUtil.dot(basis.transpose(), new RealVector(y));
    RealVector coefficients = new RealPldu(A).solve(b);
    return coefficients.vec;
  }

  /**
   * Find the polynomial that best fits (in a least-square sense) a set of data points
   *    (similar to {@code polyfit()} in Matlab)
   * <p> {@link #leastSquare()} will be equivalent to {@code polyFit(x, y, 1)}
   * @param x {@code double[]} the x coordinates of the points
   * @param y {@code double[]} the matching y coordinates of the points. {@code x.length}
   *          and {@code y.length} have to match
   * @param degree the degree of the polynomial
   * @return {@code double[degree + 1]} the coefficients in ascending order.
   * @see {@link RealMatrix#toProjectionMatrix()}, {@link RealVector#vectorProjection(RealMatrix)}
   */
  public static double[] polyFit(double[] x, double[] y, int degree) {
    if(!dimCheck(x.length, y.length)) {
      return null;
    }
    // Dont' forget the "+1", as a polynomial of degree n has (n+1) terms.
    RealMatrix basis = new RealMatrix(x.length, degree + 1);
    for (int row = 0; row < x.length; row++) {
      // You can directly fill the first column with 1.0,
      // but that's some extra code and won't save terribly much time. 
      for (int power = 0; power <= degree; power++) {
        basis.mat[row][power] = Math.pow(x[row], power);
      }
    }
    RealSquare A = RealUtil.dot(basis.transpose(), basis).toSquare();
    RealVector b = RealUtil.dot(basis.transpose(), new RealVector(y));
    RealVector coefficients = new RealPldu(A).solve(b);
    return coefficients.vec;
  }

  /**
   * Display a univariate polynomial y = P(x) in a more reader-friendly way.
   * @param coefficients {@code double} the (n+1) coefficients of the polynomial of degree n,
   * in <i>ascending</i> order.
   * @return {@code String} the polynomial to be displayed,
   *         such as "y = 2 - 3 x + 4 x^2".
   */
  public static String fittedExpression(double[] coefficients) {
    String expression = "y = " + coefficients[0];
    for (int power = 1; power < coefficients.length; power++) {
      // Feel free to add another "else if" to skip terms of zero coefficients
      if (coefficients[power] > 0.0) {
        expression += " + " + coefficients[power] + " x";
      } else {
        expression += " - " + Math.abs(coefficients[power]) + " x";
      }
      if (power > 1) {
        expression += "^" + power;
      }
    }
    return expression;
  }

  /**
   * Calculate the determinant of the matrix A using the "big" formula by adding up
   * n! terms. 
   * 
   * @param A {@code RealSquare} a square matrix
   * @return {@code double} det(A)
   * @deprecated This is a VERY inefficient algorithm.  In practice, use Gaussian elimination
   *             or LU decomposition instead.
   */
  @Deprecated
  public static double detPainful(RealSquare A) {
    if (A.size == 1) {  // When using recursion, be certain an "ending case" is set
      return A.mat[0][0];
    } else {
      double sum = 0.0;
      for (int col = 0; col < A.size; col++) {
        // Skip any zero entries that are encountered, this will save a LOT of time,
        // especially in a somewhat sparse matrix.
        if (Math.abs(A.mat[0][col]) > TOL) {
          // Implemented with recursion, this is why the method needs a parameter
          sum += Math.pow(-1, col) * A.mat[0][col] * detPainful(A.subMatrix(0, col));
        }
      }
      return sum;
    }
  }

  /**
   * Find the unique solution to the consistent linear system Ax = b using Cramer's Rule.
   * 
   * @param A {@code RealSquare} the coefficient matrix.  Has to be square and non-singular.
   * @param b {@code RealVector} the right-hand vector
   * @return {@code RealVector} the solution vector x of Ax = b
   * 
   * @deprecated This is BY FAR the LEAST efficient algorithm to solve a linear system!
   * In practice, either RREF or LU decomposition or even QR decompostion should be used.
   */
  @Deprecated
  public static RealVector cramerSolve(RealSquare A, RealVector b) {
    if (!dimCheck(A.size, b.dim)) {
      return null;
    }
    double detA = detPainful(A);
    RealVector solution = new RealVector(b.dim);
    if (Math.abs(detA) < TOL) {
      System.out.println("A is singular!");
      return null;
    }
    RealSquare B;
    for (int i = 0; i < A.size; i++) {
      B = A.copy(); // A little inefficient, only two columns have to be changed.
      B.setCol(i, b);
      solution.vec[i] = B.detPainful() / detA;
    }
    return solution;
  }

  /**
   * Returns the solution of the linear system QRx = b.
   * <p> If the system has infinitely many solutions, returns one particular solution.
   * <p> If the system has a unique solution, returns that solution.
   * <p> If the system has no solution, returns a solution that minimizes |Ax - b|.
   * @param A {@code RealMatrix} the coefficient matrix
   * @param b {@code RealVector} the right-hand vector
   * @return the solution vector of Ax = b
   */
  public static RealVector qrSolve(RealMatrix A, RealVector b) {
    String caseDep = "The system is dependent.  Returns one particular solution.";
    String caseInd = "Returns the solution vector that minimizes |Ax - b|.";

    if (!dimCheck(A.nRow, b.dim)) {
      return null;
    }

    if (A.isFullColRank()) {
      System.out.println(caseInd);
      RealQr decomposed = new RealQr(A);
      // b' = Q^T · b
      RealVector bNew = decomposed.getQ().transpose().dot(b);
      // Solve Rx = b' using back substitution (because R is upper triangular)
      return RealPldu.backSub(decomposed.getR(), bNew, false); 
    } else {
      System.out.println(caseDep);
      RealQr decomposed = new RealQr(A.transpose());
      // Solve R^T c = b using forward substitution (because R^T is lower triangular)
      RealVector c = RealPldu.forwardSub(decomposed.getR().transpose(), b, false);
      // x = Qc is the unique solution if the system is consistent,
      // Otherwise it minimizes error.
      return decomposed.getQ().dot(c);
    }
  }
}