package linalgexercise;

// import java.io.*;
// import java.util.Arrays;
import java.lang.Math;

/**
 * Stores a square matrix, as well as many methods specific to a square matrix.
 * 
 * @author Zheng Ma
 * @version %I%, %G%
 * @since 1.0
 */
public class RealSquare extends RealMatrix {
  protected int size;

  /**
   * Class constructor.  Initialize a matrix with the entries being specified with a 
   *    2d array.  
   * @param init {@code double[][]} 2d array that carries the initial values of the entries,
   *             has to be square.
   */
  public RealSquare(double[][] init) {
    super(init);
    if (!RealUtil.dimCheck(init.length, init[0].length)) {
      mat = null;
      size = 0;
    } else {
      size = init.length;
    }
  }

  /**
   * Class constructor.  Initialize an n × n square matrix will all the entries being zeros.
   * @param n {@code int} dimension of the matrix
   */
  public RealSquare(int n) {
    super(n, n);
    size = n;
  }

  /**
   * Make a deep copy of {@code this}.
   * @return {@code RealSquare} a new square matrix whose entris are identical to {@code this}.
   */
  @Override
  public RealSquare copy() {
    return new RealSquare(this.mat);
  }

  /**
   * Returns a submatrix of {@code this} by deleting one row and one column.
   * @param iRow {@code int} the index of the row to be deleted.
   * @param iCol {@code int} the index of the column to be deleted
   * @return {@code RealSquare} the (n-1) × (n-1) submatrix
   */
  @Override
  public RealSquare subMatrix(int iRow, int iCol) {
    return super.subMatrix(iRow, iCol).toSquare();
  }

  /**
   * Returns the determinant of {@code this} using a default method.
   * @return {@code this}^{-1}
   */
  public double det() {
    return new RealPldu(this).det();
  }

  /**
   * Equivalent to {@code RealUtil.detPainful(this)}
   * 
   * <p>The implementation isn't in this class, because it's implemented in an recursive manner.
   * @return {@code double} det({@code this})
   */
  @Deprecated
  public double detPainful() {
    return RealUtil.detPainful(this);
  }

  /**
   * Determines whether {@code this} is a singular matrix
   * @return {@code boolean} {@code true} iff {@code this} is singular
   */
  public boolean isSingular() {
    // As stated in "The Biggest Theorem in Introductory Linear Algebra",
    // This is equivalent to over 20 different statements.
    // But this is the easiest one to check with existing code.
    return new RealRref(this).rank != size;
  }

  /**
   * Determines whether {@code this} is a symmetric matrix.
   * @return {@code boolean} {@code true} iff {@code this} is symmetric.
   */
  public boolean isSymmetric() {
    return RealUtil.matrixCompare(this, this.transpose());
  }

  /**
   * Returns the inverse of {@code this} using a default method.
   * @return {@code this}^{-1}
   */
  public RealSquare inverse() {
    // return inverseGaussian();
    return new RealPldu(this).inverse();
  }

  /**
   * Returns the inverse of {@code this} via Gauss-Jordan elimination
   * @return {@code RealSquare} {@code this}^{-1}
   */
  public RealSquare inverseGaussian() {
    RealMatrix augmented = concatenate(RealUtil.eye(size), true);
    RealRref solved = new RealRref(augmented);
    // If this is a non-singular matrix, the left half should now be an identity matrix,...
    if (RealUtil.matrixCompare(solved.getBlock(0, 0, this.size, this.size), RealUtil.eye(size))) {
      // ... and the right half should now be this^{-1}.
      RealSquare inverted = new RealSquare(size);
      for(int col = 0; col < size; col++) {
        inverted.setCol(col, solved.getCol(size + col));
      }
      return inverted;
    } else {
      System.out.println("Singular matrix is not invertible!");
      return null;
    }
  }

  /**
   * Returns the result of LU decomposition of {@code this}
   * @return {@code RealSquare[]} the {@code [0]} and {@code [1]} elements are L and U,
   *         respectively.  All diagonal elements of L are 1.
   * @deprecated This is for practicing LU decomposition only.  Use the corresponding methods
   *             in {@link RealPldu} instead once the latter is successfully implemented.
   */
  @Deprecated
  public RealSquare[] luDecomposition() {
    double l;

    if (isSingular()) {
      System.out.println("Only deals with non-singular matrices for now!");
      return null;
    }
    RealSquare matrixL = RealUtil.eye(size); // The diagonal entries of L are 1's
    RealSquare matrixU = this.copy(); // U is the REF of the original matrix
    // Basically a Gaussian elimination with the multipliers recorded in L
    for (int pivot = 0; pivot < size; pivot++) {
      for (int row = pivot + 1; row < size; row++) {
        l = matrixU.mat[row][pivot] / matrixU.mat[pivot][pivot];
        matrixU.rowSum(row, pivot, -l);
        matrixL.mat[row][pivot] = l;
      }
    }
    RealSquare[] results = new RealSquare[2];
    results[0] = matrixL;
    results[1] = matrixU;
    return results;
  }

  /**
   * Returns the result of LDU Decomposition of {@code this}
   * @return {@code RealSquare[]} the {@code [0]}, {@code [1]} and {@code [2]} elements
   *         are the matrices L, D and U, respectively. 
   * @deprecated This is for practicing LU decomposition only.  Use the corresponding methods
   *             in {@link RealPldu} instead once the latter is successfully implemented. 
   */
  @Deprecated
  public RealSquare[] lduDecomposition()
  {
    double d;
    double[] diagonal = new double[size];
    RealSquare[] lu = luDecomposition();
    RealSquare matrixL = lu[0];
    RealSquare matrixU = lu[1]; // Again, this is a clumsy method to return two objects at once
    for (int i = 0; i < size; i++) {
      d = matrixU.get(i, i);
      matrixU.rowMulti(i, 1.0/d);
      diagonal[i] = d;
    }
    RealSquare matrixD = RealUtil.diag(diagonal);
    return new RealSquare[]{matrixL, matrixD, matrixU};
  }

  /**
   * Returns the inverse of {@code this} via LU decompostion
   * @return {@code this}^{-1}
   * @deprecated This is for practicing LU decomposition.  Use {@link RealPldu#inverse()} instead
   * once the latter is successfully implemented.
   */
  @Deprecated
  public RealSquare luInverse() {
    // An n×n identity matrix, whose ith column vectors is e_i 
    RealSquare result = RealUtil.eye(size); 
    for (int i = 0; i < size; i++) {
      // The ith column of A^{-1} is the solution vector of Ax = e_i
      result.setCol(i, RealUtil.luSolve(this, result.getCol(i)));
    }
    return result;
  }

  /**
   * Returns the eigenvalues of {@code this} calculated from QR decomposition.
   * <p> This did NOT demonstrate some more advanced numerical techniques, 
   *    such as "A_n - cI = QR -> A_{n+1} = RQ + cI" or converting to Hessenberg matrix
   * @return {@code double[]} each element is an eigenvalue of {@code this}.  If an eigenvalue
   *         has a algebraic multiplicity >1, it will be stored more than once.
   */
  public double[] eigenValue() {
    RealSquare A = this.copy();
    // Compare the diagonal between two iterations, until convergence.
    double[] diag = new double[size];
    double[] diagnew = new double[size];
    final int MAX_ITER = 1024; // If MAX-ITER is reached, the matrix probably can't be diagonized.
    // In practice, the diagonal converges within 10 iterations in most cases.
    for (int i = 0; i < MAX_ITER; i++) {
      RealQr factors = new RealQr(A, "GS"); // A_n = QR
      A = RealUtil.dot(factors.getR(), factors.getQ()).toSquare(); //A_{n+1} = RQ
      diagnew = A.getDiagonal();
      if (RealUtil.vectorCompare(new RealVector(diag), new RealVector(diagnew))) {
        System.out.println("Iterations: " + i);
        break;
      } else {
        diag = diagnew;
      }
    }
    return diagnew;
  }
  
  /**
   * Returns the eigenvector(s) of {@code this} that correspond to the eigenvalue {@code lambda}
   * @param lambda {@code double} an eigenvalue of {@code this} (prompt if it weren't)
   * @return {@code RealVector} eigenvector(s) corresponding to {@code lambda} if the latter
   *         is indeed an eigenvalue of {@code this}, otherwise prompt and return {@code null}.
   */
  public RealVector[] eigenVector(double lambda) {
    // A = this - λI
    RealSquare A = RealUtil.sub(this, RealUtil.eye(size).scalarMulti(lambda)).toSquare();
    // The eigenvectors should be the nullspace of A
    RealVector[] eigenVectors = A.nullBasis();
    // If the nullspace is {0}, you'd better double-check that lambda is indeed an eigenvalue.
    if (eigenVectors.length == 0) {
      System.out.println(lambda + " might not be an eigenvalue of this matrix.");
      return null;
    } else {
      return eigenVectors;
    }
    // This method is just imitating the procedure taught in the class.
    // There are better numerical algorithms for calculate eigenvectors efficiently and precisely.
  }
  
  /**
   * Calculate the inverse of {@code this} using the determinant and cofactors of each element.
   * @return {@code RealSquare} the inverse of {@code this}
   * @deprecated Much less efficient than any other methods (RREF, LU, PLDU, QR, etc.)
   */
  @Deprecated
  public RealSquare inversePainful() {
    double detA = this.detPainful();
    if (Math.abs(detA) < RealUtil.TOL) {
      System.out.println("Singular matrix is not invertible!");
      return null;
    }
    RealSquare inverted = new RealSquare(size);
    for (int row = 0; row < size; row++) {
      for (int col = 0; col < size; col++) {
        inverted.mat[row][col] = Math.pow(-1, row + col)
            * this.subMatrix(col, row).detPainful() / detA;
      }
    }
    return inverted;
  }
}