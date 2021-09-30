package linalgexercise;

/**
 * Stores the result of a PLDU decomposition of a non-singular square matrix A: PA = LDU,
 * and returns different parts of the result.
 * 
 * <p> Also contains the static methods {@link #forwardSub(RealMatrix, RealVector, boolean)
 * forwardSub()} and {@link #backSub(RealMatrix, RealVector, boolean) backSub()}
 * that can be called from elsewhere.
 * 
 * @author Zheng Ma
 * @version %I%, %G%
 * @since 1.0
 */
public class RealPldu {
  private RealSquare matrixLu;
  protected int dim;
  private int[] p;
  protected boolean isOddParity = false;
  protected boolean isSingular = false;

  /**
   * Class constructor.  Creat the new object by specifying {@code matrixLu} and {@code p[]}
   *    explicitly.
   * @param lu {@code RealSquare} compactly stores the result of LU decomposition
   * @param p {@code int[]} compactly stores the matrix P
   * @param isOddParity {@code boolean} {@code true} iff P is odd
   */
  public RealPldu(RealSquare lu, int[] p, boolean isOddParity) {
    this.matrixLu = lu;
    this.dim = lu.size;
    this.p = p;
    this.isOddParity = isOddParity;
  }

  /**
   * Class constructor.  Perform an PLDU decomposition on {@code A}, and store the matrices
   * L and DU compactly in {@code matrixLu}, and store the matrix P compactly in the array {@code p}
   * @param A {@code RealSquare} the matrix to be decomposed
   */
  public RealPldu(RealSquare A) {
    int tPivot;
    int tSwap;
    double l;

    dim = A.nRow;
    matrixLu = A.copy(); // Make a copy so that modifying entries won't mess up the original matrix
    p = new int[dim];
    for (int i = 0; i < dim; i++) {
      p[i] = i; // Init the permutation matrix as an identity matrix
    }

    for (int pivot = 0; pivot < dim; pivot++) {
      tPivot = findNewPivot(pivot);
      if (Math.abs(matrixLu.mat[tPivot][pivot]) < RealUtil.TOL) {
        System.out.println("Only non-singular matrix can have a PLDU decomposition!");
        isSingular = true;
      } else {
        if (tPivot != pivot) {
          matrixLu.rowSwap(pivot, tPivot); 
          tSwap = p[tPivot];
          p[tPivot] = p[pivot];
          p[pivot] = tSwap;
          isOddParity = !isOddParity; // Every swap change tha parity of the permutation
        }
        // Below is pretty much the same as RealSquare.LUDecomposition(),
        // except storing all the entries compactly in one matrix.
        for (int row = pivot + 1; row < dim; row++) {
          l = matrixLu.mat[row][pivot] / matrixLu.mat[pivot][pivot];
          for (int col = pivot +1 ; col < dim; col++) {
            matrixLu.mat[row][col] -= matrixLu.mat[pivot][col] * l;
          }
          matrixLu.mat[row][pivot] = l;
        }
      }
    }
  }

  /**
   * Returns the largest element in the current column so that it can be swapped upward
   * and be used as a pivot
   * @param pRow {@code int} the current row
   * @return {@code int} the index of the row to be swapped up
   */
  private int findNewPivot(int pRow) {
    double max = Math.abs(matrixLu.mat[pRow][pRow]);
    int index = pRow;
    for (int row = pRow + 1; row < dim; row++) {
      if (Math.abs(matrixLu.mat[row][pRow]) > max) {
        max = Math.abs(matrixLu.mat[row][pRow]);
        index = row;
      }
    }
    return index;
  }

  /**
   * Returns the matrix L in the PLDU decomposition
   * @return {@code RealSquare} lower-triangular matrix whose diagonal elements are all 1,
   *         basically the lower-triangular (bottom-left) part of the compactly stored
   *         {@code matrixLu} with diagonal elements changed to 1
   */
  public RealSquare getL() {
    RealSquare matrixL = RealUtil.eye(dim); // Fill the bottom-left part of an identity matrix
    for (int row = 1; row < dim; row++) { // Start from 1 so as not to mess up the diagonal
      for (int col = 0; col < row; col++) {
        matrixL.mat[row][col] = this.matrixLu.mat[row][col];
      }
    }
    return matrixL;
  }

  /**
   * Returns the matrix DU
   * @return {@code RealSquare} an upper-triangular matrix, basically the upper-triangular
   *         (top-right) part of the compactly stored {@code matrixLu}.
   */
  public RealSquare getDU() {
    RealSquare matrixDu = new RealSquare(dim);
    for (int row = 0; row < dim; row++) {
      for (int col = row; col < dim; col++) {
        matrixDu.mat[row][col] = matrixLu.mat[row][col];
      }
    }
    return matrixDu;
  }

  /**
   * Returns the matrix U in the PLDU decompostion
   * @return {@code RealSquare} upper-triangular matrix, all of whose diagonal elements are 1
   */
  public RealSquare getU() {
    RealSquare matrixDu = getDU();
    for (int row = 0; row < dim; row++) {
      matrixDu.rowMulti(row, 1.0 / matrixDu.mat[row][row]);
    }
    return matrixDu;
  }

  /**
   * Returns an array that contains the diagonal elements of D in top-left to bottom-right order
   * @return {@code double[]} the diagonal elements of D in the PLDU decomposition
   */
  public double[] getDCompact() {
    double[] d = new double[dim];
    for (int i = 0; i < dim; i++) {
      d[i] = this.matrixLu.mat[i][i]; 
    }
    return d;
  }

  /**
   * Returns the diagonal matrix D in the PLDU decomposition
   * @return {@code RealMatrix} the diagonal matrix D in the PLDU decomposition
   */
  public RealSquare getD() {
    return RealUtil.diag(getDCompact());
  }

  /**
   * Returns a integer array that summarizes the permutation matrix P
   * @return {@code int[]} {@code p}. {@code p[i] == n} means {@code P[i][n]} is 1, or the ith row of 
   *         PA is the nth row of A.
   */
  public int[] getPCompact() {
    return p;
  }

  /**
   * Returns the permuation matrix P
   * @return {@code RealSquare} P
   */
  public RealSquare getP() {
    return RealUtil.permutationMatrix(p);
  }

  /**
   * Returns det(A) from det(D), changing sign if necessary.
   * @return {@code double} det(A)
   */
  public double det() {
    double product = 1.0;
    for (int i = 0; i < dim; i++) {
      product *= this.matrixLu.mat[i][i];
    }
    if (isOddParity) {
      product *= -1;
    }
    return product;
  }

  /**
   * Returns the matrix LDU by operating on {@code this.matrixLU} <i>only</i>,
   *  without retrieving L and DU seperately.
   * @return {@code RealMatrix} LDU
   */
  public RealSquare getLDU() {
    int bound;
    RealSquare matrixLdu = new RealSquare(dim);
    for (int row = 0; row < dim; row++) {
      for (int col = 0; col < dim; col++) {
        bound = Math.min(row - 1, col);
        for (int i = 0; i <= bound; i++) {
          matrixLdu.mat[row][col] += matrixLu.mat[row][i] * matrixLu.mat[i][col];
        }
        if (row <= col) {
          matrixLdu.mat[row][col] += matrixLu.mat[row][col];
        }
      }
    }
    return matrixLdu;
  }

  /**
   * Returns A = P^{-1}LDU by swaping the rows of LDU around
   * @return {@code RealMatrix} the matrix A that was initially decomposed
   */
  public RealSquare restore() {
    RealSquare matrixLdu = getLDU();
    RealSquare restoredA = new RealSquare(dim);
    for (int row = 0; row < dim; row++) {
      restoredA.setRow(p[row], matrixLdu.getRow(row));
    }
    return restoredA;
  }

  /**
   * Solves a lower-triangular linear system Lx = b by forward substitution.
   * @param L {@code RealMatrix} lower-triangular matrix L, zeroes are NOT ALLOWED on the diagonal.
   * @param b {@code RealVector} the right-hand vector
   * @param normalized {@code boolean} if {@code true}, the actual diagonal of {@code L}
   *                   is ignored, assumed to be all 1's.
   * @return {@code RealVector} the solution vector
   */
  public static RealVector forwardSub(RealMatrix L, RealVector b, boolean normalized) {
    if (!RealUtil.dimCheck(L.nRow, b.dim)) {
      return null;
    }
    RealVector x = new RealVector(L.nCol);
    int range = Math.min(L.nRow, L.nCol);
    for (int i = 0; i < range; i++) {
      x.vec[i] = b.vec[i];
    }
    for (int col = 0; col < range; col++) {
      if (!normalized) {
        x.vec[col] /= L.mat[col][col];
      }
      for (int row = col + 1; row < range; row++) {
        x.vec[row] -= L.mat[row][col] * x.vec[col];
      }
    }
    return x;
  }

  /**
   * Solves an upper-triangular linear system Ux = b using back substitution
   * @param U {@code RealMatrix} upper triangular matrix, zeroes are NOT ALLOWED on the diagonal.
   *          (if it weren't, only the upper-right part is used)
   * @param b {@code RealVector} the right-hand vector
   * @param normalized {@code boolean} if {@code true}, the actual diagonal of {@code U}
   *                   is ignored, assumed to be all 1's.
   * @return {@code RealVector} the solution vector
   */
  public static RealVector backSub(RealMatrix U, RealVector b, boolean normalized) {
    if (!RealUtil.dimCheck(U.nRow, b.dim)) {
      return null;
    }
    RealVector x = new RealVector(U.nCol);
    int range = Math.min(U.nRow, U.nCol);
    for (int i = 0; i < range; i++) {
      x.vec[i] = b.vec[i];
    }
    for (int col = range - 1; col > -1; col--) {
      if (!normalized) {
        x.vec[col] /= U.mat[col][col];
      }
      for (int row = col - 1; row > -1; row--) {
        x.vec[row] -= U.mat[row][col] * x.vec[col];
      }
    }
    return x;
  }

  /**
   * Solve the consistent linear system Ax = b by solving Lc = P^{-1} b via forward substitution, 
   *    followed by solving DUx = c via back substitution.
   * @param b {@code RealVector} the right-hand vector
   * @return {@code RealVector} the solution vector
   */
  public RealVector solve(RealVector b) {
    RealVector bPermutated = new RealVector(dim);
    for (int i = 0; i < dim; i++) {
      bPermutated.vec[p[i]] = b.vec[i]; // equivalent to multiply b by P^{-1}
    }
    // Because forwardSub() does NOT care about the upper-right elements,
    // there is no need to extract L using getL().
    // Put normalized = true so that the diagonal elements will be treated as 1
    RealVector c = forwardSub(matrixLu, bPermutated, true);
    // Put normalized = false so that the diagonal elements of matrixLU will be used.
    RealVector x = backSub(matrixLu, c, false);
    return x;
  }

  /**
   * Returns A^{-1} by calculating its ith row by solve(e_i), where e_i is the ith column
   * of the identity matrix.
   * @return {@code RealSquare} A^{-1}
   * @deprecated slightly less efficient than {@link #inverse()}
   */
  @Deprecated
  public RealSquare inverseRaw() {
    RealSquare identity = RealUtil.eye(dim);
    RealSquare result = new RealSquare(dim);
    for (int i = 0; i < dim; i++) {
      result.setCol(i, solve(identity.getCol(i)));
    }
    return result;
  }

  /**
   * Returns A^{-1} by largely the same idea as {@link #inverseRaw()}, but improved efficiency
   * by a small margin via exploiting the structure of L, U and I.
   * @return {@code RealMatrix} A^{-1}
   */
  public RealSquare inverse() {
    RealSquare identity = RealUtil.eye(dim);
    RealSquare result = new RealSquare(dim);

    for (int i = 0; i < dim; i++) {
      for (int col = i; col < dim; col++) { // No need to start from the first column
        for (int row = col + 1; row < dim; row++) {
          identity.mat[row][i] -= matrixLu.mat[row][col] * identity.mat[col][i];
        }
      }
      for (int col = dim - 1; col > -1; col--) {
        identity.mat[col][i] /= matrixLu.mat[col][col];
        for (int row = col - 1; row > -1; row--) {
          identity.mat[row][i] -= matrixLu.mat[row][col] * identity.mat[col][i];
        }
      }
      result.setCol(p[i], identity.getCol(i));
    }
    return result;
  }
}