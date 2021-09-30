package linalgexercise;

/**
 * Stores the QR decomposition of a matrix (obtained by either Gram-Schimidt or Householder
 *    algorithm), and related methods.
 *
 * @author Zheng Ma
 * @version %I%, %G%
 * @since 1.0
 */
public class RealQr {
  protected RealMatrix matrixQ;
  protected RealMatrix matrixR;
  protected int nRow;
  protected int nCol;

  /**
   * Class constructor.  Perform a QR decomposition, store the matrices in {@code this.matrixQ}
   *    and {@code this.matrixR}.
   * @param A {@code RealMatrix} m × n matrix
   * @param method Either "GS" for Gram-Schmidt or "HH" for Householder
   */
  public RealQr(RealMatrix A, String method) {
    nRow = A.nRow;
    nCol = A.nCol;
    int rank = A.rank();
    if (A.isFullRowRank()) {
      matrixQ = new RealSquare(nRow);
    } else {
      matrixQ = new RealMatrix(nRow, rank);
    }
    if (A.isFullColRank()) {
      matrixR = new RealSquare(nCol);
    } else {
      matrixR = new RealMatrix(rank, nCol);
    }
    
    switch(method) {
      case "GS":
        gramSchmidt(A);
        break;
      case "HH":
        houseHolder(A, matrixQ, matrixR);
        break;
      default:
        System.out.println("Only GS and HH are available for now.");
    }
  }

  public RealQr(RealMatrix A) {
    this(A, "GS");
    // this(A, "HH");
  }
  
  /**
   * QR decomposition of matrix A based on Gram-Schmidt algorithm.
   * 
   * <p> It stores an m × r matrix Q in {@code this.matrixQ}.
   * Q will be column orthogonal: Q^T Q is the r × r identity matrix.
   * If A has full row rank, Q will be an orthogonal, square matrix.
   * 
   * <p> It stores an r × n upper-triangular matrix R in {@code this.matrixR}.
   * @param A {@code RealMatrix} m × n matrix of rank r
   */
  private void gramSchmidt(RealMatrix A) {
    matrixQ.setCol(0, A.getCol(0).normalize());
    matrixR.mat[0][0] = A.getCol(0).modulus();

    double k;
    RealVector v;
    int colQ = 1; // Avoid adding a new column in Q if a column of A is LD to those to its left
    for (int col = 1; col < nCol; col++) {
      v = A.getCol(col);
      for (int i = 0; i < colQ; i++) { // Subtract the current column by its projection...
        // on the space spanned by the existing columns of Q.
        // The result of subtraction must be orthogonal to all existing columns of Q
        k = matrixQ.getCol(i).dot(v);
        matrixR.mat[i][col] = k; // Entries of R are scalar projections of the current column...
        // on the existing columns of Q
        v = RealUtil.sub(v, matrixQ.getCol(i).scalarMulti(k));
      }
      if (!v.isZero()) { // It happens only when v is LI to the columns of A to its left.
        matrixQ.setCol(colQ, v.normalize());
        matrixR.mat[colQ][col] = v.modulus();
        colQ++;
      }
    }
  }

  /**
   * QR decomposition of matrix A based on Householder Algorithm
   * @param A {@code RealMatrix}
   * @param Q {@code RealMatrix}
   * @param R {@code RealMatrix}
   */
  private void houseHolder(RealMatrix A, RealMatrix Q, RealMatrix R) {
    RealVector u = A.getCol(0);
    if (u.vec[0] > 0.0) {
      u.vec[0] -= u.modulus();
    } else {
      u.vec[0] += u.modulus();
    }
    RealVector v = u.normalize();
    Q.setBlock(0, 0, RealUtil.sub(RealUtil.eye(A.nRow), RealUtil.outer(v, v).scalarMulti(2)));
    R.setBlock(0, 0, RealUtil.dot(Q, A));
    if (A.nCol > 2) {
      RealMatrix subQ = new RealMatrix(Q.nRow - 1, Q.nCol - 1);
      RealMatrix subR = new RealMatrix(R.nRow - 1, R.nCol - 1);
      houseHolder(R.subMatrix(0, 0), subQ, subR);
      RealMatrix restofQ = new RealMatrix(Q.nRow, Q.nCol);
      restofQ.mat[0][0] = 1;
      for (int row = 0; row < A.nRow - 1; row++) {
        for (int col = 0; col < A.nCol - 1; col++) {
          R.mat[row + 1][col + 1] = subR.mat[row][col];
          restofQ.mat[row + 1][col + 1] = subQ.mat[row][col];
        }
      }
      Q.setBlock(0, 0, RealUtil.dot(Q, restofQ));
    }
  }

  /**
   * Returns the reduced orthogonal matrix Q
   * @return {@code RealMatrix} m × r column orthogonal matrix
   */
  public RealMatrix getQ() {
    return matrixQ;
  }

  /**
   * Returns the reduced upper triangular matrix R
   * @return {@code RealMatrix} r × n upper triangular matrix
   */
  public RealMatrix getR() {
    return matrixR;
  }

  /**
   * Returns the full (not reduced) orthogonal matrix Q
   * @return {@code RealSquare} m×m orthogonal matrix.  Its column vectors form an orthonormal basis
   *         of R^m.  The first r vectors is the reduced Q that form a basis of col(A),
   *         the remaining (m-r) vectors form a basis of the orthogonal complement of col(A).
   */
  public RealSquare getQFull() {
    if (matrixQ.nRow == matrixQ.nCol) {
      return matrixQ.toSquare();
    } else {
      // The extra columns of full Q form a basis of the complement of col(A)
      RealMatrix complement = new RealMatrix(matrixQ.transpose().nullBasis());
      RealMatrix extraQ = new RealQr(complement).getQ(); // Make the extra vectors orthonormal,
      // each extra vector is already orthogonal to all the column vectors of reduced Q
      RealSquare fullQ = new RealSquare(nRow); // m × m square matrix
      fullQ.setBlock(0, 0, matrixQ);  // Original Q goes to the left
      fullQ.setBlock(0, matrixQ.nCol, extraQ); // Extra columns goes to the right
      return fullQ;
    }
  }

  /**
   * Returns a full (not reduced) upper triangular matrix R
   * @return {@code RealMatrix} m × n upper triangular matrix
   */
  public RealMatrix getRFull() {
    // Full R has the same shape as the original matrix
    RealMatrix fullR = new RealMatrix(nRow, nCol);
    // The extra bottom rows, if there is any, are all zeros,
    //    because this will be multiplied to the extra columns in the full Q
    //    that forms a basis for the orthogonal complement of col(A).
    fullR.setBlock(0, 0, matrixR);
    return fullR;
  }

  /**
   * Returns the matrix A that was fed into the constructor by multiplying Q by R
   * @return {@code RealMatrix} A = QR
   */
  public RealMatrix restore() {
    return RealUtil.dot(matrixQ, matrixR);
  }
}
