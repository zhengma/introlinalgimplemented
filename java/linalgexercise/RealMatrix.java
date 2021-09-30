package linalgexercise;

import java.util.Arrays;
// import java.util.ArrayList;
import java.lang.Math;

/**
 * Stores a matrix, relevant informations and many relevant methods.
 * 
 * @author Zheng Ma
 * @version %I%, %G%
 * @since 1.0
 */
public class RealMatrix {
  protected double[][] mat;
  protected int nRow;
  protected int nCol;
  protected RealRref rref;

  /**
   * Class constructor.  Initialize the matrix by specifying the entries with a 2d array
   * @param init {@code double[][]} the initial values of the entries
   */
  public RealMatrix(double[][] init) {
    nRow = init.length;
    nCol = init[0].length;
    mat = new double[nRow][nCol];
    for (int row = 0; row < nRow; row++) {
      for (int col = 0; col < nCol; col++) {
        mat[row][col] = init[row][col];
      }
    }
  }

  /**
   * Class constructor.  Initialize a zero matrix.
   * @param row {@code int} number of rows
   * @param col {@code int} number of columns
   */
  public RealMatrix(int row, int col) {
    this.nRow = row;
    this.nCol = col;
    mat = new double[nRow][nCol];
  }

  /**
   * Class constructor.  Arrange an array of column vectors of identical dimension
   *    from left to right into a matrix
   * @param columns {@code RealVector[]} an array of column vectors
   */
  public RealMatrix(RealVector[] columns) {
    nRow = columns[0].dim;
    nCol = columns.length;
    mat = new double[nRow][nCol];
    for (int col = 0; col < nCol; col++) {
      this.setCol(col, columns[col]);
    }
  }

  /**
   * One row per line, bound by brackets.  The entire matrix is also bound by brackets.
   * Example:
   * [[1 2 3]
   *  [4 5 6]]
   */
  public String toString() {
    String output = "[";
    for (int row = 0; row < nRow - 1; row++) {
      output += Arrays.toString(mat[row]) + "\n ";
    }
    return output + Arrays.toString(mat[nRow - 1]) + "]";
  }

  /**
   * If {@code this} is a square matrix, copy its content to a new {@link RealSquare} type
   *    object, in order to access methods specific to that class.
   * @return {@code RealSquare} the entries are identical to {@code this}.
   */
  public RealSquare toSquare() {
    if (nRow != nCol) {
      System.out.println("Not Square!");
      return null;
    }
    return new RealSquare(this.mat);
  }

  /**
   * Returns the value of a specific entry of {@code this}
   * @param row {@code int} the index of the row the entry belongs
   * @param col {@code int} the index of the column the entry belongs
   * @return {@code double} the value of the entry.
   */
  public double get(int row, int col) {
    return mat[row][col];
  }

  /**
   * Returns the value of a row vector of {@code this}
   * @param col {@code int} the index of the row to retrieve
   * @return {@code RealVector} the {@code col}th row vector of {@code this}.
   */
  public RealVector getRow(int row) {
    return new RealVector(mat[row]);
  }

  /**
   * Returns the value of a column vector of {@code this}
   * @param col {@code int} the index of the column to retrieve
   * @return {@code RealVector} the {@code col}th column vector of {@code this}.
   */
  public RealVector getCol(int col) {
    RealVector column = new RealVector(nRow);
    for(int row = 0; row < nRow; row++) {
      column.vec[row] = mat[row][col];
    }
    return column;
  }

  /**
   * Returns a block of {@code this}
   * @param initRow {@code int} index of the upper-most row of the block
   * @param initCol {@code int} index of the left-most column of the block
   * @param blockRow {@code int} number of rows of the block
   * @param blockCol {@code int} number of columns of the block
   * @return {@code RealMatrix} a block of {@code this}
   */
  public RealMatrix getBlock(int initRow, int initCol, int blockRow, int blockCol) {
    if ((initRow + blockRow > this.nRow)||(initCol + blockCol > this.nCol)) {
      System.out.println("Dimension limit exceeded.  No matrix retrieved.");
      return null;
    } else {
      RealMatrix block = new RealMatrix(blockRow, blockCol);
      for (int row = 0; row < blockRow; row++) {
        for (int col = 0; col < blockCol; col++) {
          block.mat[row][col] = this.mat[initRow+row][initCol+col];
        }
      }
      return block;
    }
  }

  /**
   * Returns a submatrix of {@code this} by deleting one row and one column
   * @param iRow {@code int} the index of the row to be deleted
   * @param iCol {@code int} the index of the column to be deleted
   * @return {@code RealMatrix} the (m-1) × (n-1) submatrix
   */
  public RealMatrix subMatrix(int iRow, int iCol) {
    int subRow;
    int subCol;
    RealMatrix sub = new RealMatrix(nRow - 1, nCol - 1);
    for (int row = 0; row < nRow; row++) {
      if (row < iRow) {
        subRow = row;
      } else if (row > iRow) {
        subRow = row - 1;
      } else {
        continue;
      }
      for (int col = 0; col < nCol; col++) {
        if (col < iCol) {
          subCol = col;
        } else if (col > iCol) {
          subCol = col - 1;
        } else {
          continue;
        }
        sub.mat[subRow][subCol] = this.mat[row][col];
      }
    }
    return sub;
  }

  /**
   * Retrieve the diagonal elements of {@code this}.
   * @return {@code double[]} an array whose length is the smaller one between the number of rows
   *         the number of columns of {@code this}.
   */
  public double[] getDiagonal() {
    int n = Math.min(nRow, nCol);
    double[] diagonal = new double[n];
    for (int i = 0; i < n; i++) {
      diagonal[i] = mat[i][i];
    }
    return diagonal;
  }

  /**
   * Returns a deep copy of {@code this}
   * @return {@code RealMatrix} a new, independent matrix object whose entries
   *         are identical to the corresponding entry of {@code this}
   */
  public RealMatrix copy() {
    return new RealMatrix(mat);
  }

  /**
   * Returns the dimension of {@code this}
   * @return {@code int[2]}, the {@code [0]} and {@code [1]} are the number of rows and of columns,
   *         respectively.
   */
  public int[] dim() {
    int[] shape = new int[2];
    shape[0] = nRow;
    shape[1] = nCol;
    return shape;
  }

  /**
   * Sets one specific entry of {@code this} to a new value.
   * @param row {@code int} the index of the row this entry belongs
   * @param col {@code int} the index of the column this entry belongs
   * @param value {@code double} the new value of the entry
   */
  public void set(int row, int col, double value) {
    mat[row][col] = value;
  }

  /**
   * Set the value of a row vector of {@code this}
   * @param col {@code int} the index of the row to be set
   * @param value {@code RealVector} the new value of the column vector of index {@code row}
   */
  public void setRow(int row, RealVector value) {
    for(int col = 0; col < nCol; col++) {
      mat[row][col] = value.vec[col];
    }
  }

  /**
   * Set the value of a column vector of {@code this}
   * @param col {@code int} the index of the column to be set
   * @param value {@code RealVector} the new value of the column vector of index {@code col}
   */
  public void setCol(int col, RealVector value) {
    for(int row = 0; row < nRow; row++) {
      mat[row][col] = value.vec[row];
    }
  }

  /**
   * Overlay a block of {@code this} with the content of another matrix.
   * @param initRow {@code int} the "upper boundary" of the block to be filled
   * @param initCol {@code int} the "left boundary" of the block to be filled
   * @param value {@code RealMatrix} dimension has to small enough so that its entirety
   *              can be filled into the block.
   */
  public void setBlock(int initRow, int initCol, RealMatrix value) {
    if ((initRow + value.nRow > this.nRow)||(initCol + value.nCol > this.nCol)) {
      System.out.println("Dimension limit exceeded.  No changes made.");
    } else {
      for (int row = 0; row < value.nRow; row++) {
        for (int col = 0; col < value.nCol; col++) {
          this.mat[initRow + row][initCol + col] = value.mat[row][col];
        }
      }
    }
  }

  /**
   * Equivalent to RealUtil.dot(this, B)
   * @param b {@code RealMatrix} {@code B.nRow == this.nCol}
   * @return {@code RealMatrix} {@code A.nRow} × {@code B.nCol}
   */
  public RealMatrix dot(RealMatrix B) {
    return RealUtil.dot(this, B);
  }

  /**
   * Equivalent to RealUtil.dot(this, b)
   * @param b {@code RealVector} {@code b.dim == this.nCol}
   * @return {@code RealVector} {@code dim == this.nRow}
   */
  public RealVector dot(RealVector b) {
    return RealUtil.dot(this, b);
  }
  
  /**
   * Transpose {@code this}
   * 
   * @return the transpose of {@code this}
   */
  public RealMatrix transpose() {
    RealMatrix result = new RealMatrix(nCol, nRow);
    for(int row = 0; row < nRow; row++) {
      result.setCol(row, getRow(row));
    }
    return result;
  }

  /**
   * Concatenate another matrix to {@code this}.  The dimension has to match.
   * @param B {@code RealMatrix} the other matrix to be concatenated
   * @param isHorizontal {@code boolean} if {@code true}, {@code B} will be concatenated to
   * the right of {@code this}.  Numbers of rows have to match.
   * <p> if {@code false}, {@code B} will be concatenated to the bottom of 
   * {@code this}.  Numbers of columns have to match
   * @return the concatenated matrix.
   */
  public RealMatrix concatenate(RealMatrix B, boolean isHorizontal) {
    RealMatrix joined;
    if(isHorizontal) {
      joined = new RealMatrix(nRow, nCol + B.nCol);
      joined.setBlock(0, 0, this);
      joined.setBlock(0, nCol, B);
    } else {
      joined = new RealMatrix(nRow + B.nRow, nCol);
      joined.setBlock(0, 0, this);
      joined.setBlock(nRow, 0, B);
    }
    return joined;
  }

  /**
   * One of the three elementary row operations: swap two rows
   * @param row1 {@code int} one of the two rows
   * @param row2 {@code int} the other of the two rows
   */
  public void rowSwap(int row1, int row2) {
    RealVector temp = getRow(row2);
    this.setRow(row2, getRow(row1));
    this.setRow(row1, temp);
  }

  /**
   * One of the three elementary row operations: multiply an entire row by a constant.
   * 
   * @param index {@code int} the index of the row to be operated on
   * @param k {@code double} the multiplier
   */
  public void rowMulti(int index, double k) {
    this.setRow(index, getRow(index).scalarMulti(k));
  }

  /**
   * One of the three elementary row operations: add a row multiplied by a constant
   * to another row
   * 
   * @param addTo {@code int} the row with this index will be added to by {@code k} times 
   * the row {@code added}.
   * @param added {@code int} the row with this index will not be changed
   * @param k {@code double} the multiplier
   */
  public void rowSum(int addTo, int added, double k) {
    for (int col = 0; col < nCol; col++)
      mat[addTo][col] += mat[added][col] * k;
  }

  /**
   * Returns {@code this} with each element multiplied by the same scalar k
   * @param k {@code double} the constant to be multiplied to
   * @return {@code RealMatrix} each entry is the entry of {@code this} multiplied by k
   */
  public RealMatrix scalarMulti(double k) {
    RealMatrix multiplied = this.copy();
    for (int row = 0; row < nRow; row++) {
      multiplied.rowMulti(row, k);
    }
    return multiplied;
  }

  /**
  * Convert the matrix to row-echelon form (REF) with Gaussian Elimination
  * algorithm.
  *
  * @return {@code RealMatrix} the REF of the current matrix.
  */
  public RealMatrix ref() {
    int tRow;   // In case a row swap is needed
    double k;
    int pRow = 0; // Location of the pivot
    int pCol = 0;
    RealMatrix echelon = copy(); // Make a copy so that the original is untouched

    while ((pRow < nRow - 1) && (pCol < nCol)) { // while is more manageable here than for
      tRow = pRow;
      while((Math.abs(echelon.get(tRow, pCol)) < RealUtil.TOL) && (tRow < nRow - 1)) {
        tRow++; // Find the first non-zero leading entry
      }
      if(tRow!=pRow)
        if(Math.abs(echelon.get(tRow, pCol)) < RealUtil.TOL) { // the rest of this column are 0
          pCol++; // Get to the next column
          continue;
        } else {
          echelon.rowSwap(pRow, tRow); // Ensure the pivot is non-zero
        }
        for(int row = pRow + 1;row < nRow; row++) { // Make the rest of the column 0
          k = -echelon.get(row, pCol) / echelon.get(pRow, pCol);
          echelon.rowSum(row, pRow, k);
        }
        pRow++;
        pCol++; // The possible leading entry of the next row
    }
    return echelon;
  }

  /**
   * Determines whether an entire row of {@code this} is made up of zeroes
   *    (within floating-point error)
   * @param row {@code int} the index of the row
   * @return {@code boolean} {@code true} iff the row vector of index {@code row} is a zero vector
   */
  public boolean isZeroRow(int row) {
    return this.getRow(row).isZero();
  }

  /**
   * Determines whether a specific column vector of {@code this} is a zero vector.
   * @param col {@code int} the index of the column
   * @return {@code boolean} {@code true} iff the column vector of index {@code col}
   *         is a zero vector. 
   */
  public boolean isZeroCol(int col) {
    return this.getCol(col).isZero();
  }

  /**
   * Determines whether {@code this} has full row rank
   * @return {@code boolean} true iff {@code this} has full row rank
   */
  public boolean isFullRowRank() {
      if (rref == null) {
        rref = new RealRref(this);
      }
      return rref.isFullRowRank();
  }

  /**
   * Determines whether {@code this} has full column rank
   * @return {@code boolean} true iff {@code this} has full column rank
   */
  public boolean isFullColRank() {
    if (rref == null) {
      rref = new RealRref(this);
    }
    return rref.isFullColRank();
  }

  /**
   * Determines whether the column vectors of {@code this} are linearly independent.
   * @return {@code boolean} {@code true} iff the column vectors of {@code this} are LI
   */
  public boolean isLinearlyIndependent() {
    return isFullColRank(); // A matrix has LI column vectors iff it has full column rank.
  }

  /**
   * Determines whether the column vectors of {@code this} are orthogonal to each other
   * @return {@code boolean} true iff any pair of column vectors of {@code this} 
   *         are orthogonal to each other.
   */
  public boolean isOrthogonal() {
    if (!isFullColRank()) {
     return false; 
    }
    for (int i = 0; i < nCol - 1; i++) {
      for (int j = i + 1; j < nCol; j++) { // Total number of pairs: n(n-1)/2
        if (!RealUtil.isOrthogonal(this.getCol(i), this.getCol(j))) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Determines whether {@code this} is column orthogonal (column vectors are orthogonal to each other
   *    and each column vector is normalized).
   * <p> To check if a matrix A is row orthogonal, use: {@code A.transpose().isOrthonormal()}
   * @return {@code boolean} {@code true} iff {@code this} is column orthogonal.
   */
  public boolean isOrthonormal() {
    // For an column orthonormal matrix A, A^{T} · A is an n × n identity matrix
    return RealUtil.matrixCompare(this.transpose().dot(this), RealUtil.eye(nCol));
  }

  /**
   * Returns the rank of {@code this}.
   * @return {@code int} rank({@code this})
   */
  public int rank() {
    if (rref == null) {
      rref = new RealRref(this);
    } // This is to ensure that the Gauss-Jordan procedure is only done once.
    return rref.rank;
  }

  /**
   * Returns the nullity of {@code this}
   * @return {@code int} the nullity of {@code this}
   */
  public int nullity() {
    if (rref == null) {
      rref = new RealRref(this);
    }
    return rref.nullity();
  }

  /**
   * Returns a set of vectors that forms a basis of the row space of {@code this}.
   * @return {@code RealVector[]} elements form a basis of the row space of {@code this}
   */
  public RealVector[] rowBasisReduced() {
    if (rref == null) {
      rref = new RealRref(this);
    }
    RealVector[] basis = new RealVector[rref.rank];
    for (int row = 0; row < rref.rank; row++) {
      basis[row] = rref.getRow(row);
    }
    return basis;
  }

  /**
   * Returns a subset of the column vectors of {@code this} that forms a basis
   *  of the column space of {@code this}.
   * @return {@code RealVector[]} elements are column vectors of {@code this}
   *         and form a basis of the column space of {@code this}
   */
  public RealVector[] colBasis() {
    // Extract the columns in the original matrix corresponding to the pivot columns of its RREF.
    if (rref == null) {
      rref = new RealRref(this);
    }
    RealVector[] basis = new RealVector[rref.rank];
    int i = 0;
    for (int col = 0; col < nCol; col++) {
      if (rref.pivots[col]) {
        basis[i] = getCol(col);
        i++;
      }
    }
    return basis;
  }

  /**
   * Returns a subset of the row vectors of {@code this} that forms a basis
   *  of the row space of {@code this}.
   * @return {@code RealVector[]} elements are row vectors of {@code this}
   *         and form a basis of the row space of {@code this}
   * @see {@link #rowBasisReduced()}, which returns a subset of the row vectors of the RREF of
   *      {@code this}, instead of the those in {@code this} itself.
   */
  public RealVector[] rowBasis() {
    return this.transpose().colBasis(); // Yup it's this simple ^_^
  }

  /**
   * Returns a basis of the nullspace of {@code this}.  The basis is probably not orthonormal.
   * @return {@code RealVector[]} elements form a basis of the nullspace of {@code this}.
   */
  public RealVector[] nullBasis() {
    if (rref == null) {
      rref = new RealRref(this);
    }
    int rank = rref.rank;
    int nullity = rref.nullity();
    RealVector[] basis = new RealVector[nullity];

    // For the m × n matrix A of rank r, a matrix N whose columns form a basis of null(A)
    //    has dimensions n × (n-r).  A convenient choice of N is constructed as follows:
    // Firstly, construct an r × (n-r) matrix U, whose columns are the first r entries of the 
    //    non-pivotal columns of the RREF of A.  Also construct an (n-r) × (n-r) identity matrix I.
    // For any 1 ≤ i ≤ n, consider whether the ith column of A is a pivot column in its RREF:
    //    1. If yes, the ith row of N is a row of U
    //    2. If no, the ith row of N is a row of I.
    // If the pivot columns of the RREF of A are exactly the first r columns,
    //    N will look like U sitting on top of an identity matrix.
    RealMatrix upper = new RealMatrix(rank, nullity);
    int i = 0;
    for (int col = 0; col < nCol; col++) {
      if (!rref.pivots[col]) {
        for (int row = 0; row < rank; row++) {
          upper.mat[row][i] = -rref.mat[row][col];
        }
        i++;
      }
    }
    RealSquare identity = RealUtil.eye(nullity);
    RealMatrix basisMatrix = new RealMatrix(nCol, nullity);
    int indexUpper = 0;
    int indexIdentity = 0;
    for (i = 0; i < nCol; i++) {
      if (rref.pivots[i]) {
        basisMatrix.setRow(i, upper.getRow(indexUpper));
        indexUpper++;
      } else {
        basisMatrix.setRow(i, identity.getRow(indexIdentity));
        indexIdentity++;
      }
    }
    for (i = 0; i < nullity; i++) {
      basis[i] = basisMatrix.getCol(i);
    }
    return basis;
  }

  /**
   * Returns the projection matrix corresponding to {@code this} P, 
   *    such that for any vector v whose dimension equals to {@code this.nCol}
   *    Pv is the projection of v on the space spanned by the column vectors of {@code this}
   * @return {@code RealSquare} n × n matrix P, where n = {@code this.nCol}
   */
  public RealSquare toProjectionMatrix() {
    RealSquare AtA = RealUtil.dot(this.transpose(), this).toSquare(); // A^T A
    return this.dot(AtA.inverse().dot(this.transpose())).toSquare(); // A (A^T A)^{-1} A^T
  }
}