package linalgexercise;

/**
 * Stores a matrix in Row Reduced Echelon Form (RREF), and extract (or allowing the extraction of)
 *    crucial properties of the original matrix, such as rank, nullity, row space, column space,
 *    null space, singularity etc., and serve as the foundation of several algorithms for
 *    solving linear systems and inverting matrices.
 * 
 * @author Zheng Ma
 * @version %I%, %G%
 * @since 1.0
 */
public class RealRref extends RealMatrix{
  protected int rank;
  // pivots[i] is true iff there exists a row in RREF for which the pivot is in the ith column.
  protected boolean[] pivots;

  /**
  * Store the reduced row-echelon form (RREF) of {@code A} in {@code this.mat}, and fill in
  * the {@code pivots[]} and {@code rank}.
  *
  * @param A {@code RealMatrix} m × n matrix
  */
  public RealRref(RealMatrix A) {
    super(A.ref().mat); // Start with the REF
    rref(); // Gauss-Jordan elimination
    pivots = new boolean[nCol]; 
    findPivots(); // Record the pivot columns
    findRank(); // Calculate and store the rank
  }

  /**
  * Convert {@code this.mat} to reduced row-echelon form (RREF) with Gauss-Jordan
  * Elimination algorithm.
  */
  private void rref() {
    int pCol;
    double k;
    for(int pRow = Math.min(nRow, nCol) - 1; pRow >= 0; pRow--) {
      pCol = pRow;
      // Find the leading entry
      while((Math.abs(mat[pRow][pCol]) < RealUtil.TOL) && (pCol < nCol - 1)) {
        pCol++;
      }
      if(Math.abs(mat[pRow][pCol]) < RealUtil.TOL) { // Move up if the entire row is 0
        continue;
      }
      k = 1.0 / mat[pRow][pCol];
      this.rowMulti(pRow, k); // The leanding entry will automatically become 1.0
      for(int row = pRow - 1; row >= 0; row--){ // Back substitution
        k = -mat[row][pCol];
        this.rowSum(row, pRow, k);
      }
    }
  }

  /**
   * Find the pivot of each row, and store them in the boolean array {@code pivots}
   */
  private void findPivots() {
    int col = -1; // Such that at the beginning of the first round of for loop col++ becomes 0
    for (int row = 0; row < Math.min(nRow, nCol); row++) {
      col++; // The pivot next row is always to the right of the pivot of the current row
      // Check (col < nCol - 1) first to avoid ArrayIndexOutOfBounds error
      while ((col < nCol - 1) && (Math.abs(mat[row][col]) < RealUtil.TOL)) {
        col++; // Search for the first non-zero entry of the row
      }
      // If the last entry of the row is still zero, this row and ANY row beneath it are zero rows.
      if (col > nCol - 1) {
        break;
      }
      pivots[col] = true;
    }
  }

  /**
   * Calculate the rank of {@code this.mat} and store it in {@code this.rank}
   */
  private void findRank() {
    this.rank = 0;
    for (int i = 0; i < nCol; i++) {
      if (pivots[i]) {
        this.rank++;
      }
    }
  }

  /**
   * Return the nullity of {@code this.mat}
   * @return {@code int} nullity of {@code this.mat}
   */
  public int nullity() {
    return nCol - rank;
  }

  /**
   * Determine whether {@code this.mat} has full row rank
   * @return {@code boolean} {@code true} iff {@code this.mat} has full row rank
   */
  public boolean isFullRowRank() {
    return rank == nRow;
  }

  /**
   * Determine whether {@code this.mat} has full column rank
   * @return {@code boolean} {@code true} iff {@code this.mat} has full column rank
   */
  public boolean isFullColRank() {
    return rank == nCol;
  }
}
