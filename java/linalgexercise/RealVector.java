package linalgexercise;

//import java.io.*;
import java.util.*;
import java.lang.Math;

public class RealVector {
  protected double[] vec;
  protected int dim;

  /**
   * Class constructor specifying the value of each component in an array
   * @param init {@code double[]} the values of the components
   */
  public RealVector(double[] init) {
    dim = init.length;
    vec = new double[dim];
    for (int i = 0; i < dim; i++) {
      vec[i] = init[i];
    }
  }

  /**
   * Class constructor specifying the dimension of the vector,
   *    and initialize each component as 0.0.
   * @param size {@int size} positive integer, the dimension of {@code this}
   */
  public RealVector(int size) {
    this.dim = size;
    vec = new double[size];
  }

  /**
   * Returns a string representing the vector in the format [1, 2, 3]
   */
  public String toString() {
    return Arrays.toString(vec);
  }

  /**
   * Returns {@code this} as a row or column matrix
   * @param isColumn {@code boolean} if {@code true}, returns a column matrix;
   *                 if {@code false}, returns a row matrix.
   * @return {@code RealMatrix} a row or column matrix
   */
  public RealMatrix toMatrix(boolean isColumn) {
    if (isColumn) {
      RealMatrix columnMatrix = new RealMatrix(dim, 1); // column matrix has 1 column
      columnMatrix.setCol(0, this);
      return columnMatrix;
    } else {
      RealMatrix rowMatrix = new RealMatrix(1, dim); // row matrix has 1 row
      rowMatrix.setRow(0, this);
      return rowMatrix;
    }
  }

  /**
   * Returns the value of a specific component of {@code this}
   * @param i {@code int} the index of that component (starting from zero)
   * @return {@code double} the value of the {@code i}th component of {@code this}
   */
  public double get(int i) {
    return vec[i];
  }

  /**
   * Returns the dimension of {@code this}
   * @return {@code int} the dimension of {@code this}
   */
  public int dim() {
    return dim;
  }

  /**
   * Modify the value of a specific component of {@code this}.
   * @param i {@code int} index of the component
   * @param value {@code double} the new value
   */
  public void set(int i, double value) {
    vec[i] = value;
  }

  /**
   * Make a deep copy of {@code this}
   * @return {@code RealVector} components are identical to {@code this}.
   */
  public RealVector copy() {
    // construct a new RealVector that initially takes the same value,
    // but is independent of the current one
    return new RealVector(vec);
  }

  /**
   * Equivalent to {@code RealUtil.dot(this, other)}
   * @param other {@code RealVector} another vector of the same dimension
   * @return {@code double} this · other
   */
  public double dot(RealVector other) {
    return RealUtil.dot(this, other);
  }

  /**
   * Equivalent to {@code RealUtil.outer(this, other)}
   * @param other {@code RealVector} another vector that does NOT necessarily have
   *              the same dimension as {@code this}
   * @return {@code RealMatrix} this (as a column vector) · other (as a row vector)
   * @see {@link RealUtil#outer(RealVector, RealVector)}
   */
  public RealMatrix outer(RealVector other) {
    return RealUtil.outer(this, other);
  }

  /**
   * Equivalent to {@code RealUtil.outer(this, other)}
   * @param other {@code RealMatrix} {@code other.nRow == this.dim}
   * @return {@code RealVector} this (as a row vector) · other
   * @see {@link RealUtil#outer(RealVector, RealMatrix)}
   */
  public RealVector outer(RealMatrix other) {
    return RealUtil.outer(this, other);
  }

  /**
   * Equivalent to {@code RealUtil.cross(this, other)}
   * @param other {@code RealVector} another vector of the same dimension as {@code this}
   * @return {@code RealVector} this × other
   * @see {@link RealUtil#cross(RealVector, RealVector)}
   */
  public RealVector cross(RealVector other) {
    return RealUtil.cross(this, other);
  }

  /**
   * Returns the modulus (magnitude, length) of {@code this}
   * @return {@code double} the modulus of {@code this}
   */
  public double modulus() {
    return Math.sqrt(RealUtil.dot(this, this));
  }

  /**
   * Determines whether {@code this} is a zero vector within floating-point error
   * @return {@code boolean} {@code}true iff {@code this} is a zero vector
   */
  public boolean isZero() {
    // No need to check each element, more code and not much more efficiency
    return this.modulus() < RealUtil.TOL;
  }

  /**
   * Returns a scalar multiplication of {@code this} with scalar k
   * @param k {@code double} the scalar to be multiplied with
   * @return {@code RealVector} k * {@code this}
   */
  public RealVector scalarMulti(double k) {
    RealVector product = this.copy();
    for(int i = 0; i < dim; i++) {
      product.vec[i] *= k;
    }
    return product;
  }

  /**
   * Normalize {@code this}.
   * @return {@code RealVector} unit vector in the direction of {@code this}
   */
  public RealVector normalize() {
    double modReciprocal = 1.0/this.modulus();
    return this.scalarMulti(modReciprocal); 
  }

  /**
   * Returns the scalar projection of {@code this} onto another vector {@code basis}
   * @param basis {@code RealVector} the vector to be projected onto
   * @return {@code double} the scalar projection of {@code this} onto {@code basis}
   */
  public double scalarProjection(RealVector basis) {
    return this.dot(basis) / basis.modulus();
  }

  /**
   * Returns the vector projection of {@code this} onto another vector {@code basis}
   * @param basis {@code RealVector} the vector to be projected onto
   * @return {@code RealVector} the vector projection of {@code this} onto {@code basis}
   */
  public RealVector vectorProjection(RealVector basis) {
    // ratio of the length of the projection to the length of the basis
    double k = RealUtil.dot(this, basis)/RealUtil.dot(basis, basis);
    return basis.copy().scalarMulti(k);
  }

  /**
   * Returns the vector projection of {@code this} onto the column space of {@code spanningSet}
   *    by turning the latter into an orthogonolized basis first, and then express the projection
   *    as a linear combination of the basis. 
   * @param spanningSet {@code RealMatrix} a matrix with {@code this.dim} rows.
   *                    Its column vectors neither have to be linearly independent
   *                    nor have to be orthogonal.
   * @return {@code RealVector} projection of {@code this} onto the column space of
   *         {@code spanning set}
   */
  public RealVector vectorProjectionRaw(RealMatrix spanningSet) {
    // spanningSet may be unorthogonolized or even LD, but it doesn't have to be normalized.
    if (!spanningSet.isOrthogonal()) {
      spanningSet = new RealQr(spanningSet).getQ();
    }
    double[] coefficients = new double[spanningSet.nCol];
    for (int i = 0; i < spanningSet.nCol; i++) {
      coefficients[i] = scalarProjection(spanningSet.getCol(i));
    }
    return RealUtil.linComb(spanningSet, coefficients);
  }

  /**
   * Returns the vector projection of v ({@code this}) onto the column space of
   *    {@code spanningSet} by solving the equation A^T A x = A^T v, where the matrix A
   *    consists of the independent columns of {@code spanningSet}.
   * @param spanningSet {@code RealMatrix} a matrix with {@code this.dim} rows.
   *                    Its column vectors does not have to be linearly independent.
   * @return {@code RealVector} projection of {@code this} onto the column space of
   *         {@code spanning set}
   */
  public RealVector vectorProjection(RealMatrix spanningSet) {
    // In case the column vectors of spanningSet are linearly dependent:
    RealMatrix basis = new RealMatrix(spanningSet.colBasis()); 
    RealSquare A = RealUtil.dot(basis.transpose(), basis).toSquare();
    RealVector b = RealUtil.dot(basis.transpose(), this);
    return RealUtil.dot(basis, new RealPldu(A).solve(b)); // Don't forget the multiplication
  }
}