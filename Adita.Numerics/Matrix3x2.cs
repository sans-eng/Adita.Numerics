//MIT License

//Copyright (c) 2022 Adita

//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

using System.Globalization;

namespace Adita.Numerics
{
    /// <summary>
    /// Represents a matrix 3x2.
    /// </summary>
    public struct Matrix3x2
    {
        #region Constructors
        /// <summary>
        /// Creates a new <see cref="Matrix3x2"/> that contains specified components.
        /// </summary>
        public Matrix3x2(double m11, double m12,
                         double m21, double m22,
                         double m31, double m32)
        {
            M11 = m11;
            M12 = m12;
            M21 = m21;
            M22 = m22;
            M31 = m31;
            M32 = m32;
        }
        #endregion Constructors


        #region Public properties
        /// <summary>
        /// Gets the value of the first row and first column of current <see cref="Matrix3x2"/>.
        /// </summary>
        public double M11 { get; }
        /// <summary>
        /// Gets the value of the first row and second column of current <see cref="Matrix3x2"/>.
        /// </summary>
        public double M12 { get; }
        /// <summary>
        /// Gets the value of the second row and first column of current <see cref="Matrix3x2"/>.
        /// </summary>
        public double M21 { get; }
        /// <summary>
        /// Gets the value of the second row and second column of current <see cref="Matrix3x2"/>.
        /// </summary>
        public double M22 { get; }
        /// <summary>
        /// Gets the value of the third row and first column of current <see cref="Matrix3x2"/>.
        /// </summary>
        public double M31 { get; private set; }
        /// <summary>
        /// Gets the value of the third row and second column of current <see cref="Matrix3x2"/>.
        /// </summary>
        public double M32 { get; private set; }

        /// <summary>
        /// Gets the multiplicative identity matrix.
        /// </summary>
        public static Matrix3x2 Identity { get; } = new(
            1d, 0d,
            0d, 1d,
            0d, 0d
        );

        /// <summary>
        /// Gets whether the matrix is the identity matrix.
        /// </summary>
        public bool IsIdentity
        {
            get
            {
                return M11 == 1f && M22 == 1f && // Check diagonal element first for early out.
                                    M12 == 0f &&
                       M21 == 0f &&
                       M31 == 0f && M32 == 0f;
            }
        }

        /// <summary>
        /// Gets or sets the translation component of this matrix.
        /// </summary>
        public Vector2 Translation
        {
            get
            {
                return new Vector2(M31, M32);
            }

            set
            {
                M31 = value.X;
                M32 = value.Y;
            }
        }
        #endregion Public properties

        #region Public methods

        /// <summary>
        /// Creates a translation matrix from the given vector.
        /// </summary>
        /// <param name="position">The translation position.</param>
        /// <returns>A translation matrix.</returns>
        public static Matrix3x2 CreateTranslation(Vector2 position)
        {
            return new Matrix3x2(1.0d, 0.0d, 0.0d, 1.0d, position.X, position.Y);
        }

        /// <summary>
        /// Creates a translation matrix from the given X and Y components.
        /// </summary>
        /// <param name="xPosition">The X position.</param>
        /// <param name="yPosition">The Y position.</param>
        /// <returns>A translation matrix.</returns>
        public static Matrix3x2 CreateTranslation(double xPosition, double yPosition)
        {
            return new Matrix3x2(1.0d, 0.0d, 0.0d, 1.0d, xPosition, yPosition);
        }

        /// <summary>
        /// Creates a scale matrix from the given X and Y components.
        /// </summary>
        /// <param name="xScale">Value to scale by on the X-axis.</param>
        /// <param name="yScale">Value to scale by on the Y-axis.</param>
        /// <returns>A scaling matrix.</returns>
        public static Matrix3x2 CreateScale(double xScale, double yScale)
        {
            return new Matrix3x2(xScale, 0.0d, 0.0d, yScale, 0.0d, 0.0d);
        }

        /// <summary>
        /// Creates a scale matrix that is offset by a given center point.
        /// </summary>
        /// <param name="xScale">Value to scale by on the X-axis.</param>
        /// <param name="yScale">Value to scale by on the Y-axis.</param>
        /// <param name="centerPoint">The center point.</param>
        /// <returns>A scaling matrix.</returns>
        public static Matrix3x2 CreateScale(double xScale, double yScale, Vector2 centerPoint)
        {
            double tx = centerPoint.X * (1 - xScale);
            double ty = centerPoint.Y * (1 - yScale);

            return new Matrix3x2(xScale, 0.0d, 0.0d, yScale, tx, ty);
        }

        /// <summary>
        /// Creates a scale matrix from the given vector scale.
        /// </summary>
        /// <param name="scales">The scale to use.</param>
        /// <returns>A scaling matrix.</returns>
        public static Matrix3x2 CreateScale(Vector2 scales)
        {
            return new Matrix3x2(scales.X, 0.0d, 0.0d, scales.Y, 0.0d, 0.0d);
        }

        /// <summary>
        /// Creates a scale matrix from the given vector scale with an offset from the given center point.
        /// </summary>
        /// <param name="scales">The scale to use.</param>
        /// <param name="centerPoint">The center offset.</param>
        /// <returns>A scaling matrix.</returns>
        public static Matrix3x2 CreateScale(Vector2 scales, Vector2 centerPoint)
        {
            double tx = centerPoint.X * (1 - scales.X);
            double ty = centerPoint.Y * (1 - scales.Y);

            return new Matrix3x2(scales.X, 0.0d, 0.0d, scales.Y, tx, ty);
        }

        /// <summary>
        /// Creates a scale matrix that scales uniformly with the given scale.
        /// </summary>
        /// <param name="scale">The uniform scale to use.</param>
        /// <returns>A scaling matrix.</returns>
        public static Matrix3x2 CreateScale(double scale)
        {
            return new Matrix3x2(scale, 0.0d, 0.0d, scale, 0.0d, 0.0d);
        }

        /// <summary>
        /// Creates a scale matrix that scales uniformly with the given scale with an offset from the given center.
        /// </summary>
        /// <param name="scale">The uniform scale to use.</param>
        /// <param name="centerPoint">The center offset.</param>
        /// <returns>A scaling matrix.</returns>
        public static Matrix3x2 CreateScale(double scale, Vector2 centerPoint)
        {
            double tx = centerPoint.X * (1 - scale);
            double ty = centerPoint.Y * (1 - scale);

            return new Matrix3x2(scale, 0.0d, 0.0d, scale, tx, ty);
        }

        /// <summary>
        /// Creates a skew matrix from the given angles in radians.
        /// </summary>
        /// <param name="radiansX">The X angle, in radians.</param>
        /// <param name="radiansY">The Y angle, in radians.</param>
        /// <returns>A skew matrix.</returns>
        public static Matrix3x2 CreateSkew(double radiansX, double radiansY)
        {
            double xTan = Math.Tan(radiansX);
            double yTan = Math.Tan(radiansY);

            return new Matrix3x2(1.0d, yTan, xTan, 1.0d, 0.0d, 0.0d);
        }

        /// <summary>
        /// Creates a skew matrix from the given angles in radians and a center point.
        /// </summary>
        /// <param name="radiansX">The X angle, in radians.</param>
        /// <param name="radiansY">The Y angle, in radians.</param>
        /// <param name="centerPoint">The center point.</param>
        /// <returns>A skew matrix.</returns>
        public static Matrix3x2 CreateSkew(double radiansX, double radiansY, Vector2 centerPoint)
        {
            double xTan = (float)Math.Tan(radiansX);
            double yTan = (float)Math.Tan(radiansY);

            double tx = -centerPoint.Y * xTan;
            double ty = -centerPoint.X * yTan;

            return new Matrix3x2(1.0d, yTan, xTan, 1.0d, tx, ty);
        }

        /// <summary>
        /// Creates a rotation matrix using the given rotation in radians.
        /// </summary>
        /// <param name="radians">The amount of rotation, in radians.</param>
        /// <returns>A rotation matrix.</returns>
        public static Matrix3x2 CreateRotation(double radians)
        {
            radians = Math.IEEERemainder(radians, Math.PI * 2);

            double c, s;

            const double epsilon = 0.001d * Math.PI / 180d;     // 0.1% of a degree

            if (radians > -epsilon && radians < epsilon)
            {
                // Exact case for zero rotation.
                c = 1;
                s = 0;
            }
            else if (radians > (Math.PI / 2) - epsilon && radians < (Math.PI / 2) + epsilon)
            {
                // Exact case for 90 degree rotation.
                c = 0;
                s = 1;
            }
            else if (radians < -Math.PI + epsilon || radians > Math.PI - epsilon)
            {
                // Exact case for 180 degree rotation.
                c = -1;
                s = 0;
            }
            else if (radians > (-Math.PI / 2) - epsilon && radians < (-Math.PI / 2) + epsilon)
            {
                // Exact case for 270 degree rotation.
                c = 0;
                s = -1;
            }
            else
            {
                // Arbitrary rotation.
                c = Math.Cos(radians);
                s = Math.Sin(radians);
            }

            // [  c  s ]
            // [ -s  c ]
            // [  0  0 ]

            return new Matrix3x2(c, s, -s, c, 0.0d, 0.0d);
        }

        /// <summary>
        /// Creates a rotation matrix using the given rotation in radians and a center point.
        /// </summary>
        /// <param name="radians">The amount of rotation, in radians.</param>
        /// <param name="centerPoint">The center point.</param>
        /// <returns>A rotation matrix.</returns>
        public static Matrix3x2 CreateRotation(double radians, Vector2 centerPoint)
        {
            radians = Math.IEEERemainder(radians, Math.PI * 2);

            double c, s;

            const double epsilon = 0.001f * Math.PI / 180f;     // 0.1% of a degree

            if (radians > -epsilon && radians < epsilon)
            {
                // Exact case for zero rotation.
                c = 1;
                s = 0;
            }
            else if (radians > (Math.PI / 2) - epsilon && radians < (Math.PI / 2) + epsilon)
            {
                // Exact case for 90 degree rotation.
                c = 0;
                s = 1;
            }
            else if (radians < -Math.PI + epsilon || radians > Math.PI - epsilon)
            {
                // Exact case for 180 degree rotation.
                c = -1;
                s = 0;
            }
            else if (radians > (-Math.PI / 2) - epsilon && radians < (-Math.PI / 2) + epsilon)
            {
                // Exact case for 270 degree rotation.
                c = 0;
                s = -1;
            }
            else
            {
                // Arbitrary rotation.
                c = Math.Cos(radians);
                s = Math.Sin(radians);
            }

            double x = (centerPoint.X * (1 - c)) + (centerPoint.Y * s);
            double y = (centerPoint.Y * (1 - c)) - (centerPoint.X * s);

            // [  c  s ]
            // [ -s  c ]
            // [  x  y ]

            return new Matrix3x2(c, s, -s, c, x, y);
        }

        /// <summary>
        /// Calculates the determinant for this matrix. 
        /// The determinant is calculated by expanding the matrix with a third column whose values are (0,0,1).
        /// </summary>
        /// <returns>The determinant.</returns>
        public double GetDeterminant()
        {
            // There isn't actually any such thing as a determinant for a non-square matrix,
            // but this 3x2 type is really just an optimization of a 3x3 where we happen to
            // know the rightmost column is always (0, 0, 1). So we expand to 3x3 format:
            //
            //  [ M11, M12, 0 ]
            //  [ M21, M22, 0 ]
            //  [ M31, M32, 1 ]
            //
            // Sum the diagonal products:
            //  (M11 * M22 * 1) + (M12 * 0 * M31) + (0 * M21 * M32)
            //
            // Subtract the opposite diagonal products:
            //  (M31 * M22 * 0) + (M32 * 0 * M11) + (1 * M21 * M12)
            //
            // Collapse out the constants and oh look, this is just a 2x2 determinant!

            return (M11 * M22) - (M21 * M12);
        }

        /// <summary>
        /// Attempts to invert the given matrix. If the operation succeeds, the inverted matrix is stored in the result parameter.
        /// </summary>
        /// <param name="matrix">The source matrix.</param>
        /// <param name="result">The output matrix.</param>
        /// <returns>True if the operation succeeded, False otherwise.</returns>
        public static bool Invert(Matrix3x2 matrix, out Matrix3x2 result)
        {
            double det = (matrix.M11 * matrix.M22) - (matrix.M21 * matrix.M12);

            if (Math.Abs(det) < double.Epsilon)
            {
                result = new Matrix3x2(double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN);
                return false;
            }

            double invDet = 1.0d / det;

            double m11 = matrix.M22 * invDet;
            double m12 = -matrix.M12 * invDet;
            double m21 = -matrix.M21 * invDet;
            double m22 = matrix.M11 * invDet;
            double m31 = ((matrix.M21 * matrix.M32) - (matrix.M31 * matrix.M22)) * invDet;
            double m32 = ((matrix.M31 * matrix.M12) - (matrix.M11 * matrix.M32)) * invDet;

            result = new Matrix3x2(m11, m12, m21, m22, m31, m32);

            return true;
        }

        /// <summary>
        /// Linearly interpolates from matrix1 to matrix2, based on the third parameter.
        /// </summary>
        /// <param name="matrix1">The first source matrix.</param>
        /// <param name="matrix2">The second source matrix.</param>
        /// <param name="amount">The relative weighting of matrix2.</param>
        /// <returns>The interpolated matrix.</returns>
        public static Matrix3x2 Lerp(Matrix3x2 matrix1, Matrix3x2 matrix2, double amount)
        {
            // First row
            double m11 = matrix1.M11 + ((matrix2.M11 - matrix1.M11) * amount);
            double m12 = matrix1.M12 + ((matrix2.M12 - matrix1.M12) * amount);

            // Second row
            double m21 = matrix1.M21 + ((matrix2.M21 - matrix1.M21) * amount);
            double m22 = matrix1.M22 + ((matrix2.M22 - matrix1.M22) * amount);

            // Third row
            double m31 = matrix1.M31 + ((matrix2.M31 - matrix1.M31) * amount);
            double m32 = matrix1.M32 + ((matrix2.M32 - matrix1.M32) * amount);

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Negates the given matrix by multiplying all values by -1.
        /// </summary>
        /// <param name="value">The source matrix.</param>
        /// <returns>The negated matrix.</returns>
        public static Matrix3x2 Negate(Matrix3x2 value)
        {
            double m11 = -value.M11;
            double m12 = -value.M12;
            double m21 = -value.M21;
            double m22 = -value.M22;
            double m31 = -value.M31;
            double m32 = -value.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Adds each matrix element in value1 with its corresponding element in value2.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>The matrix containing the summed values.</returns>
        public static Matrix3x2 Add(Matrix3x2 value1, Matrix3x2 value2)
        {
            double m11 = value1.M11 + value2.M11;
            double m12 = value1.M12 + value2.M12;
            double m21 = value1.M21 + value2.M21;
            double m22 = value1.M22 + value2.M22;
            double m31 = value1.M31 + value2.M31;
            double m32 = value1.M32 + value2.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Subtracts each matrix element in value2 from its corresponding element in value1.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>The matrix containing the resulting values.</returns>
        public static Matrix3x2 Subtract(Matrix3x2 value1, Matrix3x2 value2)
        {
            double m11 = value1.M11 - value2.M11;
            double m12 = value1.M12 - value2.M12;
            double m21 = value1.M21 - value2.M21;
            double m22 = value1.M22 - value2.M22;
            double m31 = value1.M31 - value2.M31;
            double m32 = value1.M32 - value2.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Multiplies two matrices together and returns the resulting matrix.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>The product matrix.</returns>
        public static Matrix3x2 Multiply(Matrix3x2 value1, Matrix3x2 value2)
        {

            // First row
            double m11 = (value1.M11 * value2.M11) + (value1.M12 * value2.M21);
            double m12 = (value1.M11 * value2.M12) + (value1.M12 * value2.M22);

            // Second row
            double m21 = (value1.M21 * value2.M11) + (value1.M22 * value2.M21);
            double m22 = (value1.M21 * value2.M12) + (value1.M22 * value2.M22);

            // Third row
            double m31 = (value1.M31 * value2.M11) + (value1.M32 * value2.M21) + value2.M31;
            double m32 = (value1.M31 * value2.M12) + (value1.M32 * value2.M22) + value2.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Scales all elements in a matrix by the given scalar factor.
        /// </summary>
        /// <param name="value1">The source matrix.</param>
        /// <param name="value2">The scaling value to use.</param>
        /// <returns>The resulting matrix.</returns>
        public static Matrix3x2 Multiply(Matrix3x2 value1, double value2)
        {
            double m11 = value1.M11 * value2;
            double m12 = value1.M12 * value2;
            double m21 = value1.M21 * value2;
            double m22 = value1.M22 * value2;
            double m31 = value1.M31 * value2;
            double m32 = value1.M32 * value2;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }
        /// <summary>
        /// Returns a boolean indicating whether the matrix is equal to the other given matrix.
        /// </summary>
        /// <param name="other">The other matrix to test equality against.</param>
        /// <returns>True if this matrix is equal to other; False otherwise.</returns>
        public bool Equals(Matrix3x2 other)
        {
            return M11 == other.M11 &&
                    M22 == other.M22 && // Check diagonal element first for early out.
                    M12 == other.M12 &&
                    M21 == other.M21 &&
                    M31 == other.M31 &&
                    M32 == other.M32;
        }

        /// <summary>
        /// Returns a boolean indicating whether the given Object is equal to this matrix instance.
        /// </summary>
        /// <param name="obj">The Object to compare against.</param>
        /// <returns>True if the Object is equal to this matrix; False otherwise.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is Matrix3x2 matrix3x2)
            {
                return Equals(matrix3x2);
            }

            return false;
        }

        /// <summary>
        /// Returns a String representing this matrix instance.
        /// </summary>
        /// <returns>The string representation.</returns>
        public override string ToString()
        {
            CultureInfo ci = CultureInfo.CurrentCulture;
            return String.Format(ci, "{{ {{M11:{0} M12:{1}}} {{M21:{2} M22:{3}}} {{M31:{4} M32:{5}}} }}",
                                 M11.ToString(ci), M12.ToString(ci),
                                 M21.ToString(ci), M22.ToString(ci),
                                 M31.ToString(ci), M32.ToString(ci));
        }

        /// <summary>
        /// Returns the hash code for this instance.
        /// </summary>
        /// <returns>The hash code.</returns>
        public override int GetHashCode()
        {
            return M11.GetHashCode() + M12.GetHashCode() +
                   M21.GetHashCode() + M22.GetHashCode() +
                   M31.GetHashCode() + M32.GetHashCode();
        }
        #endregion Public methods

        #region Operators

        /// <summary>
        /// Negates the given matrix by multiplying all values by -1.
        /// </summary>
        /// <param name="value">The source matrix.</param>
        /// <returns>The negated matrix.</returns>
        public static Matrix3x2 operator -(Matrix3x2 value)
        {
            double m11 = -value.M11;
            double m12 = -value.M12;
            double m21 = -value.M21;
            double m22 = -value.M22;
            double m31 = -value.M31;
            double m32 = -value.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Adds each matrix element in value1 with its corresponding element in value2.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>The matrix containing the summed values.</returns>
        public static Matrix3x2 operator +(Matrix3x2 value1, Matrix3x2 value2)
        {
            double m11 = value1.M11 + value2.M11;
            double m12 = value1.M12 + value2.M12;
            double m21 = value1.M21 + value2.M21;
            double m22 = value1.M22 + value2.M22;
            double m31 = value1.M31 + value2.M31;
            double m32 = value1.M32 + value2.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Subtracts each matrix element in value2 from its corresponding element in value1.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>The matrix containing the resulting values.</returns>
        public static Matrix3x2 operator -(Matrix3x2 value1, Matrix3x2 value2)
        {
            double m11 = value1.M11 - value2.M11;
            double m12 = value1.M12 - value2.M12;
            double m21 = value1.M21 - value2.M21;
            double m22 = value1.M22 - value2.M22;
            double m31 = value1.M31 - value2.M31;
            double m32 = value1.M32 - value2.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Multiplies two matrices together and returns the resulting matrix.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>The product matrix.</returns>
        public static Matrix3x2 operator *(Matrix3x2 value1, Matrix3x2 value2)
        {
            // First row
            double m11 = (value1.M11 * value2.M11) + (value1.M12 * value2.M21);
            double m12 = (value1.M11 * value2.M12) + (value1.M12 * value2.M22);

            // Second row
            double m21 = (value1.M21 * value2.M11) + (value1.M22 * value2.M21);
            double m22 = (value1.M21 * value2.M12) + (value1.M22 * value2.M22);

            // Third row
            double m31 = (value1.M31 * value2.M11) + (value1.M32 * value2.M21) + value2.M31;
            double m32 = (value1.M31 * value2.M12) + (value1.M32 * value2.M22) + value2.M32;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Scales all elements in a matrix by the given scalar factor.
        /// </summary>
        /// <param name="value1">The source matrix.</param>
        /// <param name="value2">The scaling value to use.</param>
        /// <returns>The resulting matrix.</returns>
        public static Matrix3x2 operator *(Matrix3x2 value1, float value2)
        {
            double m11 = value1.M11 * value2;
            double m12 = value1.M12 * value2;
            double m21 = value1.M21 * value2;
            double m22 = value1.M22 * value2;
            double m31 = value1.M31 * value2;
            double m32 = value1.M32 * value2;

            return new Matrix3x2(m11, m12, m21, m22, m31, m32);
        }

        /// <summary>
        /// Returns a boolean indicating whether the given matrices are equal.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>True if the matrices are equal; False otherwise.</returns>
        public static bool operator ==(Matrix3x2 value1, Matrix3x2 value2)
        {
            return value1.M11 == value2.M11 &&
                    value1.M22 == value2.M22 && // Check diagonal element first for early out.
                    value1.M12 == value2.M12 &&
                    value1.M21 == value2.M21 &&
                    value1.M31 == value2.M31 &&
                    value1.M32 == value2.M32;
        }

        /// <summary>
        /// Returns a boolean indicating whether the given matrices are not equal.
        /// </summary>
        /// <param name="value1">The first source matrix.</param>
        /// <param name="value2">The second source matrix.</param>
        /// <returns>True if the matrices are not equal; False if they are equal.</returns>
        public static bool operator !=(Matrix3x2 value1, Matrix3x2 value2)
        {
            return value1.M11 != value2.M11 ||
                value1.M12 != value2.M12 ||
                    value1.M21 != value2.M21 ||
                    value1.M22 != value2.M22 ||
                    value1.M31 != value2.M31 ||
                    value1.M32 != value2.M32;
        }
        #endregion Operators


    }
}
