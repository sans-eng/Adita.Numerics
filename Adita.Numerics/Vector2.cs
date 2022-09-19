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
using System.Text;

namespace Adita.Numerics
{
    /// <summary>
    /// Represents a 2d vector.
    /// </summary>
    public struct Vector2 : IEquatable<Vector2>, IFormattable
    {
        #region Constructors
        /// <summary>
        /// Creates a vector whose elements are all the single specified value.
        /// </summary>
        /// <param name="value">The element to fill the vector with.</param>
        public Vector2(double value) : this(value, value) { }

        /// <summary>
        /// Creates a vector with the given individual elements.
        /// </summary>
        /// <param name="x">The X component.</param>
        /// <param name="y">The Y component.</param>
        public Vector2(double x, double y)
        {
            X = x;
            Y = y;
        }
        #endregion Constructors

        #region Public properties
        /// <summary>
        /// Gets the X component of the vector.
        /// </summary>
        public double X { get; }
        /// <summary>
        /// Gets the Y component of the vector.
        /// </summary>
        public double Y { get; }
        #endregion Public properties

        #region Public Static Properties
        /// <summary>
        /// Returns the vector (0,0).
        /// </summary>
        public static Vector2 Zero { get { return new Vector2(); } }
        /// <summary>
        /// Returns the vector (1,1).
        /// </summary>
        public static Vector2 One { get { return new Vector2(1.0d, 1.0d); } }
        /// <summary>
        /// Returns the vector (1,0).
        /// </summary>
        public static Vector2 UnitX { get { return new Vector2(1.0d, 0.0d); } }
        /// <summary>
        /// Returns the vector (0,1).
        /// </summary>
        public static Vector2 UnitY { get { return new Vector2(0.0d, 1.0d); } }
        #endregion Public Static Properties

        #region Public instance methods
       
        /// <summary>Returns the hash code for this instance.</summary>
        /// <returns>A 32-bit signed integer that is the hash code for this instance.</returns>
        public override int GetHashCode()
        {
            int hashCode = 1764330041;
            hashCode *= -1521134295 + X.GetHashCode();
            hashCode *= -1521134295 + Y.GetHashCode();
            return hashCode;
        }

        /// <summary>
        /// Returns a boolean indicating whether the given Object is equal to this Vector2 instance.
        /// </summary>
        /// <param name="obj">The Object to compare against.</param>
        /// <returns>True if the Object is equal to this Vector2; False otherwise.</returns>
        public override bool Equals(object? obj) => obj is Vector2 vector2 && Equals(vector2);

        /// <summary>
        /// Returns a String representing this Vector2 instance.
        /// </summary>
        /// <returns>The string representation.</returns>
        public override string ToString()
        {
            return ToString("G", CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Returns a String representing this Vector2 instance, using the specified format to format individual elements.
        /// </summary>
        /// <param name="format">The format of individual elements.</param>
        /// <returns>The string representation.</returns>
        public string ToString(string format)
        {
            return ToString(format, CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Returns a String representing this Vector2 instance, using the specified format to format individual elements 
        /// and the given IFormatProvider.
        /// </summary>
        /// <param name="format">The format of individual elements.</param>
        /// <param name="formatProvider">The format provider to use when formatting elements.</param>
        /// <returns>The string representation.</returns>
        public string ToString(string? format, IFormatProvider? formatProvider)
        {
            StringBuilder sb = new();
            string separator = NumberFormatInfo.GetInstance(formatProvider).NumberGroupSeparator;
            sb.Append('<');
            sb.Append(this.X.ToString(format, formatProvider));
            sb.Append(separator);
            sb.Append(' ');
            sb.Append(this.Y.ToString(format, formatProvider));
            sb.Append('>');
            return sb.ToString();
        }

        /// <summary>
        /// Returns the length of the vector.
        /// </summary>
        /// <returns>The vector's length.</returns>
        public double Length()
        {
            double ls = Dot(this, this);
            return Math.Sqrt(ls);
        }

        /// <summary>
        /// Returns the length of the vector squared. This operation is cheaper than Length().
        /// </summary>
        /// <returns>The vector's length squared.</returns>
        public double LengthSquared()
        {
            return Dot(this, this);
        }
        /// <summary>
        /// Copies the contents of the vector into the given array.
        /// </summary>
        /// <param name="array">The destination array.</param>
        public void CopyTo(double[] array)
        {
            CopyTo(array, 0);
        }

        /// <summary>
        /// Copies the contents of the vector into the given array, starting from the given index.
        /// </summary>
        /// <exception cref="ArgumentNullException">If array is null.</exception>
        /// <exception cref="RankException">If array is multidimensional.</exception>
        /// <exception cref="ArgumentOutOfRangeException">If index is greater than end of the array or index is less than zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">If number of elements in source vector is greater than those available in destination array
        /// or if there are not enough elements to copy.</exception>
        public void CopyTo(double[] array, int index)
        {
            if (array == null)
            {
                throw new ArgumentNullException(nameof(array));
            }
            if (index < 0 || index >= array.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(array), $"{nameof(index)} is a negative number or less than {nameof(array)} length");
            }
            if ((array.Length - index) < 2)
            {
                throw new ArgumentOutOfRangeException(nameof(array), $"{nameof(array)} length is less than the source.");
            }
            array[index] = X;
            array[index + 1] = Y;
        }

        /// <summary>
        /// Returns a boolean indicating whether the given Vector2 is equal to this Vector2 instance.
        /// </summary>
        /// <param name="other">The Vector2 to compare this instance to.</param>
        /// <returns>True if the other Vector2 is equal to this instance; False otherwise.</returns>
        public bool Equals(Vector2 other)
        {
            return this.X == other.X && this.Y == other.Y;
        }
        #endregion Public Instance Methods

        #region Public Static Methods
        /// <summary>
        /// Returns the Euclidean distance between the two given points.
        /// </summary>
        /// <param name="value1">The first point.</param>
        /// <param name="value2">The second point.</param>
        /// <returns>The distance.</returns>
        public static double Distance(Vector2 value1, Vector2 value2)
        {
            Vector2 difference = value1 - value2;
            double ls = Dot(difference, difference);
            return Math.Sqrt(ls);
        }
        /// <summary>
        /// Returns the Euclidean distance squared between the two given points.
        /// </summary>
        /// <param name="value1">The first point.</param>
        /// <param name="value2">The second point.</param>
        /// <returns>The distance squared.</returns>
        public static double DistanceSquared(Vector2 value1, Vector2 value2)
        {
            Vector2 difference = value1 - value2;
            return Dot(difference, difference);
        }
        /// <summary>
        /// Returns a vector with the same direction as the given vector, but with a length of 1.
        /// </summary>
        /// <param name="value">The vector to normalize.</param>
        /// <returns>The normalized vector.</returns>
        public static Vector2 Normalize(Vector2 value)
        {
            double length = value.Length();
            return value / length;
        }
        /// <summary>
        /// Returns the reflection of a vector off a surface that has the specified normal.
        /// </summary>
        /// <param name="vector">The source vector.</param>
        /// <param name="normal">The normal of the surface being reflected off.</param>
        /// <returns>The reflected vector.</returns>
        public static Vector2 Reflect(Vector2 vector, Vector2 normal)
        {
            double dot = Dot(vector, normal);
            return vector - (2 * dot * normal);
        }
        /// <summary>
        /// Restricts a vector between a min and max value.
        /// </summary>
        /// <param name="value1">The source vector.</param>
        /// <param name="min">The minimum value.</param>
        /// <param name="max">The maximum value.</param>
        public static Vector2 Clamp(Vector2 value1, Vector2 min, Vector2 max)
        {
            // This compare order is very important!!!
            // We must follow HLSL behavior in the case user specified min value is bigger than max value.
            double x = value1.X;
            x = (x > max.X) ? max.X : x;
            x = (x < min.X) ? min.X : x;

            double y = value1.Y;
            y = (y > max.Y) ? max.Y : y;
            y = (y < min.Y) ? min.Y : y;

            return new Vector2(x, y);
        }
        /// <summary>
        /// Linearly interpolates between two vectors based on the given weighting.
        /// </summary>
        /// <param name="value1">The first source vector.</param>
        /// <param name="value2">The second source vector.</param>
        /// <param name="amount">Value between 0 and 1 indicating the weight of the second source vector.</param>
        /// <returns>The interpolated vector.</returns>
        public static Vector2 Lerp(Vector2 value1, Vector2 value2, double amount)
        {
            return new Vector2(
                value1.X + ((value2.X - value1.X) * amount),
                value1.Y + ((value2.Y - value1.Y) * amount));
        }
        /// <summary>
        /// Transforms a vector by the given matrix.
        /// </summary>
        /// <param name="position">The source vector.</param>
        /// <param name="matrix">The transformation matrix.</param>
        /// <returns>The transformed vector.</returns>
        public static Vector2 Transform(Vector2 position, Matrix3x2 matrix)
        {
            return new Vector2(
                (position.X * matrix.M11) + (position.Y * matrix.M21) + matrix.M31,
                (position.X * matrix.M12) + (position.Y * matrix.M22) + matrix.M32);
        }
        /// <summary>
        /// Transforms a vector normal by the given matrix.
        /// </summary>
        /// <param name="normal">The source vector.</param>
        /// <param name="matrix">The transformation matrix.</param>
        /// <returns>The transformed vector.</returns>
        public static Vector2 TransformNormal(Vector2 normal, Matrix3x2 matrix)
        {
            return new Vector2(
                (normal.X * matrix.M11) + (normal.Y * matrix.M21),
                (normal.X * matrix.M12) + (normal.Y * matrix.M22));
        }
        /// <summary>
        /// Returns the dot product of two vectors.
        /// </summary>
        /// <param name="value1">The first vector.</param>
        /// <param name="value2">The second vector.</param>
        /// <returns>The dot product.</returns>
        public static double Dot(Vector2 value1, Vector2 value2)
        {
            return (value1.X * value2.X) +
                   (value1.Y * value2.Y);
        }
        /// <summary>
        /// Returns a vector whose elements are the minimum of each of the pairs of elements in the two source vectors.
        /// </summary>
        /// <param name="value1">The first source vector.</param>
        /// <param name="value2">The second source vector.</param>
        /// <returns>The minimized vector.</returns>
        public static Vector2 Min(Vector2 value1, Vector2 value2)
        {
            return new Vector2(
                (value1.X < value2.X) ? value1.X : value2.X,
                (value1.Y < value2.Y) ? value1.Y : value2.Y);
        }
        /// <summary>
        /// Returns a vector whose elements are the maximum of each of the pairs of elements in the two source vectors
        /// </summary>
        /// <param name="value1">The first source vector</param>
        /// <param name="value2">The second source vector</param>
        /// <returns>The maximized vector</returns>
        public static Vector2 Max(Vector2 value1, Vector2 value2)
        {
            return new Vector2(
                (value1.X > value2.X) ? value1.X : value2.X,
                (value1.Y > value2.Y) ? value1.Y : value2.Y);
        }
        /// <summary>
        /// Returns a vector whose elements are the absolute values of each of the source vector's elements.
        /// </summary>
        /// <param name="value">The source vector.</param>
        /// <returns>The absolute value vector.</returns>
        public static Vector2 Abs(Vector2 value)
        {
            return new Vector2(Math.Abs(value.X), Math.Abs(value.Y));
        }
        /// <summary>
        /// Returns a vector whose elements are the square root of each of the source vector's elements.
        /// </summary>
        /// <param name="value">The source vector.</param>
        /// <returns>The square root vector.</returns>
        public static Vector2 SquareRoot(Vector2 value)
        {
            return new Vector2(Math.Sqrt(value.X), Math.Sqrt(value.Y));
        }
        #endregion Public Static Methods

        #region Public operator methods
        // all the below methods should be inlined as they are 
        // implemented over JIT intrinsics

        /// <summary>
        /// Adds two vectors together.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The summed vector.</returns>
        public static Vector2 Add(Vector2 left, Vector2 right)
        {
            return left + right;
        }
        /// <summary>
        /// Subtracts the second vector from the first.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The difference vector.</returns>
        public static Vector2 Subtract(Vector2 left, Vector2 right)
        {
            return left - right;
        }
        /// <summary>
        /// Multiplies two vectors together.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The product vector.</returns>
        public static Vector2 Multiply(Vector2 left, Vector2 right)
        {
            return left * right;
        }
        /// <summary>
        /// Multiplies a vector by the given scalar.
        /// </summary>
        /// <param name="left">The source vector.</param>
        /// <param name="right">The scalar value.</param>
        /// <returns>The scaled vector.</returns>
        public static Vector2 Multiply(Vector2 left, double right)
        {
            return left * right;
        }
        /// <summary>
        /// Multiplies a vector by the given scalar.
        /// </summary>
        /// <param name="left">The scalar value.</param>
        /// <param name="right">The source vector.</param>
        /// <returns>The scaled vector.</returns>
        public static Vector2 Multiply(double left, Vector2 right)
        {
            return left * right;
        }
        /// <summary>
        /// Divides the first vector by the second.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The vector resulting from the division.</returns>
        public static Vector2 Divide(Vector2 left, Vector2 right)
        {
            return left / right;
        }
        /// <summary>
        /// Divides the vector by the given scalar.
        /// </summary>
        /// <param name="left">The source vector.</param>
        /// <param name="divisor">The scalar value.</param>
        /// <returns>The result of the division.</returns>
        public static Vector2 Divide(Vector2 left, double divisor)
        {
            return left / divisor;
        }
        /// <summary>
        /// Negates a given vector.
        /// </summary>
        /// <param name="value">The source vector.</param>
        /// <returns>The negated vector.</returns>
        public static Vector2 Negate(Vector2 value)
        {
            return -value;
        }
        #endregion Public operator methods

        #region Public Operators
        /// <summary>
        /// Adds two vectors together.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The summed vector.</returns>
        public static Vector2 operator +(Vector2 left, Vector2 right)
        {
            return new Vector2(left.X + right.X, left.Y + right.Y);
        }

        /// <summary>
        /// Subtracts the second vector from the first.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The difference vector.</returns>
        public static Vector2 operator -(Vector2 left, Vector2 right)
        {
            return new Vector2(left.X - right.X, left.Y - right.Y);
        }

        /// <summary>
        /// Multiplies two vectors together.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The product vector.</returns>
        public static Vector2 operator *(Vector2 left, Vector2 right)
        {
            return new Vector2(left.X * right.X, left.Y * right.Y);
        }

        /// <summary>
        /// Multiplies a vector by the given scalar.
        /// </summary>
        /// <param name="left">The scalar value.</param>
        /// <param name="right">The source vector.</param>
        /// <returns>The scaled vector.</returns>
        public static Vector2 operator *(double left, Vector2 right)
        {
            return new Vector2(left, left) * right;
        }

        /// <summary>
        /// Multiplies a vector by the given scalar.
        /// </summary>
        /// <param name="left">The source vector.</param>
        /// <param name="right">The scalar value.</param>
        /// <returns>The scaled vector.</returns>
        public static Vector2 operator *(Vector2 left, double right)
        {
            return left * new Vector2(right, right);
        }

        /// <summary>
        /// Divides the first vector by the second.
        /// </summary>
        /// <param name="left">The first source vector.</param>
        /// <param name="right">The second source vector.</param>
        /// <returns>The vector resulting from the division.</returns>
        public static Vector2 operator /(Vector2 left, Vector2 right)
        {
            return new Vector2(left.X / right.X, left.Y / right.Y);
        }

        /// <summary>
        /// Divides the vector by the given scalar.
        /// </summary>
        /// <param name="value1">The source vector.</param>
        /// <param name="value2">The scalar value.</param>
        /// <returns>The result of the division.</returns>
        public static Vector2 operator /(Vector2 value1, double value2)
        {
            double invDiv = 1.0d / value2;
            return new Vector2(
                value1.X * invDiv,
                value1.Y * invDiv);
        }

        /// <summary>
        /// Negates a given vector.
        /// </summary>
        /// <param name="value">The source vector.</param>
        /// <returns>The negated vector.</returns>
        public static Vector2 operator -(Vector2 value)
        {
            return Zero - value;
        }

        /// <summary>
        /// Returns a boolean indicating whether the two given vectors are equal.
        /// </summary>
        /// <param name="left">The first vector to compare.</param>
        /// <param name="right">The second vector to compare.</param>
        /// <returns>True if the vectors are equal; False otherwise.</returns>
        public static bool operator ==(Vector2 left, Vector2 right)
        {
            return left.Equals(right);
        }

        /// <summary>
        /// Returns a boolean indicating whether the two given vectors are not equal.
        /// </summary>
        /// <param name="left">The first vector to compare.</param>
        /// <param name="right">The second vector to compare.</param>
        /// <returns>True if the vectors are not equal; False if they are equal.</returns>
        public static bool operator !=(Vector2 left, Vector2 right)
        {
            return !(left == right);
        }
        #endregion Public Operators
    }
}
