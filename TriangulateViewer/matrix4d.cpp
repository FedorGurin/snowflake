#include "matrix4d.h"

Matrix4D::Matrix4D(const double *values)
{
	for (int row = 0; row < 4; ++row)
	    for (int col = 0; col < 4; ++col)
	        m[col][row] = values[row * 4 + col];
	flagBits = General;
}

Matrix4D::Matrix4D(const float *values)
{
	for (int row = 0; row < 4; ++row)
	    for (int col = 0; col < 4; ++col)
	        m[col][row] = values[row * 4 + col];
	flagBits = General;
}

Matrix4D::Matrix4D(const double *values, int cols, int rows)
{
	for (int col = 0; col < 4; ++col) {
	    for (int row = 0; row < 4; ++row) {
	        if (col < cols && row < rows)
	            m[col][row] = values[col * rows + row];
	        else if (col == row)
	            m[col][row] = 1.0f;
	        else
	            m[col][row] = 0.0f;
	    }
	}
	flagBits = General;
}

// The 4x4 matrix inverse algorithm is based on that described at:
// http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q24
// Some optimization has been done to avoid making copies of 3x3
// sub-matrices and to unroll the loops.

// Calculate the determinant of a 3x3 sub-matrix.
//     | A B C |
// M = | D E F |   det(M) = A * (EI - HF) - B * (DI - GF) + C * (DH - GE)
//     | G H I |
static inline double matrixDet3(const double m[4][4], int col0, int col1, int col2, int row0, int row1, int row2)
{
	return m[col0][row0] *
	            (m[col1][row1] * m[col2][row2] -
	             m[col1][row2] * m[col2][row1]) -
	       m[col1][row0] *
	            (m[col0][row1] * m[col2][row2] -
	             m[col0][row2] * m[col2][row1]) +
	       m[col2][row0] *
	            (m[col0][row1] * m[col1][row2] -
	             m[col0][row2] * m[col1][row1]);
}

// Calculate the determinant of a 4x4 matrix.
static inline double matrixDet4(const double m[4][4])
{
	double det;
	det  = m[0][0] * matrixDet3(m, 1, 2, 3, 1, 2, 3);
	det -= m[1][0] * matrixDet3(m, 0, 2, 3, 1, 2, 3);
	det += m[2][0] * matrixDet3(m, 0, 1, 3, 1, 2, 3);
	det -= m[3][0] * matrixDet3(m, 0, 1, 2, 1, 2, 3);
	return det;
}

/*!
    Returns the determinant of this matrix.
*/
double Matrix4D::determinant() const
{
	return double(matrixDet4(m));
}

Matrix4D Matrix4D::inverted(bool *invertible) const
{
	// Handle some of the easy cases first.
	if (flagBits == Identity) {
	    if (invertible)
	        *invertible = true;
	    return Matrix4D();
	} else if (flagBits == Translation) {
	    Matrix4D inv;
	    inv.m[3][0] = -m[3][0];
	    inv.m[3][1] = -m[3][1];
	    inv.m[3][2] = -m[3][2];
	    inv.flagBits = Translation;
	    if (invertible)
	        *invertible = true;
	    return inv;
	} else if (flagBits == Rotation || flagBits == (Rotation | Translation)) {
	    if (invertible)
	        *invertible = true;
	    return orthonormalInverse();
	}
	
	Matrix4D inv(1); // The "1" says to not load the identity.
	
	double det = matrixDet4(m);
	if (det == 0.0f) {
	    if (invertible)
	        *invertible = false;
	    return Matrix4D();
	}
	det = 1.0f / det;
	
	inv.m[0][0] =  matrixDet3(m, 1, 2, 3, 1, 2, 3) * det;
	inv.m[0][1] = -matrixDet3(m, 0, 2, 3, 1, 2, 3) * det;
	inv.m[0][2] =  matrixDet3(m, 0, 1, 3, 1, 2, 3) * det;
	inv.m[0][3] = -matrixDet3(m, 0, 1, 2, 1, 2, 3) * det;
	inv.m[1][0] = -matrixDet3(m, 1, 2, 3, 0, 2, 3) * det;
	inv.m[1][1] =  matrixDet3(m, 0, 2, 3, 0, 2, 3) * det;
	inv.m[1][2] = -matrixDet3(m, 0, 1, 3, 0, 2, 3) * det;
	inv.m[1][3] =  matrixDet3(m, 0, 1, 2, 0, 2, 3) * det;
	inv.m[2][0] =  matrixDet3(m, 1, 2, 3, 0, 1, 3) * det;
	inv.m[2][1] = -matrixDet3(m, 0, 2, 3, 0, 1, 3) * det;
	inv.m[2][2] =  matrixDet3(m, 0, 1, 3, 0, 1, 3) * det;
	inv.m[2][3] = -matrixDet3(m, 0, 1, 2, 0, 1, 3) * det;
	inv.m[3][0] = -matrixDet3(m, 1, 2, 3, 0, 1, 2) * det;
	inv.m[3][1] =  matrixDet3(m, 0, 2, 3, 0, 1, 2) * det;
	inv.m[3][2] = -matrixDet3(m, 0, 1, 3, 0, 1, 2) * det;
	inv.m[3][3] =  matrixDet3(m, 0, 1, 2, 0, 1, 2) * det;
	
	if (invertible)
	    *invertible = true;
	return inv;
}

/*!
    Returns the normal matrix corresponding to this 4x4 transformation.
    The normal matrix is the transpose of the inverse of the top-left
    3x3 part of this 4x4 matrix.  If the 3x3 sub-matrix is not invertible,
    this function returns the identity.

    \sa inverted()
*/
QMatrix3x3 Matrix4D::normalMatrix() const
{
	QMatrix3x3 inv;
	
	// Handle the simple cases first.
	if (flagBits == Identity || flagBits == Translation) {
	    return inv;
	} else if (flagBits == Scale || flagBits == (Translation | Scale)) {
	    if (m[0][0] == 0.0f || m[1][1] == 0.0f || m[2][2] == 0.0f)
	        return inv;
	    inv.data()[0] = 1.0f / m[0][0];
	    inv.data()[4] = 1.0f / m[1][1];
	    inv.data()[8] = 1.0f / m[2][2];
	    return inv;
	}
	
	double det = matrixDet3(m, 0, 1, 2, 0, 1, 2);
	if (det == 0.0f)
	    return inv;
	det = 1.0f / det;
	
	float *invm = inv.data();
	
	// Invert and transpose in a single step.
	invm[0 + 0 * 3] =  (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * det;
	invm[1 + 0 * 3] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) * det;
	invm[2 + 0 * 3] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * det;
	invm[0 + 1 * 3] = -(m[0][1] * m[2][2] - m[2][1] * m[0][2]) * det;
	invm[1 + 1 * 3] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * det;
	invm[2 + 1 * 3] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) * det;
	invm[0 + 2 * 3] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * det;
	invm[1 + 2 * 3] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) * det;
	invm[2 + 2 * 3] =  (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * det;
	
	return inv;
}

/*!
    Returns this matrix, transposed about its diagonal.
*/
Matrix4D Matrix4D::transposed() const
{
	Matrix4D result(1); // The "1" says to not load the identity.
	for (int row = 0; row < 4; ++row) {
	    for (int col = 0; col < 4; ++col) {
	        result.m[col][row] = m[row][col];
	    }
	}
	return result;
}

Matrix4D Matrix4D::getOrthogonalMatrix() const
{
	Matrix4D mat(*this);
	Matrix4D it = inverted().transposed();
	while (!qFuzzyCompare(it, mat))
	{
		mat = 0.5*(mat + it);
		it = mat.inverted().transposed();
	}
	return mat;
}

void Matrix4D::rotate(float angle, double x, double y, double z)
{
	if (angle == 0.0f)
	    return;
	Matrix4D m(1); // The "1" says to not load the identity.
	double c, s, ic;
	if (angle == 90.0f || angle == -270.0f) {
	    s = 1.0f;
	    c = 0.0f;
	} else if (angle == -90.0f || angle == 270.0f) {
	    s = -1.0f;
	    c = 0.0f;
	} else if (angle == 180.0f || angle == -180.0f) {
	    s = 0.0f;
	    c = -1.0f;
	} else {
	    double a = angle * M_PI / 180.0f;
	    c = cos(a);
	    s = sin(a);
	}
	bool quick = false;
	if (x == 0.0f) {
	    if (y == 0.0f) {
	        if (z != 0.0f) {
	            // Rotate around the Z axis.
	            m.setToIdentity();
	            m.m[0][0] = c;
	            m.m[1][1] = c;
	            if (z < 0.0f) {
	                m.m[1][0] = s;
	                m.m[0][1] = -s;
	            } else {
	                m.m[1][0] = -s;
	                m.m[0][1] = s;
	            }
	            m.flagBits = General;
	            quick = true;
	        }
	    } else if (z == 0.0f) {
	        // Rotate around the Y axis.
	        m.setToIdentity();
	        m.m[0][0] = c;
	        m.m[2][2] = c;
	        if (y < 0.0f) {
	            m.m[2][0] = -s;
	            m.m[0][2] = s;
	        } else {
	            m.m[2][0] = s;
	            m.m[0][2] = -s;
	        }
	        m.flagBits = General;
	        quick = true;
	    }
	} else if (y == 0.0f && z == 0.0f) {
	    // Rotate around the X axis.
	    m.setToIdentity();
	    m.m[1][1] = c;
	    m.m[2][2] = c;
	    if (x < 0.0f) {
	        m.m[2][1] = s;
	        m.m[1][2] = -s;
	    } else {
	        m.m[2][1] = -s;
	        m.m[1][2] = s;
	    }
	    m.flagBits = General;
	    quick = true;
	}
	if (!quick) 
	{
		if (!qFuzzyIsNull(Vector3D(x,y,z)))
		{
			double len = x * x + y * y + z * z;
			if (!qFuzzyIsNull(len - 1.0))
			{
		    len = 1.0/sqrt(len);
		    x *= len;
		    y *= len;
		    z *= len;
			}
		}
		ic = 1.0f - c;
		m.m[0][0] = x * x * ic + c;
		m.m[1][0] = x * y * ic - z * s;
		m.m[2][0] = x * z * ic + y * s;
		m.m[3][0] = 0.0f;
		m.m[0][1] = y * x * ic + z * s;
		m.m[1][1] = y * y * ic + c;
		m.m[2][1] = y * z * ic - x * s;
		m.m[3][1] = 0.0f;
		m.m[0][2] = x * z * ic - y * s;
		m.m[1][2] = y * z * ic + x * s;
		m.m[2][2] = z * z * ic + c;
		m.m[3][2] = 0.0f;
		m.m[0][3] = 0.0f;
		m.m[1][3] = 0.0f;
		m.m[2][3] = 0.0f;
		m.m[3][3] = 1.0f;
	}
	int flags = flagBits;
	*this *= m;
	if (flags != Identity)
	    flagBits = flags | Rotation;
	else
	    flagBits = Rotation;
}
//
Matrix4D Matrix4D::perspective(double angle, double aspect, double nearPlane, double farPlane)
{
	// Bail out if the projection volume is zero-sized.
	if (nearPlane == farPlane || aspect == 0.0f)
	    return Matrix4D();
	
	// Construct the projection.
	Matrix4D m(1);
	double radians = (angle / 2.0f) * M_PI / 180.0f;
	double sine = sin(radians);
	if (sine == 0.0f)
	    return Matrix4D();
	double cotan = cos(radians) / sine;
	double clip = farPlane - nearPlane;
	m.m[0][0] = cotan / aspect;
	m.m[1][0] = 0.0f;
	m.m[2][0] = 0.0f;
	m.m[3][0] = 0.0f;
	m.m[0][1] = 0.0f;
	m.m[1][1] = cotan;
	m.m[2][1] = 0.0f;
	m.m[3][1] = 0.0f;
	m.m[0][2] = 0.0f;
	m.m[1][2] = 0.0f;
	m.m[2][2] = -(nearPlane + farPlane) / clip;
	m.m[3][2] = -(2.0f * nearPlane * farPlane) / clip;
	m.m[0][3] = 0.0f;
	m.m[1][3] = 0.0f;
	m.m[2][3] = -1.0f;
	m.m[3][3] = 0.0f;
	
	return m;
}

/*!
    Retrieves the 16 items in this matrix and copies them to \a values
    in row-major order.
*/
QVarLengthArray<float, 16> Matrix4D::toFloatArray() const 
{ 
	QVarLengthArray<float, 16> result(16);
  for (int index = 0; index < 16; ++index) 
		result[index] =  (*m)[index];
	return result;
}

void Matrix4D::copyDataTo(double *values) const
{
	for (int row = 0; row < 4; ++row)
		for (int col = 0; col < 4; ++col)
			values[row * 4 + col] = m[col][row];
}

void Matrix4D::copyDataTo(float *values) const
{
	for (int row = 0; row < 4; ++row)
		for (int col = 0; col < 4; ++col)
			values[row * 4 + col] = m[col][row];
}

QDebug operator<<(QDebug dbg, const Matrix4D &m)
{
	// Create a string that represents the matrix type.
	QByteArray bits;
	if ((m.flagBits & Matrix4D::Identity) != 0)
	    bits += "Identity,";
	if ((m.flagBits & Matrix4D::General) != 0)
	    bits += "General,";
	if ((m.flagBits & Matrix4D::Translation) != 0)
	    bits += "Translation,";
	if ((m.flagBits & Matrix4D::Scale) != 0)
	    bits += "Scale,";
	if ((m.flagBits & Matrix4D::Rotation) != 0)
	    bits += "Rotation,";
	if (bits.size() > 0)
	    bits = bits.left(bits.size() - 1);
	
	// Output in row-major order because it is more human-readable.
	dbg.nospace() << "Matrix4D(type:" << bits.constData() << endl
	    << qSetFieldWidth(10)
	    << m(0, 0) << m(0, 1) << m(0, 2) << m(0, 3) << endl
	    << m(1, 0) << m(1, 1) << m(1, 2) << m(1, 3) << endl
	    << m(2, 0) << m(2, 1) << m(2, 2) << m(2, 3) << endl
	    << m(3, 0) << m(3, 1) << m(3, 2) << m(3, 3) << endl
	    << qSetFieldWidth(0) << ')';
	return dbg.space();
}
/*!
    \fn QDataStream &operator<<(QDataStream &stream, const Matrix4D &matrix)
    \relates Matrix4D

    Writes the given \a matrix to the given \a stream and returns a
    reference to the stream.

    \sa {Serializing Qt Data Types}
*/
QDataStream &operator<<(QDataStream &stream, const Matrix4D &matrix)
{
	for (int row = 0; row < 4; ++row)
		for (int col = 0; col < 4; ++col)
			stream << matrix(row, col);
	return stream;
}

/*!
    \fn QDataStream &operator>>(QDataStream &stream, Matrix4D &matrix)
    \relates Matrix4D

    Reads a 4x4 matrix from the given \a stream into the given \a matrix
    and returns a reference to the stream.

    \sa {Serializing Qt Data Types}
*/
QDataStream &operator>>(QDataStream &stream, Matrix4D &matrix)
{
	double x;
	for (int row = 0; row < 4; ++row) {
		for (int col = 0; col < 4; ++col) {
		    stream >> x;
		    matrix(row, col) = double(x);
		}
	}
	matrix.optimize();
	return stream;
}
