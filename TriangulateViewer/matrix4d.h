#ifndef MATRIX4X4D_H
#define MATRIX4X4D_H

#include <QGenericMatrix>
#include <QMatrix4x4>
#include <QVarLengthArray>
#include <QRect>
#include "vector4d.h"
#include "quaterniond.h"

inline double snapAngle(double angle)
{
	if (qFuzzyIsNull((float)angle))
		return 0;
	
	double abs_angle = abs(angle);
	int sign = angle >= 0 ? 1: -1;

	if (qFuzzyIsNull((float)(abs_angle-M_PI)))
		return sign*M_PI;
	if (qFuzzyIsNull((float)(abs_angle-M_PI_2)))
		return sign*M_PI_2;

	return angle;
	//if (qFuzzyIsNull((float)(abs_angle-3*M_PI_2)))
	//	return sign*3*M_PI_2;
}

class Matrix4D
{
public:
	inline Matrix4D() { setToIdentity(); }
	explicit Matrix4D(const double *values);
	explicit Matrix4D(const float *values);
	inline Matrix4D(double m11, double m12, double m13, double m14,
	                double m21, double m22, double m23, double m24,
	                double m31, double m32, double m33, double m34,
	                double m41, double m42, double m43, double m44);
	
	template <int N, int M>
	explicit Matrix4D(const QGenericMatrix<N, M, double>& matrix);
	QMatrix4x4 toQMatrix4x4() const;// { return }
	
	Matrix4D(const double *values, int cols, int rows);
	//Matrix4D(const QTransform& transform);
	//Matrix4D(const QMatrix& matrix);
	
	inline const double& operator()(int row, int column) const;
	inline double& operator()(int row, int column);
	
	inline Vector4D column(int index) const;
	inline void setColumn(int index, const Vector4D& value);
	
	inline Vector4D row(int index) const;
	inline void setRow(int index, const Vector4D& value);
	
	inline bool isIdentity() const;
	inline void setToIdentity();
	
	inline void fill(double value);
	
	double determinant() const;
	Matrix4D inverted(bool *invertible = 0) const;
	Matrix4D transposed() const;
	QMatrix3x3 normalMatrix() const;
	//inline QGenericMatrix<3,3,double> topSubMatrix() const;
	
	inline Matrix4D& operator+=(const Matrix4D& other);
	inline Matrix4D& operator-=(const Matrix4D& other);
	inline Matrix4D& operator*=(const Matrix4D& other);
	inline Matrix4D& operator*=(double factor);
	inline Matrix4D& operator/=(float divisor);
	inline bool operator==(const Matrix4D& other) const;
	inline bool operator!=(const Matrix4D& other) const;
	
	friend inline Matrix4D operator+(const Matrix4D& m1, const Matrix4D& m2);
	friend inline Matrix4D operator-(const Matrix4D& m1, const Matrix4D& m2);
	friend inline Matrix4D operator*(const Matrix4D& m1, const Matrix4D& m2);
	friend inline Vector3D operator*(const Matrix4D& matrix, const Vector3D& vector);
	friend inline Vector3D operator*(const Vector3D& vector, const Matrix4D& matrix);
	friend inline Vector4D operator*(const Matrix4D& matrix, const Vector4D& vector);
	friend inline Vector4D operator*(const Vector4D& vector, const Matrix4D& matrix);
				 
	friend inline QVector3D operator*(const Matrix4D& matrix, const QVector3D& vector) { return (matrix*Vector3D(vector)).toQVector3D(); }
	friend inline QVector3D operator*(const QVector3D& vector, const Matrix4D& matrix) { return (Vector3D(vector)*matrix).toQVector3D(); }
	friend inline QVector4D operator*(const Matrix4D& matrix, const QVector4D& vector) { return (matrix*Vector4D(vector)).toQVector4D(); }
	friend inline QVector4D operator*(const QVector4D& vector, const Matrix4D& matrix) { return (Vector4D(vector)*matrix).toQVector4D(); }
				 
	friend inline QPoint operator*(const QPoint& point, const Matrix4D& matrix);
	friend inline QPointF operator*(const QPointF& point, const Matrix4D& matrix);
	friend inline Matrix4D operator-(const Matrix4D& matrix);
	friend inline QPoint operator*(const Matrix4D& matrix, const QPoint& point);
	friend inline QPointF operator*(const Matrix4D& matrix, const QPointF& point);
				 
	friend inline Matrix4D operator*(double factor, const Matrix4D& matrix);
	friend inline Matrix4D operator*(const Matrix4D& matrix, double factor);
	friend inline Matrix4D operator/(const Matrix4D& matrix, float divisor);
	
	friend inline bool qFuzzyCompare(const Matrix4D& m1, const Matrix4D& m2);
	
	void scale(const Vector3D& vector) { scale(vector.x, vector.y, vector.z); }
	void translate(const Vector3D& vector) { translate(vector.x, vector.y, vector.z); }
	void rotate(float angle, const Vector3D& vector) { rotate(angle, vector.x, vector.y, vector.z); }
	inline void scale(double x, double y);
	inline void scale(double x, double y, double z);
	inline void scale(float factor);
	inline void translate(double x, double y);
	inline void translate(double x, double y, double z);
	void rotate(float angle, double x, double y, double z = 0.0f);
	inline void rotate(const QuaternionD& quaternion);
	
	static Matrix4D ortho(double left, double right, double bottom, double top, double nearPlane, double farPlane);
	static Matrix4D frustum(double left, double right, double bottom, double top, double nearPlane, double farPlane);
	static Matrix4D perspective(double angle, double aspect, double nearPlane, double farPlane);
	static Matrix4D lookAt(const Vector3D& eye, const Vector3D& center, const Vector3D& up);
	
	//static Matrix4D scaleMat(double x, double y, double z);
	//static Matrix4D scaleMat(double f);
	
	void flipCoordinates();
	
	void copyDataTo(double *values) const;
	void copyDataTo(float *values) const;
	QVarLengthArray<float, 16> toFloatArray() const;
	
	//QMatrix toAffine() const;
	//QTransform toTransform() const;
	//QTransform toTransform(double distanceToPlane) const;
	
	//QPoint map(const QPoint& point) const;
	QPointF map(const QPointF& point) const { return *this * point;	}
	QVector4D map(const QVector4D& point) const { return *this * point; }
	QVector3D map(const QVector3D& point) const { return *this * point; }
	Vector4D map(const Vector4D& point) const { return *this * point; }
	Vector3D map(const Vector3D& point) const { return *this * point; }
    
	inline QRectF mapRect(const QRectF& rect) const;
	
	template <int N, int M>
	QGenericMatrix<N, M, double> toGenericMatrix() const;
	
	double *data() { flagBits = General; return *m; } //!!!
	const double *data() const { return *m; }
	const double *constData() const { return *m; }
	
	void optimize();

	friend QDebug operator<<(QDebug dbg, const Matrix4D &m);
	
	double determinant3x3() const
	{
		return Vector3D::dotProduct(Vector3D(m[0]), Vector3D::crossProduct(Vector3D(m[1]), Vector3D(m[2])));
	}
	
	Vector3D getScaleVec() const
	{
		return Vector3D(Vector3D(m[0]).length(),
										Vector3D(m[1]).length(),
										Vector3D(m[2]).length());
	}
	
	Matrix4D getScaleMat() const
	{
		Matrix4D res;
		res.scale(getScaleVec());
		return res;
	}
	
	Matrix4D getInvertedScaleMat() const
	{
		Vector3D scaleVector = getScaleVec();
		Matrix4D res;
		res.scale(1.0/scaleVector.x, 1.0/scaleVector.y, 1.0/scaleVector.z);
		return res;
	}
	
	Vector3D extractEulerAngles() const
	{
		double rotx = atan2(m[1][2],m[2][2]); 
		double roty = atan2(-m[0][2], sqrt(m[0][1]*m[0][1] + m[0][0]*m[0][0]));
		double rotz = atan2(m[0][1], m[0][0]);
		return Vector3D(snapAngle(rotx), snapAngle(roty), snapAngle(rotz));
	}
	
	Matrix4D getOrthogonalMatrix() const;
	
private:
	double m[4][4];          // Column-major order to match OpenGL.
	int flagBits;           // Flag bits from the enum below.
	
	enum {
	    Identity        = 0x0001,   // Identity matrix
	    General         = 0x0002,   // General matrix, unknown contents
	    Translation     = 0x0004,   // Contains a simple translation
	    Scale           = 0x0008,   // Contains a simple scale
	    Rotation        = 0x0010    // Contains a simple rotation
	};
	
	// Construct without initializing identity matrix.
	Matrix4D(int) { flagBits = General; }
	
	Matrix4D orthonormalInverse() const;
};

QDataStream &operator<<(QDataStream &stream, const Matrix4D &matrix);

QDataStream &operator>>(QDataStream &stream, Matrix4D &matrix);

inline Matrix4D::Matrix4D(double m11, double m12, double m13, double m14,
													double m21, double m22, double m23, double m24,
													double m31, double m32, double m33, double m34,
													double m41, double m42, double m43, double m44)
{
	m[0][0] = m11; m[0][1] = m21; m[0][2] = m31; m[0][3] = m41;
	m[1][0] = m12; m[1][1] = m22; m[1][2] = m32; m[1][3] = m42;
	m[2][0] = m13; m[2][1] = m23; m[2][2] = m33; m[2][3] = m43;
	m[3][0] = m14; m[3][1] = m24; m[3][2] = m34; m[3][3] = m44;
	flagBits = General;
}

template <int N, int M>
Matrix4D::Matrix4D(const QGenericMatrix<N, M, double>& matrix)
{
	const double *values = matrix.constData();
	for (int matrixCol = 0; matrixCol < 4; ++matrixCol) {
		for (int matrixRow = 0; matrixRow < 4; ++matrixRow) {
			if (matrixCol < N && matrixRow < M)
				m[matrixCol][matrixRow] = values[matrixCol * M + matrixRow];
			else if (matrixCol == matrixRow)
				m[matrixCol][matrixRow] = 1.0f;
			else
				m[matrixCol][matrixRow] = 0.0f;
		}
	}
	flagBits = General;
}

template <int N, int M>
QGenericMatrix<N, M, double> Matrix4D::toGenericMatrix() const
{
	QGenericMatrix<N, M, double> result;
	double *values = result.data();
	for (int matrixCol = 0; matrixCol < N; ++matrixCol) {
		for (int matrixRow = 0; matrixRow < M; ++matrixRow) {
			if (matrixCol < 4 && matrixRow < 4)
				values[matrixCol * M + matrixRow] = m[matrixCol][matrixRow];
			else if (matrixCol == matrixRow)
				values[matrixCol * M + matrixRow] = 1.0f;
			else
				values[matrixCol * M + matrixRow] = 0.0f;
		}
	}
	return result;
}

inline QMatrix4x4 Matrix4D::toQMatrix4x4() const
{
	QMatrix4x4 result;
	float *values = result.data();
	const double *data = constData();
	for (int i = 0; i < 16; ++i)
		values[i] = data[i];
	return result;
}

inline const double& Matrix4D::operator()(int aRow, int aColumn) const
{
	Q_ASSERT(aRow >= 0 && aRow < 4 && aColumn >= 0 && aColumn < 4);
	return m[aColumn][aRow];
}

inline double& Matrix4D::operator()(int aRow, int aColumn)
{
	Q_ASSERT(aRow >= 0 && aRow < 4 && aColumn >= 0 && aColumn < 4);
	flagBits = General;
	return m[aColumn][aRow];
}

inline Vector4D Matrix4D::column(int index) const
{
	Q_ASSERT(index >= 0 && index < 4);
	return Vector4D(m[index][0], m[index][1], m[index][2], m[index][3]);
}

inline void Matrix4D::setColumn(int index, const Vector4D& value)
{
	Q_ASSERT(index >= 0 && index < 4);
	m[index][0] = value.x;
	m[index][1] = value.y;
	m[index][2] = value.z;
	m[index][3] = value.w;
	flagBits = General;
}

inline Vector4D Matrix4D::row(int index) const
{
	Q_ASSERT(index >= 0 && index < 4);
	return Vector4D(m[0][index], m[1][index], m[2][index], m[3][index]);
}

inline void Matrix4D::setRow(int index, const Vector4D& value)
{
	Q_ASSERT(index >= 0 && index < 4);
	m[0][index] = value.x;
	m[1][index] = value.y;
	m[2][index] = value.z;
	m[3][index] = value.w;
	flagBits = General;
}

inline bool Matrix4D::isIdentity() const
{
	if(flagBits == Identity)
		return true;
	if(m[0][0] != 1.0f || m[0][1] != 0.0f || m[0][2] != 0.0f)
		return false;
	if(m[0][3] != 0.0f || m[1][0] != 0.0f || m[1][1] != 1.0f)
		return false;
	if(m[1][2] != 0.0f || m[1][3] != 0.0f || m[2][0] != 0.0f)
		return false;
	if(m[2][1] != 0.0f || m[2][2] != 1.0f || m[2][3] != 0.0f)
		return false;
	if(m[3][0] != 0.0f || m[3][1] != 0.0f || m[3][2] != 0.0f)
		return false;
	return (m[3][3] == 1.0f);
}

inline void Matrix4D::setToIdentity()
{
	m[0][0] = 1.0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[0][3] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = 1.0;
	m[1][2] = 0.0;
	m[1][3] = 0.0;
	m[2][0] = 0.0;
	m[2][1] = 0.0;
	m[2][2] = 1.0;
	m[2][3] = 0.0;
	m[3][0] = 0.0;
	m[3][1] = 0.0;
	m[3][2] = 0.0;
	m[3][3] = 1.0;
	flagBits = Identity;
}

inline void Matrix4D::fill(double value)
{
	m[0][0] = value;
	m[0][1] = value;
	m[0][2] = value;
	m[0][3] = value;
	m[1][0] = value;
	m[1][1] = value;
	m[1][2] = value;
	m[1][3] = value;
	m[2][0] = value;
	m[2][1] = value;
	m[2][2] = value;
	m[2][3] = value;
	m[3][0] = value;
	m[3][1] = value;
	m[3][2] = value;
	m[3][3] = value;
	flagBits = General;
}

inline Matrix4D& Matrix4D::operator+=(const Matrix4D& other)
{
	m[0][0] += other.m[0][0];
	m[0][1] += other.m[0][1];
	m[0][2] += other.m[0][2];
	m[0][3] += other.m[0][3];
	m[1][0] += other.m[1][0];
	m[1][1] += other.m[1][1];
	m[1][2] += other.m[1][2];
	m[1][3] += other.m[1][3];
	m[2][0] += other.m[2][0];
	m[2][1] += other.m[2][1];
	m[2][2] += other.m[2][2];
	m[2][3] += other.m[2][3];
	m[3][0] += other.m[3][0];
	m[3][1] += other.m[3][1];
	m[3][2] += other.m[3][2];
	m[3][3] += other.m[3][3];
	flagBits = General;
	return *this;
}

inline Matrix4D& Matrix4D::operator-=(const Matrix4D& other)
{
	m[0][0] -= other.m[0][0];
	m[0][1] -= other.m[0][1];
	m[0][2] -= other.m[0][2];
	m[0][3] -= other.m[0][3];
	m[1][0] -= other.m[1][0];
	m[1][1] -= other.m[1][1];
	m[1][2] -= other.m[1][2];
	m[1][3] -= other.m[1][3];
	m[2][0] -= other.m[2][0];
	m[2][1] -= other.m[2][1];
	m[2][2] -= other.m[2][2];
	m[2][3] -= other.m[2][3];
	m[3][0] -= other.m[3][0];
	m[3][1] -= other.m[3][1];
	m[3][2] -= other.m[3][2];
	m[3][3] -= other.m[3][3];
	flagBits = General;
	return *this;
}

inline Matrix4D& Matrix4D::operator*=(const Matrix4D& other)
{
	if(flagBits == Identity)
	{
		*this = other;
		return *this;
	}
	else if (other.flagBits == Identity)
	{
		return *this;
	}
	else
	{
		*this = *this * other;
		return *this;
	}
}

inline Matrix4D& Matrix4D::operator*=(double factor)
{
	m[0][0] *= factor;
	m[0][1] *= factor;
	m[0][2] *= factor;
	m[0][3] *= factor;
	m[1][0] *= factor;
	m[1][1] *= factor;
	m[1][2] *= factor;
	m[1][3] *= factor;
	m[2][0] *= factor;
	m[2][1] *= factor;
	m[2][2] *= factor;
	m[2][3] *= factor;
	m[3][0] *= factor;
	m[3][1] *= factor;
	m[3][2] *= factor;
	m[3][3] *= factor;
	flagBits = General;
	return *this;
}

inline bool Matrix4D::operator==(const Matrix4D& other) const
{
	return m[0][0] == other.m[0][0] &&
	       m[0][1] == other.m[0][1] &&
	       m[0][2] == other.m[0][2] &&
	       m[0][3] == other.m[0][3] &&
	       m[1][0] == other.m[1][0] &&
	       m[1][1] == other.m[1][1] &&
	       m[1][2] == other.m[1][2] &&
	       m[1][3] == other.m[1][3] &&
	       m[2][0] == other.m[2][0] &&
	       m[2][1] == other.m[2][1] &&
	       m[2][2] == other.m[2][2] &&
	       m[2][3] == other.m[2][3] &&
	       m[3][0] == other.m[3][0] &&
	       m[3][1] == other.m[3][1] &&
	       m[3][2] == other.m[3][2] &&
	       m[3][3] == other.m[3][3];
}

inline bool Matrix4D::operator!=(const Matrix4D& other) const
{
	return m[0][0] != other.m[0][0] ||
	       m[0][1] != other.m[0][1] ||
	       m[0][2] != other.m[0][2] ||
	       m[0][3] != other.m[0][3] ||
	       m[1][0] != other.m[1][0] ||
	       m[1][1] != other.m[1][1] ||
	       m[1][2] != other.m[1][2] ||
	       m[1][3] != other.m[1][3] ||
	       m[2][0] != other.m[2][0] ||
	       m[2][1] != other.m[2][1] ||
	       m[2][2] != other.m[2][2] ||
	       m[2][3] != other.m[2][3] ||
	       m[3][0] != other.m[3][0] ||
	       m[3][1] != other.m[3][1] ||
	       m[3][2] != other.m[3][2] ||
	       m[3][3] != other.m[3][3];
}

inline Matrix4D operator+(const Matrix4D& m1, const Matrix4D& m2)
{
	Matrix4D m(1);
	m.m[0][0] = m1.m[0][0] + m2.m[0][0];
	m.m[0][1] = m1.m[0][1] + m2.m[0][1];
	m.m[0][2] = m1.m[0][2] + m2.m[0][2];
	m.m[0][3] = m1.m[0][3] + m2.m[0][3];
	m.m[1][0] = m1.m[1][0] + m2.m[1][0];
	m.m[1][1] = m1.m[1][1] + m2.m[1][1];
	m.m[1][2] = m1.m[1][2] + m2.m[1][2];
	m.m[1][3] = m1.m[1][3] + m2.m[1][3];
	m.m[2][0] = m1.m[2][0] + m2.m[2][0];
	m.m[2][1] = m1.m[2][1] + m2.m[2][1];
	m.m[2][2] = m1.m[2][2] + m2.m[2][2];
	m.m[2][3] = m1.m[2][3] + m2.m[2][3];
	m.m[3][0] = m1.m[3][0] + m2.m[3][0];
	m.m[3][1] = m1.m[3][1] + m2.m[3][1];
	m.m[3][2] = m1.m[3][2] + m2.m[3][2];
	m.m[3][3] = m1.m[3][3] + m2.m[3][3];
	return m;
}

inline Matrix4D operator-(const Matrix4D& m1, const Matrix4D& m2)
{
	Matrix4D m(1);
	m.m[0][0] = m1.m[0][0] - m2.m[0][0];
	m.m[0][1] = m1.m[0][1] - m2.m[0][1];
	m.m[0][2] = m1.m[0][2] - m2.m[0][2];
	m.m[0][3] = m1.m[0][3] - m2.m[0][3];
	m.m[1][0] = m1.m[1][0] - m2.m[1][0];
	m.m[1][1] = m1.m[1][1] - m2.m[1][1];
	m.m[1][2] = m1.m[1][2] - m2.m[1][2];
	m.m[1][3] = m1.m[1][3] - m2.m[1][3];
	m.m[2][0] = m1.m[2][0] - m2.m[2][0];
	m.m[2][1] = m1.m[2][1] - m2.m[2][1];
	m.m[2][2] = m1.m[2][2] - m2.m[2][2];
	m.m[2][3] = m1.m[2][3] - m2.m[2][3];
	m.m[3][0] = m1.m[3][0] - m2.m[3][0];
	m.m[3][1] = m1.m[3][1] - m2.m[3][1];
	m.m[3][2] = m1.m[3][2] - m2.m[3][2];
	m.m[3][3] = m1.m[3][3] - m2.m[3][3];
	return m;
}

inline Matrix4D operator*(const Matrix4D& m1, const Matrix4D& m2)
{
	if (m1.flagBits == Matrix4D::Identity)
	    return m2;
	else if (m2.flagBits == Matrix4D::Identity)
	    return m1;
	
	Matrix4D m(1);
	m.m[0][0] = m1.m[0][0] * m2.m[0][0] +
	            m1.m[1][0] * m2.m[0][1] +
	            m1.m[2][0] * m2.m[0][2] +
	            m1.m[3][0] * m2.m[0][3];
	m.m[0][1] = m1.m[0][1] * m2.m[0][0] +
	            m1.m[1][1] * m2.m[0][1] +
	            m1.m[2][1] * m2.m[0][2] +
	            m1.m[3][1] * m2.m[0][3];
	m.m[0][2] = m1.m[0][2] * m2.m[0][0] +
	            m1.m[1][2] * m2.m[0][1] +
	            m1.m[2][2] * m2.m[0][2] +
	            m1.m[3][2] * m2.m[0][3];
	m.m[0][3] = m1.m[0][3] * m2.m[0][0] +
	            m1.m[1][3] * m2.m[0][1] +
	            m1.m[2][3] * m2.m[0][2] +
	            m1.m[3][3] * m2.m[0][3];
	m.m[1][0] = m1.m[0][0] * m2.m[1][0] +
	            m1.m[1][0] * m2.m[1][1] +
	            m1.m[2][0] * m2.m[1][2] +
	            m1.m[3][0] * m2.m[1][3];
	m.m[1][1] = m1.m[0][1] * m2.m[1][0] +
	            m1.m[1][1] * m2.m[1][1] +
	            m1.m[2][1] * m2.m[1][2] +
	            m1.m[3][1] * m2.m[1][3];
	m.m[1][2] = m1.m[0][2] * m2.m[1][0] +
	            m1.m[1][2] * m2.m[1][1] +
	            m1.m[2][2] * m2.m[1][2] +
	            m1.m[3][2] * m2.m[1][3];
	m.m[1][3] = m1.m[0][3] * m2.m[1][0] +
	            m1.m[1][3] * m2.m[1][1] +
	            m1.m[2][3] * m2.m[1][2] +
	            m1.m[3][3] * m2.m[1][3];
	m.m[2][0] = m1.m[0][0] * m2.m[2][0] +
	            m1.m[1][0] * m2.m[2][1] +
	            m1.m[2][0] * m2.m[2][2] +
	            m1.m[3][0] * m2.m[2][3];
	m.m[2][1] = m1.m[0][1] * m2.m[2][0] +
	            m1.m[1][1] * m2.m[2][1] +
	            m1.m[2][1] * m2.m[2][2] +
	            m1.m[3][1] * m2.m[2][3];
	m.m[2][2] = m1.m[0][2] * m2.m[2][0] +
	            m1.m[1][2] * m2.m[2][1] +
	            m1.m[2][2] * m2.m[2][2] +
	            m1.m[3][2] * m2.m[2][3];
	m.m[2][3] = m1.m[0][3] * m2.m[2][0] +
	            m1.m[1][3] * m2.m[2][1] +
	            m1.m[2][3] * m2.m[2][2] +
	            m1.m[3][3] * m2.m[2][3];
	m.m[3][0] = m1.m[0][0] * m2.m[3][0] +
	            m1.m[1][0] * m2.m[3][1] +
	            m1.m[2][0] * m2.m[3][2] +
	            m1.m[3][0] * m2.m[3][3];
	m.m[3][1] = m1.m[0][1] * m2.m[3][0] +
	            m1.m[1][1] * m2.m[3][1] +
	            m1.m[2][1] * m2.m[3][2] +
	            m1.m[3][1] * m2.m[3][3];
	m.m[3][2] = m1.m[0][2] * m2.m[3][0] +
	            m1.m[1][2] * m2.m[3][1] +
	            m1.m[2][2] * m2.m[3][2] +
	            m1.m[3][2] * m2.m[3][3];
	m.m[3][3] = m1.m[0][3] * m2.m[3][0] +
	            m1.m[1][3] * m2.m[3][1] +
	            m1.m[2][3] * m2.m[3][2] +
	            m1.m[3][3] * m2.m[3][3];
	return m;
}

inline Vector3D operator*(const Vector3D& vector, const Matrix4D& matrix)
{
	double x, y, z, w;
	x = vector.x * matrix.m[0][0] +
	    vector.y * matrix.m[0][1] +
	    vector.z * matrix.m[0][2] +
	    matrix.m[0][3];
	y = vector.x * matrix.m[1][0] +
	    vector.y * matrix.m[1][1] +
	    vector.z * matrix.m[1][2] +
	    matrix.m[1][3];
	z = vector.x * matrix.m[2][0] +
	    vector.y * matrix.m[2][1] +
	    vector.z * matrix.m[2][2] +
	    matrix.m[2][3];
	w = vector.x * matrix.m[3][0] +
	    vector.y * matrix.m[3][1] +
	    vector.z * matrix.m[3][2] +
	    matrix.m[3][3];
	if (w == 1.0f)
	    return Vector3D(x, y, z);
	else
	    return Vector3D(x / w, y / w, z / w);
}

inline Vector3D operator*(const Matrix4D& matrix, const Vector3D& vector)
{
	double x, y, z, w;
	if (matrix.flagBits == Matrix4D::Identity) 
	{
	    return vector;
	} 
	else if (matrix.flagBits == Matrix4D::Translation) 
	{
	    return  Vector3D(vector.x + matrix.m[3][0],
	                     vector.y + matrix.m[3][1],
	                     vector.z + matrix.m[3][2]);
	} else if (matrix.flagBits ==
	                (Matrix4D::Translation | Matrix4D::Scale)) {
	    return Vector3D(vector.x * matrix.m[0][0] + matrix.m[3][0],
	                    vector.y * matrix.m[1][1] + matrix.m[3][1],
	                    vector.z * matrix.m[2][2] + matrix.m[3][2]);
	} 
	else if (matrix.flagBits == Matrix4D::Scale) 
	{
	    return Vector3D(vector.x * matrix.m[0][0],
	                    vector.y * matrix.m[1][1],
	                    vector.z * matrix.m[2][2]);
	} 
	else 
	{
	    x = vector.x * matrix.m[0][0] +
	        vector.y * matrix.m[1][0] +
	        vector.z * matrix.m[2][0] +
	        matrix.m[3][0];
	    y = vector.x * matrix.m[0][1] +
	        vector.y * matrix.m[1][1] +
	        vector.z * matrix.m[2][1] +
	        matrix.m[3][1];
	    z = vector.x * matrix.m[0][2] +
	        vector.y * matrix.m[1][2] +
	        vector.z * matrix.m[2][2] +
	        matrix.m[3][2];
	    w = vector.x * matrix.m[0][3] +
	        vector.y * matrix.m[1][3] +
	        vector.z * matrix.m[2][3] +
	        matrix.m[3][3];
	    if (w == 1.0f)
	        return Vector3D(x, y, z);
	    else
	        return Vector3D(x / w, y / w, z / w);
	}
}

inline Vector4D operator*(const Vector4D& vector, const Matrix4D& matrix)
{
	double x, y, z, w;
	x = vector.x * matrix.m[0][0] +
	    vector.y * matrix.m[0][1] +
	    vector.z * matrix.m[0][2] +
	    vector.w * matrix.m[0][3];
	y = vector.x * matrix.m[1][0] +
	    vector.y * matrix.m[1][1] +
	    vector.z * matrix.m[1][2] +
	    vector.w * matrix.m[1][3];
	z = vector.x * matrix.m[2][0] +
	    vector.y * matrix.m[2][1] +
	    vector.z * matrix.m[2][2] +
	    vector.w * matrix.m[2][3];
	w = vector.x * matrix.m[3][0] +
	    vector.y * matrix.m[3][1] +
	    vector.z * matrix.m[3][2] +
	    vector.w * matrix.m[3][3];
	return Vector4D(x, y, z, w);
}

inline Vector4D operator*(const Matrix4D& matrix, const Vector4D& vector)
{
	double x, y, z, w;
	x = vector.x * matrix.m[0][0] +
	    vector.y * matrix.m[1][0] +
	    vector.z * matrix.m[2][0] +
	    vector.w * matrix.m[3][0];
	y = vector.x * matrix.m[0][1] +
	    vector.y * matrix.m[1][1] +
	    vector.z * matrix.m[2][1] +
	    vector.w * matrix.m[3][1];
	z = vector.x * matrix.m[0][2] +
	    vector.y * matrix.m[1][2] +
	    vector.z * matrix.m[2][2] +
	    vector.w * matrix.m[3][2];
	w = vector.x * matrix.m[0][3] +
	    vector.y * matrix.m[1][3] +
	    vector.z * matrix.m[2][3] +
	    vector.w * matrix.m[3][3];
	return Vector4D(x, y, z, w);
}

inline QPoint operator*(const QPoint& point, const Matrix4D& matrix)
{
	double xin, yin;
	double x, y, w;
	xin = point.x();
	yin = point.y();
	x = xin * matrix.m[0][0] +
	    yin * matrix.m[0][1] +
	    matrix.m[0][3];
	y = xin * matrix.m[1][0] +
	    yin * matrix.m[1][1] +
	    matrix.m[1][3];
	w = xin * matrix.m[3][0] +
	    yin * matrix.m[3][1] +
	    matrix.m[3][3];
	if (w == 1.0f)
	    return QPoint(qRound(x), qRound(y));
	else
	    return QPoint(qRound(x / w), qRound(y / w));
}

inline QPointF operator*(const QPointF& point, const Matrix4D& matrix)
{
	double xin, yin;
	double x, y, w;
	xin = point.x();
	yin = point.y();
	x = xin * matrix.m[0][0] +
	    yin * matrix.m[0][1] +
	    matrix.m[0][3];
	y = xin * matrix.m[1][0] +
	    yin * matrix.m[1][1] +
	    matrix.m[1][3];
	w = xin * matrix.m[3][0] +
	    yin * matrix.m[3][1] +
	    matrix.m[3][3];
	if (w == 1.0f) {
	    return QPointF(double(x), double(y));
	} else {
	    return QPointF(double(x / w), double(y / w));
	}
}

inline QPoint operator*(const Matrix4D& matrix, const QPoint& point)
{
	double xin, yin;
	double x, y, w;
	xin = point.x();
	yin = point.y();
	if (matrix.flagBits == Matrix4D::Identity) {
	    return point;
	} else if (matrix.flagBits == Matrix4D::Translation) {
	    return QPoint(qRound(xin + matrix.m[3][0]),
	                  qRound(yin + matrix.m[3][1]));
	} else if (matrix.flagBits ==
	                (Matrix4D::Translation | Matrix4D::Scale)) {
	    return QPoint(qRound(xin * matrix.m[0][0] + matrix.m[3][0]),
	                  qRound(yin * matrix.m[1][1] + matrix.m[3][1]));
	} else if (matrix.flagBits == Matrix4D::Scale) {
	    return QPoint(qRound(xin * matrix.m[0][0]),
	                  qRound(yin * matrix.m[1][1]));
	} else {
	    x = xin * matrix.m[0][0] +
	        yin * matrix.m[1][0] +
	        matrix.m[3][0];
	    y = xin * matrix.m[0][1] +
	        yin * matrix.m[1][1] +
	        matrix.m[3][1];
	    w = xin * matrix.m[0][3] +
	        yin * matrix.m[1][3] +
	        matrix.m[3][3];
	    if (w == 1.0f)
	        return QPoint(qRound(x), qRound(y));
	    else
	        return QPoint(qRound(x / w), qRound(y / w));
	}
}

inline QPointF operator*(const Matrix4D& matrix, const QPointF& point)
{
	double xin, yin;
	double x, y, w;
	xin = point.x();
	yin = point.y();
	if (matrix.flagBits == Matrix4D::Identity) {
	    return point;
	} else if (matrix.flagBits == Matrix4D::Translation) {
	    return QPointF(xin + matrix.m[3][0],
	                   yin + matrix.m[3][1]);
	} else if (matrix.flagBits ==
	                (Matrix4D::Translation | Matrix4D::Scale)) {
	    return QPointF(xin * matrix.m[0][0] + matrix.m[3][0],
	                   yin * matrix.m[1][1] + matrix.m[3][1]);
	} else if (matrix.flagBits == Matrix4D::Scale) {
	    return QPointF(xin * matrix.m[0][0],
	                   yin * matrix.m[1][1]);
	} else {
	    x = xin * matrix.m[0][0] +
	        yin * matrix.m[1][0] +
	        matrix.m[3][0];
	    y = xin * matrix.m[0][1] +
	        yin * matrix.m[1][1] +
	        matrix.m[3][1];
	    w = xin * matrix.m[0][3] +
	        yin * matrix.m[1][3] +
	        matrix.m[3][3];
	    if (w == 1.0f) {
	        return QPointF(double(x), double(y));
	    } else {
	        return QPointF(double(x / w), double(y / w));
	    }
	}
}

inline Matrix4D operator-(const Matrix4D& matrix)
{
	Matrix4D m(1);
	m.m[0][0] = -matrix.m[0][0];
	m.m[0][1] = -matrix.m[0][1];
	m.m[0][2] = -matrix.m[0][2];
	m.m[0][3] = -matrix.m[0][3];
	m.m[1][0] = -matrix.m[1][0];
	m.m[1][1] = -matrix.m[1][1];
	m.m[1][2] = -matrix.m[1][2];
	m.m[1][3] = -matrix.m[1][3];
	m.m[2][0] = -matrix.m[2][0];
	m.m[2][1] = -matrix.m[2][1];
	m.m[2][2] = -matrix.m[2][2];
	m.m[2][3] = -matrix.m[2][3];
	m.m[3][0] = -matrix.m[3][0];
	m.m[3][1] = -matrix.m[3][1];
	m.m[3][2] = -matrix.m[3][2];
	m.m[3][3] = -matrix.m[3][3];
	return m;
}

inline Matrix4D operator*(double factor, const Matrix4D& matrix)
{
	Matrix4D m(1);
	m.m[0][0] = matrix.m[0][0] * factor;
	m.m[0][1] = matrix.m[0][1] * factor;
	m.m[0][2] = matrix.m[0][2] * factor;
	m.m[0][3] = matrix.m[0][3] * factor;
	m.m[1][0] = matrix.m[1][0] * factor;
	m.m[1][1] = matrix.m[1][1] * factor;
	m.m[1][2] = matrix.m[1][2] * factor;
	m.m[1][3] = matrix.m[1][3] * factor;
	m.m[2][0] = matrix.m[2][0] * factor;
	m.m[2][1] = matrix.m[2][1] * factor;
	m.m[2][2] = matrix.m[2][2] * factor;
	m.m[2][3] = matrix.m[2][3] * factor;
	m.m[3][0] = matrix.m[3][0] * factor;
	m.m[3][1] = matrix.m[3][1] * factor;
	m.m[3][2] = matrix.m[3][2] * factor;
	m.m[3][3] = matrix.m[3][3] * factor;
	return m;
}

inline Matrix4D operator*(const Matrix4D& matrix, double factor)
{
	Matrix4D m(1);
	m.m[0][0] = matrix.m[0][0] * factor;
	m.m[0][1] = matrix.m[0][1] * factor;
	m.m[0][2] = matrix.m[0][2] * factor;
	m.m[0][3] = matrix.m[0][3] * factor;
	m.m[1][0] = matrix.m[1][0] * factor;
	m.m[1][1] = matrix.m[1][1] * factor;
	m.m[1][2] = matrix.m[1][2] * factor;
	m.m[1][3] = matrix.m[1][3] * factor;
	m.m[2][0] = matrix.m[2][0] * factor;
	m.m[2][1] = matrix.m[2][1] * factor;
	m.m[2][2] = matrix.m[2][2] * factor;
	m.m[2][3] = matrix.m[2][3] * factor;
	m.m[3][0] = matrix.m[3][0] * factor;
	m.m[3][1] = matrix.m[3][1] * factor;
	m.m[3][2] = matrix.m[3][2] * factor;
	m.m[3][3] = matrix.m[3][3] * factor;
	return m;
}

inline bool qFuzzyCompare(const Matrix4D& m1, const Matrix4D& m2)
{
	return qFuzzyCompare(m1.m[0][0], m2.m[0][0]) &&
	       qFuzzyCompare(m1.m[0][1], m2.m[0][1]) &&
	       qFuzzyCompare(m1.m[0][2], m2.m[0][2]) &&
	       qFuzzyCompare(m1.m[0][3], m2.m[0][3]) &&
	       qFuzzyCompare(m1.m[1][0], m2.m[1][0]) &&
	       qFuzzyCompare(m1.m[1][1], m2.m[1][1]) &&
	       qFuzzyCompare(m1.m[1][2], m2.m[1][2]) &&
	       qFuzzyCompare(m1.m[1][3], m2.m[1][3]) &&
	       qFuzzyCompare(m1.m[2][0], m2.m[2][0]) &&
	       qFuzzyCompare(m1.m[2][1], m2.m[2][1]) &&
	       qFuzzyCompare(m1.m[2][2], m2.m[2][2]) &&
	       qFuzzyCompare(m1.m[2][3], m2.m[2][3]) &&
	       qFuzzyCompare(m1.m[3][0], m2.m[3][0]) &&
	       qFuzzyCompare(m1.m[3][1], m2.m[3][1]) &&
	       qFuzzyCompare(m1.m[3][2], m2.m[3][2]) &&
	       qFuzzyCompare(m1.m[3][3], m2.m[3][3]);
}

inline Matrix4D& Matrix4D::operator/=(float divisor)
{
	if (divisor == 1.0f)
		return *this;

  m[0][0] /= divisor;
  m[0][1] /= divisor;
  m[0][2] /= divisor;
  m[0][3] /= divisor;
  m[1][0] /= divisor;
  m[1][1] /= divisor;
  m[1][2] /= divisor;
  m[1][3] /= divisor;
  m[2][0] /= divisor;
  m[2][1] /= divisor;
  m[2][2] /= divisor;
  m[2][3] /= divisor;
  m[3][0] /= divisor;
  m[3][1] /= divisor;
  m[3][2] /= divisor;
  m[3][3] /= divisor;
  flagBits = General;
  return *this;
}

inline Matrix4D operator/(const Matrix4D& matrix, float divisor)
{
	if (divisor == 1.0f)
		return matrix;
  Matrix4D m(1); // The "1" says to not load the identity.
  m.m[0][0] = matrix.m[0][0] / divisor;
  m.m[0][1] = matrix.m[0][1] / divisor;
  m.m[0][2] = matrix.m[0][2] / divisor;
  m.m[0][3] = matrix.m[0][3] / divisor;
  m.m[1][0] = matrix.m[1][0] / divisor;
  m.m[1][1] = matrix.m[1][1] / divisor;
  m.m[1][2] = matrix.m[1][2] / divisor;
  m.m[1][3] = matrix.m[1][3] / divisor;
  m.m[2][0] = matrix.m[2][0] / divisor;
  m.m[2][1] = matrix.m[2][1] / divisor;
  m.m[2][2] = matrix.m[2][2] / divisor;
  m.m[2][3] = matrix.m[2][3] / divisor;
  m.m[3][0] = matrix.m[3][0] / divisor;
  m.m[3][1] = matrix.m[3][1] / divisor;
  m.m[3][2] = matrix.m[3][2] / divisor;
  m.m[3][3] = matrix.m[3][3] / divisor;
  return m;
}

/*!
    \fn bool qFuzzyCompare(const Matrix4D& m1, const Matrix4D& m2)
    \relates Matrix4D

    Returns true if \a m1 and \a m2 are equal, allowing for a small
    fuzziness factor for floating-point comparisons; false otherwise.
*/

inline void Matrix4D::scale(double x, double y)
{
	if (x == 1.0 && y == 1.0)
		return;
	if (flagBits == Identity) {
	    m[0][0] = x;
	    m[1][1] = y;
	    flagBits = Scale;
	} else if (flagBits == Scale || flagBits == (Scale | Translation)) {
	    m[0][0] *= x;
	    m[1][1] *= y;
	} else if (flagBits == Translation) {
	    m[0][0] = x;
	    m[1][1] = y;
	    flagBits |= Scale;
	} else {
	    m[0][0] *= x;
	    m[0][1] *= x;
	    m[0][2] *= x;
	    m[0][3] *= x;
	    m[1][0] *= y;
	    m[1][1] *= y;
	    m[1][2] *= y;
	    m[1][3] *= y;
	    flagBits = General;
	}
}

inline void Matrix4D::scale(double x, double y, double z)
{
	if (x == 1.0 && y == 1.0 && z == 1.0)
		return;
	if (flagBits == Identity) {
	    m[0][0] = x;
	    m[1][1] = y;
	    m[2][2] = z;
	    flagBits = Scale;
	} else if (flagBits == Scale || flagBits == (Scale | Translation)) {
	    m[0][0] *= x;
	    m[1][1] *= y;
	    m[2][2] *= z;
	} else if (flagBits == Translation) {
	    m[0][0] = x;
	    m[1][1] = y;
	    m[2][2] = z;
	    flagBits |= Scale;
	} else {
	    m[0][0] *= x;
	    m[0][1] *= x;
	    m[0][2] *= x;
	    m[0][3] *= x;
	    m[1][0] *= y;
	    m[1][1] *= y;
	    m[1][2] *= y;
	    m[1][3] *= y;
	    m[2][0] *= z;
	    m[2][1] *= z;
	    m[2][2] *= z;
	    m[2][3] *= z;
	    flagBits = General;
	}
}

inline void Matrix4D::scale(float factor)
{
	if (factor == 1.0f)
		return;

	if (flagBits == Identity) {
	    m[0][0] = factor;
	    m[1][1] = factor;
	    m[2][2] = factor;
	    flagBits = Scale;
	} else if (flagBits == Scale || flagBits == (Scale | Translation)) {
	    m[0][0] *= factor;
	    m[1][1] *= factor;
	    m[2][2] *= factor;
	} else if (flagBits == Translation) {
	    m[0][0] = factor;
	    m[1][1] = factor;
	    m[2][2] = factor;
	    flagBits |= Scale;
	} else {
	    m[0][0] *= factor;
	    m[0][1] *= factor;
	    m[0][2] *= factor;
	    m[0][3] *= factor;
	    m[1][0] *= factor;
	    m[1][1] *= factor;
	    m[1][2] *= factor;
	    m[1][3] *= factor;
	    m[2][0] *= factor;
	    m[2][1] *= factor;
	    m[2][2] *= factor;
	    m[2][3] *= factor;
	    flagBits = General;
	}
}

inline void Matrix4D::translate(double x, double y)
{
	if (x == 0.0 && y == 0.0)
		return;
	if (flagBits == Identity) {
	    m[3][0] = x;
	    m[3][1] = y;
	    flagBits = Translation;
	} else if (flagBits == Translation) {
	    m[3][0] += x;
	    m[3][1] += y;
	} else if (flagBits == Scale) {
	    m[3][0] = m[0][0] * x;
	    m[3][1] = m[1][1] * y;
	    m[3][2] = 0.;
	    flagBits |= Translation;
	} else if (flagBits == (Scale | Translation)) {
	    m[3][0] += m[0][0] * x;
	    m[3][1] += m[1][1] * y;
	} else {
	    m[3][0] += m[0][0] * x + m[1][0] * y;
	    m[3][1] += m[0][1] * x + m[1][1] * y;
	    m[3][2] += m[0][2] * x + m[1][2] * y;
	    m[3][3] += m[0][3] * x + m[1][3] * y;
	    if (flagBits == Rotation)
	        flagBits |= Translation;
	    else if (flagBits != (Rotation | Translation))
	        flagBits = General;
	}
}

inline void Matrix4D::translate(double x, double y, double z)
{
	if (x == 0.0 && y == 0.0 && z == 0.0)
		return;
	if (flagBits == Identity) {
	    m[3][0] = x;
	    m[3][1] = y;
	    m[3][2] = z;
	    flagBits = Translation;
	} else if (flagBits == Translation) {
	    m[3][0] += x;
	    m[3][1] += y;
	    m[3][2] += z;
	} else if (flagBits == Scale) {
	    m[3][0] = m[0][0] * x;
	    m[3][1] = m[1][1] * y;
	    m[3][2] = m[2][2] * z;
	    flagBits |= Translation;
	} else if (flagBits == (Scale | Translation)) {
	    m[3][0] += m[0][0] * x;
	    m[3][1] += m[1][1] * y;
	    m[3][2] += m[2][2] * z;
	} else {
	    m[3][0] += m[0][0] * x + m[1][0] * y + m[2][0] * z;
	    m[3][1] += m[0][1] * x + m[1][1] * y + m[2][1] * z;
	    m[3][2] += m[0][2] * x + m[1][2] * y + m[2][2] * z;
	    m[3][3] += m[0][3] * x + m[1][3] * y + m[2][3] * z;
	    if (flagBits == Rotation)
	        flagBits |= Translation;
	    else if (flagBits != (Rotation | Translation))
	        flagBits = General;
	}
}

/*!
    Multiples this matrix by another that rotates coordinates according
    to a specified \a quaternion.  The \a quaternion is assumed to have
    been normalized.

    \sa scale(), translate(), QQuaternion
*/
inline void Matrix4D::rotate(const QuaternionD& quaternion)
{
	// Algorithm from:
	// http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q54
	Matrix4D m(1);
	double xx = quaternion.x() * quaternion.x();
	double xy = quaternion.x() * quaternion.y();
	double xz = quaternion.x() * quaternion.z();
	double xw = quaternion.x() * quaternion.scalar();
	double yy = quaternion.y() * quaternion.y();
	double yz = quaternion.y() * quaternion.z();
	double yw = quaternion.y() * quaternion.scalar();
	double zz = quaternion.z() * quaternion.z();
	double zw = quaternion.z() * quaternion.scalar();
	m.m[0][0] = 1.0f - 2 * (yy + zz);
	m.m[1][0] =        2 * (xy - zw);
	m.m[2][0] =        2 * (xz + yw);
	m.m[3][0] = 0.0f;
	m.m[0][1] =        2 * (xy + zw);
	m.m[1][1] = 1.0f - 2 * (xx + zz);
	m.m[2][1] =        2 * (yz - xw);
	m.m[3][1] = 0.0f;
	m.m[0][2] =        2 * (xz - yw);
	m.m[1][2] =        2 * (yz + xw);
	m.m[2][2] = 1.0f - 2 * (xx + yy);
	m.m[3][2] = 0.0f;
	m.m[0][3] = 0.0f;
	m.m[1][3] = 0.0f;
	m.m[2][3] = 0.0f;
	m.m[3][3] = 1.0f;
	int flags = flagBits;
	*this *= m;
	if (flags != Identity)
	    flagBits = flags | Rotation;
	else
	    flagBits = Rotation;
}

//
inline Matrix4D Matrix4D::ortho(double left, double right, double bottom, double top, double nearPlane, double farPlane)
{
	// Bail out if the projection volume is zero-sized.
	if (left == right || bottom == top || nearPlane == farPlane)
	    return Matrix4D();
	
	// Construct the projection.
	double width = right - left;
	double invheight = top - bottom;
	double clip = farPlane - nearPlane;
	if (clip == 2.0f && (nearPlane + farPlane) == 0.0f) 
	{
	    // We can express this projection as a translate and scale
	    // which will be more efficient to modify with further
	    // transformations than producing a "General" matrix
			Matrix4D m;
	    m.translate(-(left + right) / width, -(top + bottom) / invheight, 0.0f);
	    m.scale(2.0f / width, 2.0f / invheight, -1.0f);
	    return m;
	}
	Matrix4D m(1);
	m.m[0][0] = 2.0f / width;
	m.m[1][0] = 0.0f;
	m.m[2][0] = 0.0f;
	m.m[3][0] = -(left + right) / width;
	m.m[0][1] = 0.0f;
	m.m[1][1] = 2.0f / invheight;
	m.m[2][1] = 0.0f;
	m.m[3][1] = -(top + bottom) / invheight;
	m.m[0][2] = 0.0f;
	m.m[1][2] = 0.0f;
	m.m[2][2] = -2.0f / clip;
	m.m[3][2] = -(nearPlane + farPlane) / clip;
	m.m[0][3] = 0.0f;
	m.m[1][3] = 0.0f;
	m.m[2][3] = 0.0f;
	m.m[3][3] = 1.0f;
	
	return m;
}

inline Matrix4D Matrix4D::frustum(double left, double right, double bottom, double top, double nearPlane, double farPlane)
{
	// Bail out if the projection volume is zero-sized.
	if (left == right || bottom == top || nearPlane == farPlane)
	    return Matrix4D();
	
	// Construct the projection.
	Matrix4D m(1);
	double width = right - left;
	double invheight = top - bottom;
	double clip = farPlane - nearPlane;
	m.m[0][0] = 2.0f * nearPlane / width;
	m.m[1][0] = 0.0f;
	m.m[2][0] = (left + right) / width;
	m.m[3][0] = 0.0f;
	m.m[0][1] = 0.0f;
	m.m[1][1] = 2.0f * nearPlane / invheight;
	m.m[2][1] = (top + bottom) / invheight;
	m.m[3][1] = 0.0f;
	m.m[0][2] = 0.0f;
	m.m[1][2] = 0.0f;
	m.m[2][2] = -(nearPlane + farPlane) / clip;
	m.m[3][2] = -2.0f * nearPlane * farPlane / clip;
	m.m[0][3] = 0.0f;
	m.m[1][3] = 0.0f;
	m.m[2][3] = -1.0f;
	m.m[3][3] = 0.0f;
	
	return m;
}

/*!
    Multiplies this matrix by another that applies a perspective
    projection.  The field of view will be \a angle degrees within
    a window with a given \a aspect ratio.  The projection will
    have the specified \a nearPlane and \a farPlane clipping planes.

    \sa ortho(), frustum()
*/
inline Matrix4D Matrix4D::lookAt(const Vector3D& eye, const Vector3D& center, const Vector3D& up)
{
	Vector3D forward = (center - eye).normalized();
	Vector3D side = Vector3D::crossProduct(forward, up).normalized();
	Vector3D upVector = Vector3D::crossProduct(side, forward);
	
	Matrix4D m(1);
	
	m.m[0][0] = side.x;
	m.m[1][0] = side.y;
	m.m[2][0] = side.z;
	m.m[3][0] = 0.0;
	m.m[0][1] = upVector.x;
	m.m[1][1] = upVector.y;
	m.m[2][1] = upVector.z;
	m.m[3][1] = 0.0;
	m.m[0][2] = -forward.x;
	m.m[1][2] = -forward.y;
	m.m[2][2] = -forward.z;
	m.m[3][2] = 0.0;
	m.m[0][3] = 0.0;
	m.m[1][3] = 0.0;
	m.m[2][3] = 0.0;
	m.m[3][3] = 1.0;
	
	m.translate(-eye);
	return m;
}
//
/*!
    Flips between right-handed and left-handed coordinate systems
    by multiplying the y and z co-ordinates by -1.  This is normally
    used to create a left-handed orthographic view without scaling
    the viewport as ortho() does.

    \sa ortho()
*/
inline void Matrix4D::flipCoordinates()
{
	if (flagBits == Scale || flagBits == (Scale | Translation)) {
	    m[1][1] = -m[1][1];
	    m[2][2] = -m[2][2];
	} else if (flagBits == Translation) {
	    m[1][1] = -m[1][1];
	    m[2][2] = -m[2][2];
	    flagBits |= Scale;
	} else if (flagBits == Identity) {
	    m[1][1] = -1.0f;
	    m[2][2] = -1.0f;
	    flagBits = Scale;
	} else {
	    m[1][0] = -m[1][0];
	    m[1][1] = -m[1][1];
	    m[1][2] = -m[1][2];
	    m[1][3] = -m[1][3];
	    m[2][0] = -m[2][0];
	    m[2][1] = -m[2][1];
	    m[2][2] = -m[2][2];
	    m[2][3] = -m[2][3];
	    flagBits = General;
	}
}

inline Matrix4D Matrix4D::orthonormalInverse() const
{
	Matrix4D result(1);  // The '1' says not to load identity
	
	result.m[0][0] = m[0][0];
	result.m[1][0] = m[0][1];
	result.m[2][0] = m[0][2];
	
	result.m[0][1] = m[1][0];
	result.m[1][1] = m[1][1];
	result.m[2][1] = m[1][2];
	
	result.m[0][2] = m[2][0];
	result.m[1][2] = m[2][1];
	result.m[2][2] = m[2][2];
	
	result.m[0][3] = 0.0f;
	result.m[1][3] = 0.0f;
	result.m[2][3] = 0.0f;
	
	result.m[3][0] = -(result.m[0][0] * m[3][0] + result.m[1][0] * m[3][1] + result.m[2][0] * m[3][2]);
	result.m[3][1] = -(result.m[0][1] * m[3][0] + result.m[1][1] * m[3][1] + result.m[2][1] * m[3][2]);
	result.m[3][2] = -(result.m[0][2] * m[3][0] + result.m[1][2] * m[3][1] + result.m[2][2] * m[3][2]);
	result.m[3][3] = 1.0f;
	
	return result;
}

/*!
    Optimize the usage of this matrix from its current elements.

    Some operations such as translate(), scale(), and rotate() can be
    performed more efficiently if the matrix being modified is already
    known to be the identity, a previous translate(), a previous
    scale(), etc.

    Normally the Matrix4D class keeps track of this special type internally
    as operations are performed.  However, if the matrix is modified
    directly with operator()() or data(), then Matrix4D will lose track of
    the special type and will revert to the safest but least efficient
    operations thereafter.

    By calling optimize() after directly modifying the matrix,
    the programmer can force Matrix4D to recover the special type if
    the elements appear to conform to one of the known optimized types.

    \sa operator()(), data(), translate()
*/
inline void Matrix4D::optimize()
{
	// If the last element is not 1, then it can never be special.
	if (m[3][3] != 1.0f) {
	    flagBits = General;
	    return;
	}
	
	// If the upper three elements m12, m13, and m21 are not all zero,
	// or the lower elements below the diagonal are not all zero, then
	// the matrix can never be special.
	if (m[1][0] != 0.0f || m[2][0] != 0.0f || m[2][1] != 0.0f) {
	    flagBits = General;
	    return;
	}
	if (m[0][1] != 0.0f || m[0][2] != 0.0f || m[0][3] != 0.0f ||
	    m[1][2] != 0.0f || m[1][3] != 0.0f || m[2][3] != 0.0f) {
	    flagBits = General;
	    return;
	}
	
	// Determine what we have in the remaining regions of the matrix.
	bool identityAlongDiagonal
	    = (m[0][0] == 1.0f && m[1][1] == 1.0f && m[2][2] == 1.0f);
	bool translationPresent
	    = (m[3][0] != 0.0f || m[3][1] != 0.0f || m[3][2] != 0.0f);
	
	// Now determine the special matrix type.
	if (translationPresent && identityAlongDiagonal)
	    flagBits = Translation;
	else if (translationPresent)
	    flagBits = (Translation | Scale);
	else if (identityAlongDiagonal)
	    flagBits = Identity;
	else
	    flagBits = Scale;
}

inline QRectF Matrix4D::mapRect(const QRectF& rect) const
{
	if (flagBits == (Translation | Scale) || flagBits == Scale) {
	    double x = rect.x() * m[0][0] + m[3][0];
	    double y = rect.y() * m[1][1] + m[3][1];
	    double w = rect.width() * m[0][0];
	    double h = rect.height() * m[1][1];
	    if (w < 0) {
	        w = -w;
	        x -= w;
	    }
	    if (h < 0) {
	        h = -h;
	        y -= h;
	    }
	    return QRectF(x, y, w, h);
	} else if (flagBits == Translation) {
	    return rect.translated(m[3][0], m[3][1]);
	}
	
	QPointF tl = map(rect.topLeft()); QPointF tr = map(rect.topRight());
	QPointF bl = map(rect.bottomLeft()); QPointF br = map(rect.bottomRight());
	
	double xmin = qMin(qMin(tl.x(), tr.x()), qMin(bl.x(), br.x()));
	double xmax = qMax(qMax(tl.x(), tr.x()), qMax(bl.x(), br.x()));
	double ymin = qMin(qMin(tl.y(), tr.y()), qMin(bl.y(), br.y()));
	double ymax = qMax(qMax(tl.y(), tr.y()), qMax(bl.y(), br.y()));
	
	return QRectF(QPointF(xmin, ymin), QPointF(xmax, ymax));
}

#endif
