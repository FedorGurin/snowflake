#ifndef VECTOR3D_H
#define VECTOR3D_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <QVector3D>
#include <qmath.h>

class Vector3D
{
public:
	union{
		struct
		{
			double x,y,z;
		};
		struct
		{
			double r,g,b;
		};
		double v[3];
	};

	Vector3D():x(0),y(0),z(0) {}
	Vector3D(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
	Vector3D(const Vector3D &vec): x(vec.x), y(vec.y), z(vec.z) {}
	explicit Vector3D(const QVector3D &vec): x(vec.x()), y(vec.y()), z(vec.z()) {}
	explicit Vector3D(const double vec[3]): x(vec[0]), y(vec[1]), z(vec[2]) {}
	QVector3D toQVector3D() const { return QVector3D(x,y,z); }

	inline Vector3D normalized() const;
  inline void normalize();

	double length() const { return sqrt(x*x+y*y+z*z); }
	double squaredLength() const { return x*x+y*y+z*z; }
	double &operator[](int a) { return v[a]; }
	double operator[](int a) const { return v[a]; }

	static double dotProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	static Vector3D crossProduct(const Vector3D& v1, const Vector3D& v2)
	{
		return Vector3D(v1.y * v2.z - v1.z * v2.y,
										v1.z * v2.x - v1.x * v2.z,
										v1.x * v2.y - v1.y * v2.x);
	}

	inline Vector3D &operator+=(const Vector3D &vector);
  inline Vector3D &operator-=(const Vector3D &vector);
  inline Vector3D &operator*=(double factor);
  inline Vector3D &operator*=(const Vector3D& vector);
  inline Vector3D &operator/=(double divisor);

	friend inline bool operator==(const Vector3D &v1, const Vector3D &v2);
  friend inline bool operator!=(const Vector3D &v1, const Vector3D &v2);
  friend inline const Vector3D operator+(const Vector3D &v1, const Vector3D &v2);
  friend inline const Vector3D operator-(const Vector3D &v1, const Vector3D &v2);
  friend inline const Vector3D operator*(double factor, const Vector3D &vector);
  friend inline const Vector3D operator*(const Vector3D &vector, double factor);
  friend const Vector3D operator*(const Vector3D &v1, const Vector3D& v2);
  friend inline const Vector3D operator-(const Vector3D &vector);
  friend inline const Vector3D operator/(const Vector3D &vector, double divisor);

  friend inline bool qFuzzyCompare(const Vector3D& v1, const Vector3D& v2);
};

QDataStream &operator<<(QDataStream &stream, const Vector3D &vector);

QDataStream &operator>>(QDataStream &stream, Vector3D &vector);

inline Vector3D &Vector3D::operator+=(const Vector3D &vector)
{
	x += vector.x;
	y += vector.y;
	z += vector.z;
	return *this;
}

inline Vector3D &Vector3D::operator-=(const Vector3D &vector)
{
	x -= vector.x;
	y -= vector.y;
	z -= vector.z;
	return *this;
}

inline Vector3D &Vector3D::operator*=(double factor)
{
	x *= factor;
	y *= factor;
	z *= factor;
	return *this;
}

inline Vector3D &Vector3D::operator*=(const Vector3D& vector)
{
	x *= vector.x;
	y *= vector.y;
	z *= vector.z;
	return *this;
}

inline Vector3D &Vector3D::operator/=(double divisor)
{
	x /= divisor;
	y /= divisor;
	z /= divisor;
	return *this;
}

inline bool operator==(const Vector3D &v1, const Vector3D &v2)
{
	return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline bool operator!=(const Vector3D &v1, const Vector3D &v2)
{
	return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
}

inline const Vector3D operator+(const Vector3D &v1, const Vector3D &v2)
{
	return Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline const Vector3D operator-(const Vector3D &v1, const Vector3D &v2)
{
	return Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline const Vector3D operator*(double factor, const Vector3D &vector)
{
	return Vector3D(vector.x * factor, vector.y * factor, vector.z * factor);
}

inline const Vector3D operator*(const Vector3D &vector, double factor)
{
	return Vector3D(vector.x * factor, vector.y * factor, vector.z * factor);
}

inline const Vector3D operator*(const Vector3D &v1, const Vector3D& v2)
{
	return Vector3D(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline const Vector3D operator-(const Vector3D &vector)
{
	return Vector3D(-vector.x, -vector.y, -vector.z);
}

inline const Vector3D operator/(const Vector3D &vector, double divisor)
{
	return Vector3D(vector.x / divisor, vector.y / divisor, vector.z / divisor);
}

inline bool qFuzzyIsNull(const Vector3D& v1)
{
	return qFuzzyIsNull(v1.x) &&
         qFuzzyIsNull(v1.y) &&
         qFuzzyIsNull(v1.z);
}

inline bool qFuzzyCompare(const Vector3D& v1, const Vector3D& v2)
{
	return qFuzzyCompare(v1.x, v2.x) &&
	       qFuzzyCompare(v1.y, v2.y) &&
	       qFuzzyCompare(v1.z, v2.z);
}

inline Vector3D qMin(const Vector3D& v1, const Vector3D& v2)
{
	return Vector3D(qMin(v1.x, v2.x), qMin(v1.y, v2.y), qMin(v1.z, v2.z));
}

inline Vector3D qMax(const Vector3D& v1, const Vector3D& v2)
{
	return Vector3D(qMax(v1.x, v2.x), qMax(v1.y, v2.y), qMax(v1.z, v2.z));
}

inline QVector3D qMin(const QVector3D& v1, const QVector3D& v2)
{
	return QVector3D(qMin(v1.x(), v2.x()), qMin(v1.y(), v2.y()), qMin(v1.z(), v2.z()));
}

inline QVector3D qMax(const QVector3D& v1, const QVector3D& v2)
{
	return QVector3D(qMax(v1.x(), v2.x()), qMax(v1.y(), v2.y()), qMax(v1.z(), v2.z()));
}

inline Vector3D Vector3D::normalized() const
{
	if (qFuzzyIsNull(*this))
		return *this;
  double len = squaredLength();

  if (qFuzzyIsNull(len - 1.0))
      return *this;
   
	return *this / sqrt(len);
}

void Vector3D::normalize()
{
	if (qFuzzyIsNull(*this))
		return;
  
	double len = squaredLength();
  if (qFuzzyIsNull(len-1.0))
		return;

  len = 1.0/sqrt(len);

  x *= len;
  y *= len;
  z *= len;
}

#endif
