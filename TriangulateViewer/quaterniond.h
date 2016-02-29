#ifndef QUATERNIOND_H
#define QUATERNIOND_H

#include "vector3d.h"
#include "vector4d.h"
#include <QtCore/qmath.h>

#include <QQuaternion>

class QuaternionD
{
public:
	QuaternionD(): wp(1.0f), xp(0.0f), yp(0.0f), zp(0.0f) {}
	QuaternionD(double scalar, double xpos, double ypos, double zpos) : wp(scalar), xp(xpos), yp(ypos), zp(zpos) {}
	QuaternionD(double scalar, const Vector3D& vector): wp(scalar), xp(vector.x), yp(vector.y), zp(vector.z) {}
	explicit QuaternionD(const Vector4D& vector): wp(vector.w), xp(vector.x), yp(vector.y), zp(vector.z) {}
	
	bool isNull() const;
	bool isIdentity() const;
	
	Vector3D vector() const { return Vector3D(xp, yp, zp); }
	void setVector(const Vector3D& vector);
	void setVector(double x, double y, double z) { xp = x; yp = y; zp = z; } 
	
	Vector4D toVector4D() const { return Vector4D(xp, yp, zp, wp); }
	QQuaternion toQQuaternion() const { return QQuaternion(wp, xp, yp, zp); }
	
	double x() const { return xp; }
	double y() const { return yp; }
	double z() const { return zp; }
	double scalar() const { return wp; }
	
	void setX(double x) { xp = x; }
	void setY(double y) { yp = y; }
	void setZ(double z) { zp = z; }
	void setScalar(double scalar) { wp = scalar; }
	
	double length() const  { return sqrt(xp * xp + yp * yp + zp * zp + wp * wp); }
	double lengthSquared() const { return xp * xp + yp * yp + zp * zp + wp * wp; }
	
	inline QuaternionD normalized() const;
	inline void normalize();
	
	inline QuaternionD conjugate() const { return QuaternionD(wp, -xp, -yp, -zp); }
	
	inline Vector3D rotatedVector(const Vector3D& vector) const;
	
	inline QuaternionD &operator+=(const QuaternionD &quaternion);
	inline QuaternionD &operator-=(const QuaternionD &quaternion);
	inline QuaternionD &operator*=(double factor);
	inline QuaternionD &operator*=(const QuaternionD &quaternion);
	inline QuaternionD &operator/=(double divisor);
	
	friend inline bool operator==(const QuaternionD &q1, const QuaternionD &q2);
	friend inline bool operator!=(const QuaternionD &q1, const QuaternionD &q2);
	friend inline const QuaternionD operator+(const QuaternionD &q1, const QuaternionD &q2);
	friend inline const QuaternionD operator-(const QuaternionD &q1, const QuaternionD &q2);
	friend inline const QuaternionD operator*(double factor, const QuaternionD &quaternion);
	friend inline const QuaternionD operator*(const QuaternionD &quaternion, double factor);
	friend inline const QuaternionD operator*(const QuaternionD &q1, const QuaternionD& q2);
	friend inline const QuaternionD operator-(const QuaternionD &quaternion);
	friend inline const QuaternionD operator/(const QuaternionD &quaternion, double divisor);
	
	friend inline bool qFuzzyCompare(const QuaternionD& q1, const QuaternionD& q2);
	friend inline bool qFuzzyIsNull(const QuaternionD& q1);
	
	    //operator QVariant() const;
	
	static inline QuaternionD fromAxisAndAngle(const Vector3D& axis, double angle);
	static inline QuaternionD fromAxisAndAngle(double x, double y, double z, double angle);
	
	static inline QuaternionD slerp(const QuaternionD& q1, const QuaternionD& q2, double t);
	static inline QuaternionD nlerp(const QuaternionD& q1, const QuaternionD& q2, double t);
	
private:
	double wp, xp, yp, zp;
};

inline bool QuaternionD::isNull() const
{
	return qIsNull(xp) && qIsNull(yp) && qIsNull(zp) && qIsNull(wp);
}

inline bool QuaternionD::isIdentity() const
{
	return qIsNull(xp) && qIsNull(yp) && qIsNull(zp) && wp == 1.0f;
}

inline QuaternionD &QuaternionD::operator+=(const QuaternionD &quaternion)
{
	xp += quaternion.xp;
	yp += quaternion.yp;
	zp += quaternion.zp;
	wp += quaternion.wp;
	return *this;
}

inline QuaternionD &QuaternionD::operator-=(const QuaternionD &quaternion)
{
	xp -= quaternion.xp;
	yp -= quaternion.yp;
	zp -= quaternion.zp;
	wp -= quaternion.wp;
	return *this;
}

inline QuaternionD &QuaternionD::operator*=(double factor)
{
	xp *= factor;
	yp *= factor;
	zp *= factor;
	wp *= factor;
	return *this;
}

inline const QuaternionD operator*(const QuaternionD &q1, const QuaternionD& q2)
{
	double ww = (q1.zp + q1.xp) * (q2.xp + q2.yp);
	double yy = (q1.wp - q1.yp) * (q2.wp + q2.zp);
	double zz = (q1.wp + q1.yp) * (q2.wp - q2.zp);
	double xx = ww + yy + zz;
	double qq = 0.5 * (xx + (q1.zp - q1.xp) * (q2.xp - q2.yp));
	
	double w = qq - ww + (q1.zp - q1.yp) * (q2.yp - q2.zp);
	double x = qq - xx + (q1.xp + q1.wp) * (q2.xp + q2.wp);
	double y = qq - yy + (q1.wp - q1.xp) * (q2.yp + q2.zp);
	double z = qq - zz + (q1.zp + q1.yp) * (q2.wp - q2.xp);
	
	return QuaternionD(w, x, y, z);
}

inline QuaternionD &QuaternionD::operator*=(const QuaternionD &quaternion)
{
	*this = *this * quaternion;
	return *this;
}

inline QuaternionD &QuaternionD::operator/=(double divisor)
{
	xp /= divisor;
	yp /= divisor;
	zp /= divisor;
	wp /= divisor;
	return *this;
}

inline bool operator==(const QuaternionD &q1, const QuaternionD &q2)
{
	return q1.xp == q2.xp && q1.yp == q2.yp && q1.zp == q2.zp && q1.wp == q2.wp;
}

inline bool operator!=(const QuaternionD &q1, const QuaternionD &q2)
{
	return q1.xp != q2.xp || q1.yp != q2.yp || q1.zp != q2.zp || q1.wp != q2.wp;
}

inline const QuaternionD operator+(const QuaternionD &q1, const QuaternionD &q2)
{
	return QuaternionD(q1.wp + q2.wp, q1.xp + q2.xp, q1.yp + q2.yp, q1.zp + q2.zp);
}

inline const QuaternionD operator-(const QuaternionD &q1, const QuaternionD &q2)
{
	return QuaternionD(q1.wp - q2.wp, q1.xp - q2.xp, q1.yp - q2.yp, q1.zp - q2.zp);
}

inline const QuaternionD operator*(double factor, const QuaternionD &quaternion)
{
	return QuaternionD(quaternion.wp * factor, quaternion.xp * factor, quaternion.yp * factor, quaternion.zp * factor);
}

inline const QuaternionD operator*(const QuaternionD &quaternion, double factor)
{
	return QuaternionD(quaternion.wp * factor, quaternion.xp * factor, quaternion.yp * factor, quaternion.zp * factor);
}

inline const QuaternionD operator-(const QuaternionD &quaternion)
{
	return QuaternionD(-quaternion.wp, -quaternion.xp, -quaternion.yp, -quaternion.zp);
}

inline const QuaternionD operator/(const QuaternionD &quaternion, double divisor)
{
	return QuaternionD(quaternion.wp / divisor, quaternion.xp / divisor, quaternion.yp / divisor, quaternion.zp / divisor);
}

inline bool qFuzzyCompare(const QuaternionD& q1, const QuaternionD& q2)
{
	return qFuzzyCompare(q1.xp, q2.xp) &&
	       qFuzzyCompare(q1.yp, q2.yp) &&
	       qFuzzyCompare(q1.zp, q2.zp) &&
	       qFuzzyCompare(q1.wp, q2.wp);
}

inline bool qFuzzyIsNull(const QuaternionD& q1)
{
	return qFuzzyIsNull(q1.xp) &&
	       qFuzzyIsNull(q1.yp) &&
	       qFuzzyIsNull(q1.zp) &&
	       qFuzzyIsNull(q1.wp);
}

inline void QuaternionD::setVector(const Vector3D& aVector)
{
	xp = aVector.x;
	yp = aVector.y;
	zp = aVector.z;
}
   
inline QuaternionD QuaternionD::normalized() const
{
	if (qFuzzyIsNull(*this))
		return QuaternionD(0.0f, 0.0f, 0.0f, 0.0f);
	double len = lengthSquared();
	if (qFuzzyIsNull(len - 1.0))
	    return *this;
	else
		return *this * (1.0/sqrt(len));	    
}

/*!
    Normalizes the currect quaternion in place.  Nothing happens if this
    is a null quaternion or the length of the quaternion is very close to 1.

    \sa length(), normalized()
*/
inline void QuaternionD::normalize()
{
	if (qFuzzyIsNull(*this))
		return;
	
	double len = lengthSquared();
	if (qFuzzyIsNull(len - 1.0))
	    return;
	
	len = 1.0/sqrt(len);
	
	xp *= len;
	yp *= len;
	zp *= len;
	wp *= len;
}

inline Vector3D QuaternionD::rotatedVector(const Vector3D& vector) const
{
	return (*this * QuaternionD(0, vector) * conjugate()).vector();
}

inline QuaternionD QuaternionD::fromAxisAndAngle(const Vector3D& axis, double angle)
{
	// Algorithm from:
	// http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q56
	// We normalize the result just in case the values are close
	// to zero, as suggested in the above FAQ.
	double a = (angle / 2.0f) * M_PI / 180.0f;
	double s = sin(a);
	double c = cos(a);
	Vector3D ax = axis.normalized();
	return QuaternionD(c, ax.x * s, ax.y * s, ax.z * s).normalized();
}

/*!
    Creates a normalized quaternion that corresponds to rotating through
    \a angle degrees about the 3D axis (\a x, \a y, \a z).
*/
inline QuaternionD QuaternionD::fromAxisAndAngle(double x, double y, double z, double angle)
{
	if (!qFuzzyIsNull(Vector3D(x,y,z)))
	{
		double length = x * x + y * y + z * z;
		if (!qFuzzyIsNull(length - 1.0)) 
		{
			length = 1.0/sqrt(length);
			x *= length;
			y *= length;
			z *= length;
		}
	}
	double a = (angle / 2.0) * M_PI / 180.0;
	double s = sin(a);
	double c = cos(a);
	return QuaternionD(c, x * s, y * s, z * s).normalized();
}

QuaternionD QuaternionD::slerp(const QuaternionD& q1, const QuaternionD& q2, double t)
{
	// Handle the easy cases first.
	if (t <= 0.0f)
	    return q1;
	else if (t >= 1.0f)
	    return q2;
	
	// Determine the angle between the two quaternions.
	QuaternionD q2b;
	double dot;
	dot = q1.xp * q2.xp + q1.yp * q2.yp + q1.zp * q2.zp + q1.wp * q2.wp;
	if (dot >= 0.0f) {
	    q2b = q2;
	} else {
	    q2b = -q2;
	    dot = -dot;
	}
	
	// Get the scale factors.  If they are too small,
	// then revert to simple linear interpolation.
	double factor1 = 1.0f - t;
	double factor2 = t;
	if ((1.0 - dot) > 0.0000001) {
	    double angle = acos(dot);
	    double sinOfAngle = sin(angle);
	    if (sinOfAngle > 0.0000001) 
			{
	        factor1 = sin((1.0 - t) * angle) / sinOfAngle;
	        factor2 = sin(t * angle) / sinOfAngle;
	    }
	}
	
	// Construct the result quaternion.
	return q1 * factor1 + q2b * factor2;
}

/*!
    Interpolates along the shortest linear path between the rotational
    positions \a q1 and \a q2.  The value \a t should be between 0 and 1,
    indicating the distance to travel between \a q1 and \a q2.
    The result will be normalized().

    If \a t is less than or equal to 0, then \a q1 will be returned.
    If \a t is greater than or equal to 1, then \a q2 will be returned.

    The nlerp() function is typically faster than slerp() and will
    give approximate results to spherical interpolation that are
    good enough for some applications.

    \sa slerp()
*/
inline QuaternionD QuaternionD::nlerp(const QuaternionD& q1, const QuaternionD& q2, double t)
{
	// Handle the easy cases first.
	if (t <= 0.0f)
	    return q1;
	else if (t >= 1.0f)
	    return q2;
	
	// Determine the angle between the two quaternions.
	QuaternionD q2b;
	double dot;
	dot = q1.xp * q2.xp + q1.yp * q2.yp + q1.zp * q2.zp + q1.wp * q2.wp;
	if (dot >= 0.0f)
	    q2b = q2;
	else
	    q2b = -q2;
	
	// Perform the linear interpolation.
	return (q1 * (1.0f - t) + q2b * t).normalized();
}

#endif
