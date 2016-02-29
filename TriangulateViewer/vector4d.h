#ifndef VECTOR4D_H
#define VECTOR4D_H

#include <QVector4D>
#include "vector3d.h"

class Vector4D
{
	public:
	union
	{
		struct
		{
			double x,y,z,w;
		};
		struct
		{
			double r,g,b,a;
		};
		double v[4];
	};
	Vector4D(): x(0), y(0), z(0), w(1.0) {}
	Vector4D(double _x, double _y, double _z, double _w = 1.0): x(_x), y(_y), z(_z), w(_w) {}
	Vector4D(const Vector4D &vec): x(vec.x), y(vec.y), z(vec.z), w(vec.w) {}
	explicit Vector4D(const QVector4D &vec): x(vec.x()), y(vec.y()), z(vec.z()), w(vec.w()) {}
	Vector4D(const Vector3D &vec, double _w = 1.0): x(vec.x), y(vec.y), z(vec.z), w(_w) {}
	explicit Vector4D(const double vec[4]): x(vec[0]), y(vec[1]), z(vec[2]), w(vec[3]) {}	
	Vector3D toVector3D() { return Vector3D(x,y,z); }
	QVector4D toQVector4D() { return QVector4D(x,y,z,w); }
};

#endif
