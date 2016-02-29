#include <QDataStream>
#include "vector3d.h"

QDataStream &operator<<(QDataStream &stream, const Vector3D &vector)
{
	stream << vector.x << vector.y << vector.z;
	return stream;
}

QDataStream &operator>>(QDataStream &stream, Vector3D &vector)
{
	stream >> vector.x >> vector.y >> vector.z;
	return stream;
}
