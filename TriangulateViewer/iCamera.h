#ifndef _I_CAMERA_H_
#define _I_CAMERA_H_

#include <QPointF>
#include <QRectF>
#include "matrix4d.h"

enum { _FLAT_CAMERA_ = 0, _POLAR_CAMERA_, _FREE_CAMERA_ };
enum {_Right_Eye_ = 1, _Left_Eye_ = 2 };

class ICamera
{
public:
	virtual double aspect() const															  = 0;
	virtual int screenWidth() const															= 0;
	virtual int screenHeight() const														= 0;
	virtual void setScreenWidthAndHeight(int width, int height) = 0;
};

class ICamera3D: public ICamera
{
public:
	virtual double latitude() const												= 0; //радианы
	virtual double longitude() const											= 0; //радианы
	virtual double altitude() const												= 0; //метры над уровнем моря	
	virtual double normalizedHeight() const								= 0; //altitude/earthRadius


	virtual double getMinVisibleRange() const { return 0; }
	virtual double getMaxVisibleRange() const { return 0; }
	virtual double getLevel() const { return -1; }

	virtual double fov() const														{ return 0; } //градусы
	virtual double nearPlane() const											= 0;
	virtual double farPlane() const												= 0;
	virtual const Matrix4D &getProjectionMatrix() const	= 0;
	virtual const Matrix4D &getViewMatrix() const				= 0;
	virtual double sceneRadius() const										{ return 0;}

	virtual void showAll() {}
	virtual void setSceneRadius(double radius)						{}
	virtual void zoom(int delta)													{}
	virtual void rotate(const QPoint &delta) {}
	virtual void rotateHead(const QPoint &delta) {}
	
	virtual void moveSide(float delta) {}
	virtual void moveVert(float delta) {}
	virtual void moveForward(float delta) {}

	virtual void setStereo(int mode)											{}
	virtual void setEye(int eye)													{}

	virtual void setPosition(const Vector3D &pos) {}
	virtual void setViewPoint(const Vector3D &pos) {} 
	virtual void setUpDir(const Vector3D &pos)				{}

	virtual void setLatLonAlt(double lat, double lon, double h) {} //градусы, градусы, метры
	virtual void setHeadPitchRoll(double head, double pitch, double roll) {}//градусы, градусы, градусы
	virtual void update() {}

	virtual Vector3D position() const { return Vector3D(); }
	virtual Vector3D viewPoint() const { return Vector3D(); }
	virtual Vector3D upDir() const { return Vector3D(); }
	virtual void setOrtho(double, double, double, double, double, double) {}
	virtual void setPerspective(double, double, double) {}
};

#endif
