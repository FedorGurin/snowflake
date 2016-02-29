#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "iCamera.h"
 
class CPolarCamera: public ICamera3D
{
	double _latitude; 
	double _longitude;
	double _altitude;  // расстояние от камеры до центра сцены

	//glm::dvec3 pos; //радианы, высота в радиусах земли
	//glm::dvec3 rap;
	
	Matrix4D projection;	
	Matrix4D view;
	Matrix4D leftView;
	Matrix4D rightView;

	//glm::dmat4 projection_x_view;
	
	double zNear;
	double zFar;
	
	int _screenWidth;
	int _screenHeight;
	double _aspect;
	
	double fovY; //радианы

	double _sceneRadius;
	float separation;
	float realEyeSeparation;

public:
	CPolarCamera();
	virtual ~CPolarCamera() {}

	virtual double aspect() const { return _aspect; }
	virtual double latitude() const { return _latitude; }
	virtual double longitude() const { return _longitude; } 
	virtual double altitude() const { return _altitude*_sceneRadius; }
	virtual double normalizedHeight() const { return _altitude; }
	virtual double fov() const { return fovY; }
	virtual double nearPlane() const { return zNear; }
	virtual double farPlane() const { return zFar; }
	virtual int screenWidth() const { return _screenWidth; }
	virtual int screenHeight() const { return _screenHeight; }

	virtual const Matrix4D &getProjectionMatrix() const { return projection; }
	virtual const Matrix4D &getViewMatrix() const { return view; }

	virtual void setScreenWidthAndHeight(int width, int height);
	virtual void setSceneRadius(double radius);
	virtual double sceneRadius() const { return _sceneRadius; }
	virtual void zoom(int delta);
	virtual void rotate(QPoint delta) { rotateAround(delta); }
	void rotateAround(QPoint delta);
	
private:
	void calcAspect() { _aspect = _screenWidth/double(_screenHeight); }
	void calcNearFar();
	void calcProjectionMatrix();
	void calcViewMatrix();
	void updateView();
	void updateProjection();

	friend class CViewer3D;	
};

class CFreeCamera: public ICamera3D
{
	Vector3D _viewPoint;
	Vector3D _position;
	Vector3D _up;

	Matrix4D projection;	
	Matrix4D rightProjection;
	Matrix4D leftProjection;
	Matrix4D view;
	
	double zNear;
	double zFar;
	
	int _screenWidth;
	int _screenHeight;
	double _aspect;
	
	double fovY; //радианы

	double _sceneRadius;
	Vector3D _sceneCenter;
	
	int stereoMode;
	bool rightEye;

	float separation;
	float realEyeSeparation;
	bool isOrtho;

public:
	CFreeCamera();
	virtual ~CFreeCamera() {}

	virtual double aspect() const { return _aspect; }
	virtual double latitude() const { return 0; }
	virtual double longitude() const { return 0; } 
	virtual double altitude() const { return 0; }
	virtual double normalizedHeight() const { return 0; }

	virtual double fov() const { return fovY; }
	virtual double nearPlane() const { return zNear; }
	virtual double farPlane() const { return zFar; }
	virtual int screenWidth() const { return _screenWidth; }
	virtual int screenHeight() const { return _screenHeight; }

	virtual const Matrix4D &getProjectionMatrix() const { return stereoMode ? (rightEye ? rightProjection : leftProjection) : projection; }
	virtual const Matrix4D &getViewMatrix() const { return view; }

	virtual void setScreenWidthAndHeight(int width, int height);
	virtual void setSceneRadius(double radius);
	virtual double sceneRadius() const { return _sceneRadius; }
	virtual void zoom(int delta);
	virtual void rotate(const QPoint &delta) { rotateAround(delta); }
	virtual void rotateHead(const QPoint &delta);
	virtual void moveSide(float delta);
	virtual void moveVert(float delta);
	virtual void moveForward(float delta);
	virtual void showAll();

	virtual void setPosition(const Vector3D &pos) { _position = pos; updateNearFar(); updateView(); }
	virtual void setViewPoint(const Vector3D &pos) { _viewPoint = pos; updateNearFar();  updateView(); }
	virtual void setUpDir(const Vector3D &pos)					{ _up = pos; updateNearFar(); updateView(); }

	virtual Vector3D position() const { return _position; }
	virtual Vector3D viewPoint() const { return _viewPoint; }
	virtual Vector3D upDir() const { return _up; }

	virtual void setStereo(int mode) { stereoMode = mode; }
	virtual void setEye(int eye) { rightEye = (eye == _Right_Eye_); }
	virtual void setOrtho(double, double, double, double, double, double);
	virtual void setPerspective(double, double, double);
	
	void changeViewPoint(const Vector3D &p);
	void unproject(const QPoint &pnt);
	
private:
	void rotateAround(QPoint delta);
	void calcAspect() { _aspect = _screenHeight ? _screenWidth/double(_screenHeight) : 1; }
	void calcNearFar();
	void calcProjectionMatrix();
	void calcViewMatrix();
	void updateProjection() { calcProjectionMatrix(); }
	void updateView() { calcViewMatrix(); }
	void updateNearFar() { calcNearFar(); calcProjectionMatrix(); }

	friend class CViewer3D;	
};

#endif
