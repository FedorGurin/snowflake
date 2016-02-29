#include "camera.h"

#define _USE_MATH_DEFINES
#include <cmath>

float separationKoef = 0.8;
float realScreenWidth = 2.0f;
float convergence = 1.0f;
const float interocular = 0.065f;

CPolarCamera::CPolarCamera()
{
	_latitude = 30; //градусы
	_longitude = 45; //градусы
	_altitude = 2.0; 

	_screenWidth = 600;
	_screenHeight = 400;
	fovY = 45.0;
	_sceneRadius = 1;
	
	calcNearFar();
	calcAspect();
	updateProjection();
	updateView();

	//glm::dmat4 projection;
	//glm::dmat4 view;

	//glm::dmat4 projection_x_view;
	
}

void CPolarCamera::updateView()
{
	calcViewMatrix();
}

void CPolarCamera::updateProjection()
{
	calcProjectionMatrix();
}

void CPolarCamera::setScreenWidthAndHeight(int width, int height) 
{ 
	_screenWidth = width; 
	_screenHeight = height; 
	calcAspect();
	updateProjection(); 
}

void CPolarCamera::setSceneRadius(double radius)
{
	_sceneRadius = radius;
	calcNearFar();
	calcProjectionMatrix();
	calcViewMatrix();
}

void CPolarCamera::calcNearFar()
{
	zFar = (_altitude+1.0)*_sceneRadius;
	zNear = zFar/10000.0;
}

void CPolarCamera::calcViewMatrix()
{
	view.setToIdentity();
	view.translate(0.0, 0.0, - _altitude*_sceneRadius);
	view.rotate(_latitude-90.0, 1.0, 0.0, 0.0);
	view.rotate(-_longitude, 0.0, 0.0, 1.0);
	//rightView = leftView = 
	
	//rightView.rotate(-0.5*eyeDistance/(_altitude*_sceneRadius), 0.0, 0.0, 1.0);
	//leftView.rotate(+0.5*eyeDistance/(_altitude*_sceneRadius), 0.0, 0.0, 1.0);
}

void CPolarCamera::calcProjectionMatrix()
{ 
	projection = Matrix4D::perspective(fovY, _aspect, zNear, zFar); 
}

void CPolarCamera::zoom(int delta)
{
	const double wheelSensitivityCoef = 8E-4;
	const double coef = _altitude;
	_altitude += coef*delta*wheelSensitivityCoef;
	calcNearFar();

	updateView();
	updateProjection();
}

void CPolarCamera::rotateAround(QPoint delta)
{
	const double mouseSensitivityCoef = 0.1;
	const double coef = _altitude;

	_latitude += delta.y()*coef*mouseSensitivityCoef;
	_longitude -= delta.x()*coef*mouseSensitivityCoef;
	while (_latitude > 180.0) _latitude -= 360.0;
	while (_longitude > 180.0) _longitude -= 360.0;
	
	while (_latitude < -180.0) _latitude += 360.0;
	while (_longitude < -180.0) _longitude += 360.0;
	//if (_latitude < -90.0) _latitude = -90.0;
	//while (_longitude > 180.0) _longitude -= 360;
	//while (_longitude < -180.0) _longitude += 360;

	updateView();
}

///////////////////////////////
CFreeCamera::CFreeCamera():isOrtho(false)
{
	_viewPoint = Vector3D(0,0,0);
	_position = Vector3D(10,10,10);
	_up = Vector3D(0,0,1);
	Vector3D side = Vector3D::crossProduct(_viewPoint-_position, _up).normalized();
	_up = Vector3D::crossProduct(side, _viewPoint-_position).normalized();

	_screenWidth = 600;
	_screenHeight = 400;
	fovY = 45.0;
	_sceneRadius = 100;
	_sceneCenter = Vector3D(0,0,0);
	stereoMode = 0;
	realEyeSeparation = interocular/realScreenWidth;
	separation = separationKoef*realEyeSeparation;

	
	calcAspect();
	updateNearFar();
	updateView();
}

void CFreeCamera::setScreenWidthAndHeight(int width, int height) 
{ 
	_screenWidth = width; 
	_screenHeight = height; 
	calcAspect();
	updateProjection(); 
}

void CFreeCamera::setSceneRadius(double radius)
{
	_sceneRadius = radius;
	//qDebug() << "sceneRadius" << _sceneCenter << _sceneRadius;
	updateNearFar();
}

void CFreeCamera::calcNearFar()
{
	if (isOrtho)
		return;
	double distToCenter = (_position - _sceneCenter).length();
	zFar = distToCenter+_sceneRadius;
	zNear = qMax(zFar/1000.0, distToCenter-_sceneRadius);
}

void CFreeCamera::setOrtho(double left, double right, double bottom, double top, double near, double far)
{
	isOrtho = true;
	zNear = near;
	zFar = far;
	projection = Matrix4D::ortho(left, right, bottom, top, near, far);
}
void CFreeCamera::setPerspective(double fovy, double near, double far)
{
	isOrtho = false;
	calcNearFar();
	calcProjectionMatrix();
}

void CFreeCamera::calcViewMatrix()
{
	view = Matrix4D::lookAt(_position, _viewPoint, _up);
	
	//{
	//	QVector3D dir = _viewPoint - _position;
	//	QVector3D side = QVector3D::crossProduct(dir, _up);
	//	side.normalize();
	//	side *= eyeDistance/2.0;
	//	QVector3D leftPos = _position - side;
	//	QVector3D rightPos = _position + side;
	//	
	//	rightView.setToIdentity();
	//	leftView.setToIdentity();
	//	rightView.lookAt(_position + side, _viewPoint, _up);
	//	leftView.lookAt(_position - side, _viewPoint, _up);
	//}
	//else
	//{

	//}
}

void CFreeCamera::calcProjectionMatrix()
{ 
	if (isOrtho)
		return;
	projection = Matrix4D::perspective(fovY, _aspect, zNear, zFar); 
	if (stereoMode)
	{
		leftProjection = rightProjection = projection;
		//int side = rightEye ? 1 : -1; 
		leftProjection(0, 2) += separation;
		leftProjection(0, 3) = separation*convergence;
		leftProjection.optimize();

		rightProjection(0, 2) -= separation;
		rightProjection(0, 3) = -separation*convergence;
		rightProjection.optimize();
	}
	
}

void CFreeCamera::zoom(int delta)
{
	const double wheelSensitivityCoef = 8E-4;
	Vector3D viewDir = _viewPoint-_position;
	double dist = viewDir.length();
	viewDir /= dist;
	const double coef = dist;
	dist += coef*delta*wheelSensitivityCoef;
	viewDir *= dist;
	_position = _viewPoint - viewDir;
	
	updateNearFar();
	updateView();
}

void CFreeCamera::rotateAround(QPoint delta)
{
	const double mouseSensitivityCoef = 0.05;	
	
	Vector3D dir = _position-_viewPoint;
	double dist = dir.length();
	//const double coef = dist;
	Vector3D side = Vector3D::crossProduct(dir, _up);

	QuaternionD hRot = QuaternionD::fromAxisAndAngle(Vector3D(0,0,1), -mouseSensitivityCoef*delta.x()).normalized();

	dir = hRot.rotatedVector(dir);
	side = hRot.rotatedVector(side);
	//_up = hRot.rotatedVector(_up);

	QuaternionD vRot = QuaternionD::fromAxisAndAngle(side, mouseSensitivityCoef*delta.y()).normalized();

	dir = vRot.rotatedVector(dir);
	_up = Vector3D::crossProduct(side, dir).normalized();

	_position = _viewPoint + dir;

	updateNearFar();
	updateView();
}

void CFreeCamera::rotateHead(const QPoint &delta)
{
	const double mouseSensitivityCoef = 0.05;	

	Vector3D dir = _position-_viewPoint;
	double dist = dir.length();
	//const double coef = dist;
	Vector3D side = Vector3D::crossProduct(dir, _up);

	QuaternionD hRot = QuaternionD::fromAxisAndAngle(Vector3D(0,0,1), -mouseSensitivityCoef*delta.x()).normalized();

	dir = hRot.rotatedVector(dir);
	side = hRot.rotatedVector(side);
	//_up = hRot.rotatedVector(_up);

	QuaternionD vRot = QuaternionD::fromAxisAndAngle(side, mouseSensitivityCoef*delta.y()).normalized();

	dir = vRot.rotatedVector(dir);
	_up = Vector3D::crossProduct(side, dir).normalized();

	_viewPoint	= _position - dir;

	updateView();

}

void CFreeCamera::changeViewPoint(const Vector3D &p)
{
	_viewPoint = p;
	_up = Vector3D(0,0,1);
	Vector3D side = Vector3D::crossProduct(_viewPoint-_position, _up).normalized();
	_up = Vector3D::crossProduct(side, _viewPoint-_position).normalized();
	updateView();
}

void CFreeCamera::unproject(const QPoint &pnt)
{

}

void CFreeCamera::moveSide(float delta)
{
	const double mouseSensitivityCoef = 0.01;	

	Vector3D dir = _viewPoint - _position;
	double dist = dir.length();
	dir = dir/dist;
	const double coef = dist;
	Vector3D side = Vector3D::crossProduct(dir, _up);
	_viewPoint += side*delta*mouseSensitivityCoef*dist;
	_position += side*delta*mouseSensitivityCoef*dist;
	
	updateNearFar();
	updateView();
}

void CFreeCamera::showAll()
{
	_viewPoint = _sceneCenter;
	_position = _sceneCenter+ Vector3D(2,2,2)*_sceneRadius;

	_up = Vector3D(0,0,1);
	Vector3D side = Vector3D::crossProduct(_viewPoint-_position, _up).normalized();
	_up = Vector3D::crossProduct(side, _viewPoint-_position).normalized();

	updateNearFar();
	updateView();
}

void CFreeCamera::moveVert(float delta)
{
	const double mouseSensitivityCoef = 0.01;	

	Vector3D dir = _viewPoint - _position;
	double dist = dir.length();
	dir = dir/dist;
	const double coef = dist;
	_viewPoint += _up*delta*mouseSensitivityCoef*dist;
	_position += _up*delta*mouseSensitivityCoef*dist;
	updateNearFar();
	updateView();
}
void CFreeCamera::moveForward(float delta)
{
	const double mouseSensitivityCoef = 0.01;	

	Vector3D dir = _viewPoint - _position;
	double dist = dir.length();
	dir = dir/dist;
	const double coef = dist;
	_viewPoint += dir*delta*mouseSensitivityCoef*dist;
	_position += dir*delta*mouseSensitivityCoef*dist;
	updateNearFar();
	updateView();
}