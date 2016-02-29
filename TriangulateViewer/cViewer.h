#ifndef _GLWINDOW_H
#define _GLWINDOW_H

#include "iviewer.h"
#include "camera.h"

class QOpenGLShaderProgram;

class CViewer3D : public ISynchroViewer
{
	IScene* scene;
	CFreeCamera freeCamera;
	ICamera3D* camera;

	QPoint lastPos;	
	QColor backgroundColor;
	bool stereo;

public:
	CViewer3D(const QSurfaceFormat &fmt, QWidget* parent);
	virtual ~CViewer3D();

	virtual void setBackgroundColor(const QColor &col) {backgroundColor = col;}

	virtual void setScene(IScene *sc, int cType = _FLAT_CAMERA_){ scene = sc; }	
	virtual void setStereo(int mode) {}// if (format().stereo()) stereo = on; }
	virtual bool isStereo(){return stereo;}

	virtual ICamera3D *getCamera3D() { return camera; }
	virtual void setSceneRadius(double radius) { camera->setSceneRadius(radius); }

	virtual QWidget* widget() { return this; }

	virtual void setCamera(int type) {}

protected:
	virtual void wheelEvent(QWheelEvent *e);
	virtual void keyPressEvent(QKeyEvent *e);
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseDoubleClickEvent(QMouseEvent *e);
	
	virtual void paintGL();
	virtual void resizeGL(int width, int height);
	virtual void initializeGL();

	Vector3D unproject(const QPoint &pnt);
};

#endif
