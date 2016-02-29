#ifndef IVIEWER_H
#define IVIEWER_H

#include "iCamera.h"
#include <QOpenGLWidget>
#include "gl_version.h"


class QOpenGLShaderProgram;
class Matrix4D;
class IScene;
class QSurfaceFormat;


class IViewer
{
public:
	static IViewer* create(QWidget* p, bool stereo = false, bool stencil = false);
	virtual ~IViewer(){}	

	virtual void setScene(IScene* sc, int cType = _FLAT_CAMERA_)			= 0;
	virtual void setStereo(int mode)  																= 0;
	virtual bool isStereo()																						= 0;
	virtual bool stereoFormat() const																	= 0;

	virtual ICamera3D *getCamera3D()                     							= 0;

	virtual void setCamera(int type)                                  = 0;
	virtual void setSceneRadius(double radius)												= 0;
	virtual QWidget *getWidget()																			= 0;
	//virtual QGLContext * context() const                              = 0;	
	virtual void setBackgroundColor(const QColor &col)								= 0;	
	virtual void makeCurrent()																				= 0;
	virtual void updateGL()																						= 0;
};

class ISynchroViewer : public QOpenGLWidget,  public IViewer
{
	Q_OBJECT

protected:
	bool invalid;
	QAbstractOpenGLFunctions *glFunc;

public:
	ISynchroViewer(const QSurfaceFormat &fmt, QWidget *parent, int ms = 40);
	virtual ~ISynchroViewer() {}

	virtual void updateGL() { invalid = true; }
	virtual void makeCurrent() { gl = static_cast<Current_OpenGL_Version_Profile *>(glFunc); QOpenGLWidget::makeCurrent(); }

protected:
	virtual QWidget* getWidget() { return this; }
	//virtual QGLContext* context() const { return QGLWidget::context(); }
	virtual bool stereoFormat() const { return format().stereo(); }
	
private slots:
	void onTimerTick();
};

#endif
