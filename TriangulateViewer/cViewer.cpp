#include <QtGui>
#include <QMessageBox>
#include <QOpenGLShaderProgram>
#include <GL/glu.h>
#include "cViewer.h"
#include "iScene.h"

Current_OpenGL_Version_Profile *gl = 0;

IViewer* IViewer::create(QWidget* parent, bool stereo, bool stencil)
{
	QSurfaceFormat fmt;
	//fmt.setRenderableType(QSurfaceFormat::OpenGL);
	//fmt.setDepthBufferSize(24);

	if (stereo)
		fmt.setStereo(true);
	if (stencil)
		fmt.setStencilBufferSize(8);

	fmt.setSamples(4);
	setVersionProfile(fmt);
	QSurfaceFormat::setDefaultFormat(fmt);
	//qDebug() << fmt;
	return new CViewer3D(fmt, parent);
}


ISynchroViewer::ISynchroViewer(const QSurfaceFormat &fmt, QWidget *parent, int msec): QOpenGLWidget(parent), invalid(false)
{
	//setFormat(fmt);
	QTimer* timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(onTimerTick()));
  timer->start(msec);
}

void ISynchroViewer::onTimerTick()
{
	if (invalid)
	{
		if (isVisible())
			update();			
		invalid = false;
	}
}

/////////////////////////////////////////////
CViewer3D::CViewer3D(const QSurfaceFormat &fmt, QWidget *parent): ISynchroViewer(fmt, parent), scene(NULL)
{
	qDebug() << format();
	setFocusPolicy(Qt::StrongFocus);
	//setAutoBufferSwap(true);
	setMouseTracking(true);
	backgroundColor = Qt::darkGray;
	stereo = format().stereo();
	camera = &freeCamera;
}

CViewer3D::~CViewer3D()
{
}

void CViewer3D::initializeGL()
{
	glFunc = context()->versionFunctions<Current_OpenGL_Version_Profile>();
	if(!glFunc)
	{
		QMessageBox::information(this, ("Предупреждение"), ("Видеоадаптер не поддерживает OpenGL 4.3. Приложение будет закрыто"));
    exit(1);
	}
	glFunc->initializeOpenGLFunctions();
	
	gl = static_cast<Current_OpenGL_Version_Profile *>(glFunc);

	
	gl->glEnable(GL_DEPTH_TEST);
	gl->glDepthFunc(GL_LEQUAL);
	gl->glDisable(GL_CULL_FACE);
	gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	gl->glClearColor(backgroundColor.redF(), backgroundColor.greenF(), backgroundColor.blueF(), 1);
	
	//gl->glEnable(GL_DEPTH_CLAMP);
	//gl->glDisable(GL_DEPTH_CLAMP);

	if (scene)
		scene->initGL(this);
}

void CViewer3D::resizeGL(int width, int height)
{
	gl = static_cast<Current_OpenGL_Version_Profile *>(glFunc);
	camera->setScreenWidthAndHeight(width, height);	
}

void CViewer3D::paintGL()
{	
	//gl->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());
	gl = static_cast<Current_OpenGL_Version_Profile *>(glFunc);
	gl->glClearColor(backgroundColor.redF(), backgroundColor.greenF(), backgroundColor.blueF(), 1);
	gl->glClearDepth(1.0);
	gl->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//gl->glViewport(0, 0, camera->screenWidth(), camera->screenHeight());	
	if(!scene)
		return;

	if (!stereo)
	{
		scene->draw();
	}
	else
	{
		gl->glDrawBuffer(GL_RIGHT);		
		scene->draw(1);		
	  
	  gl->glDrawBuffer(GL_LEFT);
	  scene->draw(2);		
	}
	//swapBuffers();
}


void CViewer3D::wheelEvent(QWheelEvent *e)
{
	camera->zoom(e->delta());
	updateGL();
}

void CViewer3D::mousePressEvent(QMouseEvent *e)
{
  lastPos = e->pos();
}

void CViewer3D::mouseMoveEvent(QMouseEvent *e)
{
	if (e->buttons() & Qt::LeftButton)
	{
		camera->rotate(e->pos()-lastPos);
		//polarCamera.rotateAround(e->pos()-lastPos);
		lastPos = e->pos();
		
		updateGL();			
	}	
	else if (e->buttons() & Qt::RightButton)
	{
		//camera->rotate(e->pos()-lastPos);
		freeCamera.rotateHead(e->pos()-lastPos);
		lastPos = e->pos();
		
		updateGL();			
	}	
}

Vector3D CViewer3D::unproject(const QPoint &pnt)
{  
	GLint viewport[4] = { 0, 0, camera->screenWidth(), camera->screenHeight() };
  //GLdouble mvmatrix[16], projmatrix[16];
  //GLfloat color[3];

  GLint vx, vy; 
  GLfloat vz;
  makeCurrent();
            
  vx = pnt.x();
  vy = viewport[3] - 1 - pnt.y();
  gl->glReadPixels(vx, vy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &vz);

	double x,y,z;

	gluUnProject((GLdouble) vx, (GLdouble) vy, (GLdouble) vz, camera->getViewMatrix().constData(), camera->getProjectionMatrix().constData(), viewport, &x, &y, &z);
  return Vector3D(x,y,z);   
}

void CViewer3D::mouseDoubleClickEvent(QMouseEvent *e)
{
	freeCamera.changeViewPoint(unproject(e->pos()));
	updateGL();
}

void CViewer3D::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_Left)
	{
		freeCamera.moveSide(-1);
		updateGL();
	}
	else if (e->key() == Qt::Key_Right)
	{
		freeCamera.moveSide(1);
		updateGL();
	}
	else if (e->key() == Qt::Key_Up)
	{
		freeCamera.moveVert(1);
		updateGL();
	}
	else if (e->key() == Qt::Key_Down)
	{
		freeCamera.moveVert(-1);
		updateGL();
	}
	else if (e->key() == Qt::Key_Home)
	{
		freeCamera.showAll();
		updateGL();
	}
}
