#ifndef _SCENE3D_H_
#define _SCENE3D_H_

#include "iScene.h"
#include "matrix4d.h"
#include "vector3d.h"
#include "timeMeasurer.h"
#include <QTime>

#pragma pack(push)
#pragma pack(1)

class QOpenGLShaderProgram;

class CScene3D : public IScene
{
	//Служебные
	QOpenGLShaderProgram *prg;
	quint32 vb;
	quint32 ib;
	GLuint tb0, tb1, tb2, tb3, tb4;
	quint32 vaoObject;
	IViewer* viewer;
	float sceneRadius;
	QVector4D lightDirection;
	int noiseSize;
	QTime time1;

	TimeMeasurer * timeMeasurer;

public:
	CScene3D(): viewer(0),sceneRadius(1),lightDirection(0, 0, 1, 0),prg(0),vb(0),ib(0){}
	
	virtual ~CScene3D();

	virtual void draw(int flags = 0);

	virtual void showBestView();
	virtual void update(){}
	virtual void initGL(IViewer *view);	
	virtual void clearGL();
	virtual bool ready() { return viewer; }
	virtual void clear();
	virtual void setParm(const QVector <float> &parms, const QString &str);

private:
	CScene3D(const CScene3D& sc);
};

#pragma pack(pop)

#endif
