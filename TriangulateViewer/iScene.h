#ifndef _I_SCENE_H_
#define _I_SCENE_H_

#include <QVector>
#include <QPoint>
#include "iGLResource.h"

class IScene : public IGLResource
{
public:	
	virtual ~IScene() {}

	virtual void draw(int flags = 0) = 0;
	virtual void update()						 = 0;
	virtual void clear()						 = 0;
	virtual bool ready()						 = 0;
	virtual void showBestView()			 = 0;
	virtual void setParm(const QVector <float> &parms, const QString &str) = 0;
};

#endif
