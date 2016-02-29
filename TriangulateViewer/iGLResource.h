#ifndef _I_GL_RESOURCE_H_
#define _I_GL_RESOURCE_H_

class IViewer;

class IGLResource
{
public:
	virtual void initGL(IViewer *viewer) = 0;
	virtual void clearGL()							 = 0;
};

#endif