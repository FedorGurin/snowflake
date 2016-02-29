#pragma once

#include "gl_version.h"
#include "qdebug.h"

class TimeMeasurer
{
private:

	Current_OpenGL_Version_Profile *gl;
	int frameCounter;
	double timeBuffer;
	int frameQuantity;
	GLuint64 startTime, stopTime;
	GLuint queryID[2];

public:
	TimeMeasurer(): frameCounter(0), timeBuffer(0), frameQuantity(0), gl(0){}
	TimeMeasurer(Current_OpenGL_Version_Profile *_gl, int _frameQuantity);
	~TimeMeasurer(void);
	void startMeasuring();
	void stopMeasuring();
};

