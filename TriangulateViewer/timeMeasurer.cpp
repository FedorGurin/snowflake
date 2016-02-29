#include "timeMeasurer.h"


TimeMeasurer::TimeMeasurer(Current_OpenGL_Version_Profile *_gl,  int _frameQuantity)
{
	gl=_gl;
	frameQuantity=_frameQuantity;
	frameCounter=0;
	timeBuffer=0;
}


TimeMeasurer::~TimeMeasurer(void)
{

}

void TimeMeasurer::startMeasuring()
{
	// generate two queries
	gl->glGenQueries(2, queryID);
	// issue the first query
	// Records the time only after all previous
	// commands have been completed
	gl->glQueryCounter(queryID[0], GL_TIMESTAMP);

	
}
void TimeMeasurer::stopMeasuring()
{
	// issue the second query
	// records the time when the sequence of OpenGL
	// commands has been fully executed
	gl->glQueryCounter(queryID[1], GL_TIMESTAMP);
	// wait until the results are available
	GLint stopTimerAvailable = 0;
	while (!stopTimerAvailable) 
	{
		gl->glGetQueryObjectiv(queryID[1], GL_QUERY_RESULT_AVAILABLE, &stopTimerAvailable);
	}
 
	// get query results
	gl->glGetQueryObjectui64v(queryID[0], GL_QUERY_RESULT, &startTime);
	gl->glGetQueryObjectui64v(queryID[1], GL_QUERY_RESULT, &stopTime);
 /*
	qDebug()<<"Time spent on the GPU:" << ((stopTime - startTime) / 1000000.0) << "ms\n";
	*/
	timeBuffer+=((stopTime - startTime) / 1000000.0);
	frameCounter++;
	if(frameCounter==frameQuantity)
	{
		qDebug()<<"Average time spent on the GPU per "<< frameQuantity<< " frames:" << (timeBuffer/(frameQuantity)) << "ms\n";
		qDebug()<<"FPS:" << ((1000/(timeBuffer/(frameQuantity)))) << "frames\n";
		frameCounter=0;
		timeBuffer=0;
	}
}