#include "scene3d.h"
#include "iViewer.h"
#include <QOpenGLShaderProgram>
#include <QFile>
#include "libgpujpeg/gpujpeg.h"
#include "omp.h"
#include "stdio.h"
#include <qthread.h>

extern QString dataFile;
/*
#define PREFIX "1024_1024"
#define FILE_TEST "./input_test_image_1024_1024.jpg"
#define IMAGE_WIDTH 1024
#define IMAGE_HEIGH 1024
#define RESTART_INTERVAL 8
*/

/*#define PREFIX "512_512"
#define FILE_TEST "./input_test_image_512_512.jpg"
#define IMAGE_WIDTH 512
#define IMAGE_HEIGH 512
#define RESTART_INTERVAL 0*/


#define PREFIX "256_256"
#define FILE_TEST "./input_test_image_256_256.jpg"
#define IMAGE_WIDTH 256
#define IMAGE_HEIGH 256
#define RESTART_INTERVAL 8

///////////////////////////////
CScene3D::~CScene3D()
{
	if(viewer)
	{
		viewer->makeCurrent();
		viewer->setScene(0);
		clearGL();
	}
	clear();
}

float distance(QPointF p1, QPointF p2)
{
	return qSqrt((p1.x()-p2.x())*(p1.x()-p2.x())+(p1.y()-p2.y())*(p1.y()-p2.y()));
}

void CScene3D::setParm(const QVector <float> &parms, const QString &str)
{

	update();

}

struct cudaGraphicsResource *cuda_vbo_resource=0;
GLuint texture_vbo_id;
void CScene3D::initGL(IViewer *view)
{
	timeMeasurer = new TimeMeasurer(gl, 10);
	sceneRadius = 700;
	time1.start();

		QFile file("mesaure.txt");
	file.open(QFile::WriteOnly | QIODevice::Text |QIODevice::Append);
	QTextStream out(&file);
	
	if(viewer)
		clearGL();
	viewer = view;
	viewer->setBackgroundColor(Qt::gray);

	prg = new QOpenGLShaderProgram();

	prg->addShaderFromSourceFile(QOpenGLShader::Vertex,   ":/test.vsh.h");
	prg->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/test.fsh.h");

	prg->link();	

	
	//! размеры изображения
	//int image_width = 256;
	//int image_height = 256;

	int image_width = IMAGE_WIDTH;
	int image_height = IMAGE_HEIGH;
	// создаем текстуру для отрисовки
	gl->glGenTextures(1, &tb0);
  gl->glBindTexture(GL_TEXTURE_2D, tb0);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	gl->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16.0f);
  gl->glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image_width,image_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  gl->glBindTexture(GL_TEXTURE_2D, 0);	
	

	//! инициализация устройства
//	if ( gpujpeg_init_device(0, GPUJPEG_OPENGL_INTEROPERABILITY ))
	if ( gpujpeg_init_device(0,GPUJPEG_VERBOSE))
	{
		printf("init_device faild");
	}	
	//! содание декодера
	struct gpujpeg_decoder* decoder = gpujpeg_decoder_create();
	if ( decoder == NULL )
	{
		qDebug("Error in decoder\n");
	}
	// загрузка jpg из файла
	int image_size = 0;
	uint8_t* image = NULL;
	//настройки по умолчанию
	struct gpujpeg_image_parameters param_image;
  gpujpeg_image_set_default_parameters(&param_image);
  param_image.width = image_width;
  param_image.height = image_height;
	

	//настройки по умолчанию
  struct gpujpeg_parameters param;
  gpujpeg_set_default_parameters(&param);
	param.restart_interval = RESTART_INTERVAL;//0 - CPU, restart_interval>0 - GPU (MCU = restart_interval)
  param.interleaved = 1;//чередование компанент
	param.verbose =1;
	param.segment_info=1;

	//! стандартный subSampling
  gpujpeg_parameters_chroma_subsampling(&param);	  
	//! инициализация декодера
	gpujpeg_decoder_init(decoder, &param, &param_image);	
	
	int isOpen = gpujpeg_image_load_from_file(FILE_TEST, &image, &image_size);
	if(isOpen!=0)
	{
		qDebug("Can`t open jpeg file\n");
	}
	
	//! связывание с текстурой
	struct gpujpeg_opengl_texture *texture=NULL;
	texture = gpujpeg_opengl_texture_register(tb0,GPUJPEG_OPENGL_TEXTURE_WRITE);
	// настройка декодера
	struct gpujpeg_decoder_output decoder_output;
	double t1_s = omp_get_wtime();
	for(int i=0;i<1;i++)
	{		
		//! изображение в буфер opengl	
		gpujpeg_decoder_output_set_texture(&decoder_output,texture);
		//gpujpeg_decoder_output_set_cuda_buffer(&decoder_output);
		// декодирование данных
		if ( gpujpeg_decoder_decode(decoder, image, image_size, &decoder_output) != 0 )
		{
			qDebug("Error in decoder_decode\n");
		} 	
	}

	double t1_e = omp_get_wtime();

	double period1=t1_e-t1_s;
	printf("Period=%f/n",period1);
	gpujpeg_opengl_texture_unregister(texture);
	gpujpeg_image_destroy(image);
	gpujpeg_decoder_destroy(decoder);
		
	//QImage img1(FILE_TEST);
	//QImage img0(FILE_TEST);	

	 QThread::sleep(10);
	double t2_s = omp_get_wtime();
	/*for(int i=0;i<1000;i++)
	{	
	
	
	img0 = QGLWidget::convertToGLFormat(img1);
	if(!img0.isNull())
	{
	
		gl->glActiveTexture(GL_TEXTURE0);
		gl->glBindTexture(GL_TEXTURE_2D, tb0);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	gl->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	gl->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16.0f);
  gl->glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image_width,image_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	gl->glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img0.width(), img0.height(), GL_RGBA, GL_UNSIGNED_BYTE, img0.constBits());
  gl->glBindTexture(GL_TEXTURE_2D, 0);
	
	}
	}*/
		double t2_e = omp_get_wtime();
		double period2=t2_e-t2_s;
	printf("Period=%f/n",period2);
	

	out<<PREFIX<<"QImage t="<<period2<<"\n";
	out<<PREFIX<<"CPU+GPU t="<<period1<<"\n";
//	delay(20);
	if(sceneRadius)
	{
		viewer->setSceneRadius(sceneRadius);
	}
	showBestView();	
}

void CScene3D::clearGL()
{
	if (prg)
		delete prg;

	if (vb)
		gl->glDeleteBuffers(1, &vb);
	ib = vb = 0;

	viewer = NULL;
}


void CScene3D::clear() 
{
	sceneRadius = 1;//distance(sceneRect.center(),sceneRect.topLeft());
}

void CScene3D::showBestView()
{
	if (!viewer)
		return;

	viewer->getCamera3D()->showAll();

}

QString checkGlError(const QString &message)
{
	GLenum err = gl->glGetError();

	if(err == GL_NO_ERROR)
		return "";

	QString str;
	switch(err)
	{
	case 0x0506: str = "invalid frambuffer operation "; break;
	case GL_INVALID_ENUM: str = "GL_INVALID_ENUM "; break;
	case GL_INVALID_VALUE: str = "GL_INVALID_VALUE "; break;
	case GL_INVALID_OPERATION: str = "GL_INVALID_OPERATION "; break;
	case GL_OUT_OF_MEMORY: str = "GL_OUT_OF_MEMORY "; break;
	case GL_STACK_OVERFLOW: str = "GL_STACK_OVERFLOW "; break;
	case GL_STACK_UNDERFLOW: str = "GL_STACK_UNDERFLOW "; break;
	default: str = "unknown GL error "; break;
		qDebug() << message << str;
	}
	return str + message;
}

void CScene3D::draw(int)
{
	timeMeasurer->startMeasuring();
	viewer->makeCurrent();

	gl->glEnable(GL_MULTISAMPLE);
	gl->glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE);
	gl->glDisable(GL_STENCIL_TEST);
	gl->glEnable(GL_DEPTH_TEST);
	gl->glDisable(GL_BLEND);
	//gl->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	gl->glDisable(GL_CULL_FACE);

	//gl->glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	gl->glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	Matrix4D viewMatrix = viewer->getCamera3D()->getViewMatrix();

	Matrix4D projectionMatrix = viewer->getCamera3D()->getProjectionMatrix();



	float time = time1.elapsed()/100;

	prg->bind();
	prg->setUniformValue("modelViewMatrix", viewMatrix.toQMatrix4x4());
	prg->setUniformValue("modelViewMatrixInverted", viewMatrix.toQMatrix4x4().inverted());
	prg->setUniformValue("projectionMatrix", projectionMatrix.toQMatrix4x4());
	prg->setUniformValue("timeStamp", time);
	prg->setUniformValue("pictureTex", 0);
	/*gl->glBindBuffer(GL_ARRAY_BUFFER, 0);
	gl->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	gl->glDrawArrays(GL_TRIANGLE_STRIP,0,4);*/
	//checkGlError("draw");

	gl->glBindVertexArray(0);

	gl->glActiveTexture(GL_TEXTURE0);
	gl->glBindTexture(GL_TEXTURE_2D, tb0);

	gl->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

	timeMeasurer->stopMeasuring();

}
