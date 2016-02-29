#include "mainwindow.h"
#include "scene3d.h"
#include "iViewer.h"

CTriangulateViewer::CTriangulateViewer(QWidget *parent)	: QMainWindow(parent)
{
	setupUi(this);

	scene = new CScene3D();

	viewer = IViewer::create(this, false);
	
	sceneLayout->addWidget(viewer->getWidget());

	viewer->setScene(scene);
}

CTriangulateViewer::~CTriangulateViewer()
{
	viewer->setScene(0);
	delete scene;	
}

void CTriangulateViewer::sendParms()
{

}