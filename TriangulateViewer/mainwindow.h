#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "ui_mainwindow.h"
class IScene;
class IViewer;

class CTriangulateViewer : public QMainWindow, private Ui::UITriangulateViewer
{
	Q_OBJECT

	IViewer* viewer;
	IScene* scene;

public:
	CTriangulateViewer(QWidget *parent = 0);
	~CTriangulateViewer();
private slots:
	void sendParms();
};

#endif // MAINWINDOW_H
