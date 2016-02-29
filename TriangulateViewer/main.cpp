#include <QFileInfo>
#include <QDir>
#include "mainwindow.h"

QString dataFile = "triangulation.dat";

int main(int argc, char *argv[])
{
	QDir subdir = QFileInfo(argv[0]).absoluteDir();

	if(argc == 2)
		dataFile = subdir.absoluteFilePath(QString(argv[1]).simplified());
		
	QApplication a(argc, argv);

	CTriangulateViewer w;
	w.show();

	return a.exec();
}
