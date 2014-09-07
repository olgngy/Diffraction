#include "mainwindow.h"
#include "laserenv.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    // test comment
    QApplication a(argc, argv);
    //MainWindow w;
    LaserEnv w;
    w.show();
    return a.exec();
}
