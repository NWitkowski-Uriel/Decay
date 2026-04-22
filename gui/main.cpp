// gui/main.cpp
#include <QApplication>
#include <TApplication.h>
#include "mainwindow.h"

int main(int argc, char *argv[]) {
    // ROOT application musi być utworzone przed QApplication, aby działały canvasy
    TApplication rootApp("DeltaAnalysis", &argc, argv);

    QApplication app(argc, argv);

    MainWindow w;
    w.show();

    return app.exec();
}