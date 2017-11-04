#include "svmgsucanvas.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    SvmGsuCanvas w;
    w.show();

    return a.exec();
}
