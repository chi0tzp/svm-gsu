#ifndef SVMGSUCANVAS_H
#define SVMGSUCANVAS_H

#include <QMainWindow>

namespace Ui {
class SvmGsuCanvas;
}

class SvmGsuCanvas : public QMainWindow
{
    Q_OBJECT

public:
    explicit SvmGsuCanvas(QWidget *parent = 0);
    ~SvmGsuCanvas();

private:
    Ui::SvmGsuCanvas *ui;
};

#endif // SVMGSUCANVAS_H
