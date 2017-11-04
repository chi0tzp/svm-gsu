#include "svmgsucanvas.h"
#include "ui_svmgsucanvas.h"

SvmGsuCanvas::SvmGsuCanvas(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::SvmGsuCanvas)
{
    ui->setupUi(this);
}

SvmGsuCanvas::~SvmGsuCanvas()
{
    delete ui;
}
