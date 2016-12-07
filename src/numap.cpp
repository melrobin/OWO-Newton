#include "numap.h"
#include "numap_ui.h"
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>
#include <QtCore/QDir>
NuMap::NuMap(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::NuMap)
{
    ui->setupUi(this);
}

NuMap::~NuMap()
{
    delete ui;
}

void NuMap::on_pushButton_clicked()
{
   // QMessageBox("We're going to start training here!")
}

void NuMap::on_pushButton_2_clicked()
{
   // QMessageBox("We're going to start validing here!")
}

void NuMap::on_actionExit_triggered()
{

}

void NuMap::on_actionSet_Training_File_triggered()
{
    QDir directory;
    QString path = QFileDialog::getExistingDirectory (this, tr("Directory"), directory.path());
    if ( path.isNull() == false )
    {
        directory.setPath(path);
    }
}
