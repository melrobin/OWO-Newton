#ifndef NUMAP_H
#define NUMAP_H

#include <QtGui/QMainWindow>

namespace Ui {
class NuMap;
}

class NuMap : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit NuMap(QWidget *parent = 0);
    ~NuMap();
    
private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_actionExit_triggered();

    void on_actionSet_Training_File_triggered();

private:
    Ui::NuMap *ui;
};

#endif // NUMAP_H
