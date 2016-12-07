/********************************************************************************
** Form generated from reading UI file 'numap.ui'
**
** Created by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_NUMAP_H
#define UI_NUMAP_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLCDNumber>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_NuMap
{
public:
    QAction *actionNew;
    QAction *actionSave;
    QAction *actionExit;
    QAction *actionIndex;
    QAction *actionAbout;
    QAction *actionFile_Statistics;
    QAction *actionCount_Patterns;
    QAction *actionSet_Training_File;
    QAction *actionSet_Validation_File;
    QWidget *centralWidget;
    QComboBox *comboBox;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QCheckBox *checkBox;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QSpinBox *spinBox;
    QSpinBox *spinBox_2;
    QSpinBox *spinBox_3;
    QSpinBox *spinBox_4;
    QLCDNumber *lcdNumber;
    QLabel *label_6;
    QLabel *label_7;
    QLCDNumber *lcdNumber_2;
    QLabel *label_8;
    QProgressBar *progressBar;
    QLabel *label_9;
    QLabel *label_10;
    QMenuBar *menuBar;
    QMenu *menu_File;
    QMenu *menuHelp;
    QMenu *menuUtility;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *NuMap)
    {
        if (NuMap->objectName().isEmpty())
            NuMap->setObjectName(QString::fromUtf8("NuMap"));
        NuMap->resize(612, 323);
        actionNew = new QAction(NuMap);
        actionNew->setObjectName(QString::fromUtf8("actionNew"));
        actionSave = new QAction(NuMap);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        actionExit = new QAction(NuMap);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionIndex = new QAction(NuMap);
        actionIndex->setObjectName(QString::fromUtf8("actionIndex"));
        actionAbout = new QAction(NuMap);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        actionFile_Statistics = new QAction(NuMap);
        actionFile_Statistics->setObjectName(QString::fromUtf8("actionFile_Statistics"));
        actionCount_Patterns = new QAction(NuMap);
        actionCount_Patterns->setObjectName(QString::fromUtf8("actionCount_Patterns"));
        actionSet_Training_File = new QAction(NuMap);
        actionSet_Training_File->setObjectName(QString::fromUtf8("actionSet_Training_File"));
        actionSet_Validation_File = new QAction(NuMap);
        actionSet_Validation_File->setObjectName(QString::fromUtf8("actionSet_Validation_File"));
        centralWidget = new QWidget(NuMap);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        comboBox = new QComboBox(centralWidget);
        comboBox->setObjectName(QString::fromUtf8("comboBox"));
        comboBox->setGeometry(QRect(150, 140, 121, 25));
        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(400, 240, 84, 26));
        pushButton_2 = new QPushButton(centralWidget);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        pushButton_2->setGeometry(QRect(500, 240, 84, 26));
        checkBox = new QCheckBox(centralWidget);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));
        checkBox->setGeometry(QRect(10, 210, 87, 21));
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 140, 111, 31));
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 20, 59, 15));
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(10, 50, 59, 15));
        label_4 = new QLabel(centralWidget);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(10, 80, 91, 16));
        label_5 = new QLabel(centralWidget);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(10, 110, 71, 16));
        spinBox = new QSpinBox(centralWidget);
        spinBox->setObjectName(QString::fromUtf8("spinBox"));
        spinBox->setGeometry(QRect(150, 10, 55, 25));
        spinBox_2 = new QSpinBox(centralWidget);
        spinBox_2->setObjectName(QString::fromUtf8("spinBox_2"));
        spinBox_2->setGeometry(QRect(150, 40, 55, 25));
        spinBox_3 = new QSpinBox(centralWidget);
        spinBox_3->setObjectName(QString::fromUtf8("spinBox_3"));
        spinBox_3->setGeometry(QRect(150, 70, 55, 25));
        spinBox_4 = new QSpinBox(centralWidget);
        spinBox_4->setObjectName(QString::fromUtf8("spinBox_4"));
        spinBox_4->setGeometry(QRect(150, 100, 55, 25));
        lcdNumber = new QLCDNumber(centralWidget);
        lcdNumber->setObjectName(QString::fromUtf8("lcdNumber"));
        lcdNumber->setGeometry(QRect(520, 170, 64, 23));
        label_6 = new QLabel(centralWidget);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(378, 170, 91, 20));
        label_7 = new QLabel(centralWidget);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(380, 210, 111, 16));
        lcdNumber_2 = new QLCDNumber(centralWidget);
        lcdNumber_2->setObjectName(QString::fromUtf8("lcdNumber_2"));
        lcdNumber_2->setGeometry(QRect(520, 200, 64, 23));
        label_8 = new QLabel(centralWidget);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(10, 190, 81, 16));
        progressBar = new QProgressBar(centralWidget);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(160, 240, 118, 23));
        progressBar->setValue(24);
        label_9 = new QLabel(centralWidget);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(10, 240, 111, 21));
        label_10 = new QLabel(centralWidget);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setGeometry(QRect(150, 190, 59, 15));
        NuMap->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(NuMap);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 612, 23));
        menu_File = new QMenu(menuBar);
        menu_File->setObjectName(QString::fromUtf8("menu_File"));
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        menuUtility = new QMenu(menuBar);
        menuUtility->setObjectName(QString::fromUtf8("menuUtility"));
        NuMap->setMenuBar(menuBar);
        mainToolBar = new QToolBar(NuMap);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        NuMap->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(NuMap);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        NuMap->setStatusBar(statusBar);

        menuBar->addAction(menu_File->menuAction());
        menuBar->addAction(menuUtility->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menu_File->addAction(actionNew);
        menu_File->addAction(actionSave);
        menu_File->addAction(actionExit);
        menu_File->addAction(actionSet_Training_File);
        menu_File->addAction(actionSet_Validation_File);
        menuHelp->addAction(actionIndex);
        menuHelp->addAction(actionAbout);
        menuUtility->addAction(actionFile_Statistics);
        menuUtility->addAction(actionCount_Patterns);

        retranslateUi(NuMap);

        QMetaObject::connectSlotsByName(NuMap);
    } // setupUi

    void retranslateUi(QMainWindow *NuMap)
    {
        NuMap->setWindowTitle(QApplication::translate("NuMap", "NuMap", 0, QApplication::UnicodeUTF8));
        actionNew->setText(QApplication::translate("NuMap", "New", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("NuMap", "Save", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("NuMap", "Exit", 0, QApplication::UnicodeUTF8));
        actionIndex->setText(QApplication::translate("NuMap", "Index", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("NuMap", "About", 0, QApplication::UnicodeUTF8));
        actionFile_Statistics->setText(QApplication::translate("NuMap", "File Statistics", 0, QApplication::UnicodeUTF8));
        actionCount_Patterns->setText(QApplication::translate("NuMap", "Count Patterns", 0, QApplication::UnicodeUTF8));
        actionSet_Training_File->setText(QApplication::translate("NuMap", "Set Training File...", 0, QApplication::UnicodeUTF8));
        actionSet_Validation_File->setText(QApplication::translate("NuMap", "Set Validation File...", 0, QApplication::UnicodeUTF8));
        comboBox->clear();
        comboBox->insertItems(0, QStringList()
         << QApplication::translate("NuMap", "BP", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "OWO-BP", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "OIG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "MOLF", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "LM", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "CG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "OWO-Newton", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "OIT", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("NuMap", "ONT", 0, QApplication::UnicodeUTF8)
        );
        pushButton->setText(QApplication::translate("NuMap", "Train", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("NuMap", "Test", 0, QApplication::UnicodeUTF8));
        checkBox->setText(QApplication::translate("NuMap", "Validate", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("NuMap", "Training Method", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("NuMap", "Inputs", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("NuMap", "Outputs", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("NuMap", "Hidden Units", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("NuMap", "Iterations", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("NuMap", "Training Error", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("NuMap", "Validation Error", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("NuMap", "Training File", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("NuMap", "Training progress", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("NuMap", "TextLabel", 0, QApplication::UnicodeUTF8));
        menu_File->setTitle(QApplication::translate("NuMap", "&File", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("NuMap", "Help", 0, QApplication::UnicodeUTF8));
        menuUtility->setTitle(QApplication::translate("NuMap", "Utility", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class NuMap: public Ui_NuMap {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_NUMAP_H
