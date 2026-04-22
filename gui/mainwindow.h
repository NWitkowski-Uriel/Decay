#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QCheckBox>
#include <QPushButton>
#include <QTabWidget>
#include <QGroupBox>
#include <QRadioButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPlainTextEdit>

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void onComputeClicked();
    void onDecayCheckChanged(int state);
    void onRunTestsClicked();

private:
    QWidget *centralWidget;
    QVBoxLayout *mainLayout;
    QHBoxLayout *selectionLayout;
    QHBoxLayout *buttonLayout;

    QGroupBox *particleGroup;
    QGroupBox *distGroup;
    QGroupBox *modelGroup;

    QCheckBox *cbProton;
    QCheckBox *cbNeutron;
    QCheckBox *cbPiPlus;
    QCheckBox *cbPiMinus;
    QCheckBox *cbPiZero;

    QCheckBox *cbMt;
    QCheckBox *cbRapidity;
    QCheckBox *cbRatio;

    QCheckBox *cbIncludeDecays;
    QRadioButton *rbDirac;
    QRadioButton *rbBW;
    QRadioButton *rbPS;

    QPushButton *computeButton;
    QPushButton *testButton;
    QTabWidget *tabWidget;
    QPlainTextEdit *testOutput;

    void enableModelGroup(bool enable);
};

#endif // MAINWINDOW_H
