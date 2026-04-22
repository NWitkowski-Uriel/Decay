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
#include <QTextEdit>
#include <QLabel>
#include <QProcess>
#include <QList>
#include <QProgressBar>

struct ProcessStep {
    QString program;
    QStringList args;
    QString description;
};

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void onComputeClicked();
    void onDecayCheckChanged(int state);
    void onRunTestsClicked();
    void onProcessFinished(int exitCode, QProcess::ExitStatus status);
    void onStdout();
    void onStderr();

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
    QTextEdit *logOutput;
    QTextEdit *testOutput;
    QLabel *statusLabel;
    QProgressBar *progressBar;

    QProcess *process;
    QList<ProcessStep> steps;
    int currentStep;
    bool runningTests;

    void startPipeline();
    void startNextStep();
    void appendLog(const QString& text, const QString& color, bool toTests = false);
    void setUiBusy(bool busy);
    void enableModelGroup(bool enable);
};

#endif
