#include "mainwindow.h"

#include <QMessageBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QPushButton>
#include <QTabWidget>
#include <QDir>
#include <QCoreApplication>
#include <QFileInfo>

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    setWindowTitle("Delta Analysis – GUI");
    resize(900, 650);

    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    mainLayout = new QVBoxLayout(centralWidget);
    selectionLayout = new QHBoxLayout();
    buttonLayout = new QHBoxLayout();

    particleGroup = new QGroupBox("Particles");
    QVBoxLayout *partLayout = new QVBoxLayout(particleGroup);
    cbProton = new QCheckBox("Proton");
    cbNeutron = new QCheckBox("Neutron");
    cbPiPlus = new QCheckBox("#pi^{+}");
    cbPiMinus = new QCheckBox("#pi^{-}");
    cbPiZero = new QCheckBox("#pi^{0}");
    partLayout->addWidget(cbProton);
    partLayout->addWidget(cbNeutron);
    partLayout->addWidget(cbPiPlus);
    partLayout->addWidget(cbPiMinus);
    partLayout->addWidget(cbPiZero);

    distGroup = new QGroupBox("Distribution type");
    QVBoxLayout *distLayout = new QVBoxLayout(distGroup);
    cbMt = new QCheckBox("Transverse mass (m_t)");
    cbRapidity = new QCheckBox("Rapidity (y)");
    cbRatio = new QCheckBox("Total/Primordial ratio");
    distLayout->addWidget(cbMt);
    distLayout->addWidget(cbRapidity);
    distLayout->addWidget(cbRatio);

    modelGroup = new QGroupBox("Delta decay model");
    QVBoxLayout *modelLayout = new QVBoxLayout(modelGroup);
    cbIncludeDecays = new QCheckBox("Include Delta decays");
    cbIncludeDecays->setChecked(true);
    connect(cbIncludeDecays, &QCheckBox::stateChanged, this, &MainWindow::onDecayCheckChanged);
    rbDirac = new QRadioButton("Dirac");
    rbBW = new QRadioButton("Breit‑Wigner");
    rbPS = new QRadioButton("Phase Shift");
    rbDirac->setChecked(true);
    modelLayout->addWidget(cbIncludeDecays);
    modelLayout->addWidget(rbDirac);
    modelLayout->addWidget(rbBW);
    modelLayout->addWidget(rbPS);

    selectionLayout->addWidget(particleGroup);
    selectionLayout->addWidget(distGroup);
    selectionLayout->addWidget(modelGroup);

    computeButton = new QPushButton("Compute & Plot");
    testButton = new QPushButton("Run decay tests");

    connect(computeButton, &QPushButton::clicked, this, &MainWindow::onComputeClicked);
    connect(testButton, &QPushButton::clicked, this, &MainWindow::onRunTestsClicked);

    buttonLayout->addWidget(computeButton);
    buttonLayout->addWidget(testButton);

    tabWidget = new QTabWidget();
    logOutput = new QTextEdit();
    logOutput->setReadOnly(true);
    testOutput = new QTextEdit();
    testOutput->setReadOnly(true);
    tabWidget->addTab(logOutput, "Logs");
    tabWidget->addTab(testOutput, "Decay tests");

    statusLabel = new QLabel("Idle");
    progressBar = new QProgressBar();
    progressBar->setMinimum(0);
    progressBar->setValue(0);

    mainLayout->addLayout(selectionLayout);
    mainLayout->addLayout(buttonLayout);
    mainLayout->addWidget(tabWidget);
    mainLayout->addWidget(progressBar);
    mainLayout->addWidget(statusLabel);

    process = new QProcess(this);
    connect(process, &QProcess::readyReadStandardOutput, this, &MainWindow::onStdout);
    connect(process, &QProcess::readyReadStandardError, this, &MainWindow::onStderr);
    connect(process, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
            this, &MainWindow::onProcessFinished);

    currentStep = 0;
    runningTests = false;
}

void MainWindow::appendLog(const QString& text, const QString& color, bool toTests) {
    QString html = "<span style='color:" + color + "'>" + text.toHtmlEscaped() + "</span>";
    if (toTests)
        testOutput->append(html);
    else
        logOutput->append(html);
}

void MainWindow::setUiBusy(bool busy) {
    computeButton->setEnabled(!busy);
    testButton->setEnabled(!busy);
}

void MainWindow::onStdout() {
    appendLog(QString::fromLocal8Bit(process->readAllStandardOutput()), "black", runningTests);
}

void MainWindow::onStderr() {
    appendLog(QString::fromLocal8Bit(process->readAllStandardError()), "red", runningTests);
}

void MainWindow::onProcessFinished(int exitCode, QProcess::ExitStatus) {
    progressBar->setValue(currentStep);

    if (exitCode != 0) {
        statusLabel->setText("Error");
        setUiBusy(false);
        return;
    }

    startNextStep();
}

void MainWindow::startPipeline() {
    currentStep = 0;
    progressBar->setMaximum(steps.size());
    progressBar->setValue(0);
    logOutput->clear();
    setUiBusy(true);
    startNextStep();
}

void MainWindow::startNextStep() {
    if (currentStep >= steps.size()) {
        statusLabel->setText("Done");
        setUiBusy(false);
        return;
    }

    auto step = steps[currentStep++];
    statusLabel->setText(step.description);
    appendLog("\n=== " + step.description + " ===", "blue", runningTests);
    process->start(step.program, step.args);
}

void MainWindow::onComputeClicked() {
    runningTests = false;

    QStringList particles;
    if (cbProton->isChecked()) particles << "proton";
    if (cbNeutron->isChecked()) particles << "neutron";
    if (cbPiPlus->isChecked()) particles << "piplus";
    if (cbPiMinus->isChecked()) particles << "piminus";
    if (cbPiZero->isChecked()) particles << "pi0";

    if (particles.isEmpty()) {
        QMessageBox::warning(this, "Selection", "Select particles");
        return;
    }

    QString appPath = QCoreApplication::applicationDirPath();

    steps.clear();

    steps.append({appPath + "/compute_spectra",
                  {"--particles=" + particles.join(",")},
                  "Running compute_spectra"});

    if (cbMt->isChecked() || cbRatio->isChecked()) {
        steps.append({appPath + "/export_publication_plots",
                      {"--particles=" + particles.join(",")},
                      "Generating plots"});
    }

    startPipeline();
}

void MainWindow::onRunTestsClicked() {
    runningTests = true;

    QString exec = QCoreApplication::applicationDirPath() + "/test_decay_multiplicities";
    steps.clear();
    steps.append({exec, {}, "Running decay tests"});
    startPipeline();
}

void MainWindow::onDecayCheckChanged(int state) {
    bool enable = (state == Qt::Checked);
    rbDirac->setEnabled(enable);
    rbBW->setEnabled(enable);
    rbPS->setEnabled(enable);
}
