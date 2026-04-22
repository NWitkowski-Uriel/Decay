// gui/mainwindow.cpp
#include "mainwindow.h"

#include <QMessageBox>
#include <QProcess>
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
    resize(900, 550);

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
    rbDirac = new QRadioButton("Dirac (fixed mass)");
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
    testOutput = new QPlainTextEdit();
    testOutput->setReadOnly(true);
    tabWidget->addTab(testOutput, "Decay tests");
    tabWidget->setVisible(true);

    mainLayout->addLayout(selectionLayout);
    mainLayout->addLayout(buttonLayout);
    mainLayout->addWidget(tabWidget);
}

MainWindow::~MainWindow() {}

void MainWindow::onDecayCheckChanged(int state) {
    bool enable = (state == Qt::Checked);
    rbDirac->setEnabled(enable);
    rbBW->setEnabled(enable);
    rbPS->setEnabled(enable);
}

void MainWindow::onRunTestsClicked() {
    QString appPath = QCoreApplication::applicationDirPath();
    QString exec = appPath + "/test_decay_multiplicities";

    if (!QFileInfo::exists(exec)) {
        QMessageBox::critical(this, "Error", "test_decay_multiplicities not found");
        return;
    }

    QProcess proc;
    proc.start(exec, QStringList());
    proc.waitForFinished(-1);

    testOutput->setPlainText(proc.readAllStandardOutput());
}

void MainWindow::onComputeClicked() {
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

    QString execPath = QCoreApplication::applicationDirPath();

    // 1. Compute spectra (ROOT file)
    {
        QProcess proc;
        proc.start(execPath + "/compute_spectra",
                   QStringList() << "--particles=" + particles.join(","));
        proc.waitForFinished(-1);

        if (proc.exitCode() != 0) {
            QMessageBox::critical(this, "Error", "compute_spectra failed");
            return;
        }
    }

    // 2. Export publication plots
    if (cbMt->isChecked() || cbRatio->isChecked()) {
        QProcess proc;
        proc.start(execPath + "/export_publication_plots",
                   QStringList() << "--particles=" + particles.join(","));
        proc.waitForFinished(-1);

        if (proc.exitCode() != 0) {
            QMessageBox::critical(this, "Error", "plot export failed");
            return;
        }
    }

    QMessageBox::information(this, "Done", "Plots saved in output/");
}
