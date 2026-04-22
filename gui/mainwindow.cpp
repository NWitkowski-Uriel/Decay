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
    resize(800, 450);

    centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    mainLayout = new QVBoxLayout(centralWidget);
    selectionLayout = new QHBoxLayout();

    // --- Grupa cząstek ---
    particleGroup = new QGroupBox("Particles");
    QVBoxLayout *partLayout = new QVBoxLayout(particleGroup);
    cbProton   = new QCheckBox("Proton");
    cbNeutron  = new QCheckBox("Neutron");
    cbPiPlus   = new QCheckBox("#pi^{+}");
    cbPiMinus  = new QCheckBox("#pi^{-}");
    cbPiZero   = new QCheckBox("#pi^{0}");
    partLayout->addWidget(cbProton);
    partLayout->addWidget(cbNeutron);
    partLayout->addWidget(cbPiPlus);
    partLayout->addWidget(cbPiMinus);
    partLayout->addWidget(cbPiZero);

    // --- Grupa typów rozkładów ---
    distGroup = new QGroupBox("Distribution type");
    QVBoxLayout *distLayout = new QVBoxLayout(distGroup);
    cbMt       = new QCheckBox("Transverse mass (m_t)");
    cbRapidity = new QCheckBox("Rapidity (y)");
    cbRatio    = new QCheckBox("Total/Primordial ratio");
    distLayout->addWidget(cbMt);
    distLayout->addWidget(cbRapidity);
    distLayout->addWidget(cbRatio);

    // --- Grupa modeli (rozpady Δ) ---
    modelGroup = new QGroupBox("Delta decay model");
    QVBoxLayout *modelLayout = new QVBoxLayout(modelGroup);
    cbIncludeDecays = new QCheckBox("Include Delta decays");
    cbIncludeDecays->setChecked(true);
    connect(cbIncludeDecays, &QCheckBox::stateChanged, this, &MainWindow::onDecayCheckChanged);
    rbDirac = new QRadioButton("Dirac (fixed mass)");
    rbBW    = new QRadioButton("Breit‑Wigner (Γ=120 MeV)");
    rbPS    = new QRadioButton("Phase Shift");
    rbDirac->setChecked(true);
    modelLayout->addWidget(cbIncludeDecays);
    modelLayout->addWidget(rbDirac);
    modelLayout->addWidget(rbBW);
    modelLayout->addWidget(rbPS);

    selectionLayout->addWidget(particleGroup);
    selectionLayout->addWidget(distGroup);
    selectionLayout->addWidget(modelGroup);
    selectionLayout->addStretch();

    computeButton = new QPushButton("Compute & Plot");
    connect(computeButton, &QPushButton::clicked, this, &MainWindow::onComputeClicked);

    tabWidget = new QTabWidget();
    tabWidget->setVisible(false);

    mainLayout->addLayout(selectionLayout);
    mainLayout->addWidget(computeButton);
    mainLayout->addWidget(tabWidget);
}

MainWindow::~MainWindow() {}

void MainWindow::onDecayCheckChanged(int state) {
    bool enable = (state == Qt::Checked);
    rbDirac->setEnabled(enable);
    rbBW->setEnabled(enable);
    rbPS->setEnabled(enable);
}

void MainWindow::onComputeClicked() {
    // Zbierz wybory
    QStringList particles;
    if (cbProton->isChecked())   particles << "proton";
    if (cbNeutron->isChecked())  particles << "neutron";
    if (cbPiPlus->isChecked())   particles << "piplus";
    if (cbPiMinus->isChecked())  particles << "piminus";
    if (cbPiZero->isChecked())   particles << "pi0";

    QStringList distributions;
    if (cbMt->isChecked())       distributions << "mt";
    if (cbRapidity->isChecked()) distributions << "rapidity";
    if (cbRatio->isChecked())    distributions << "ratio";

    QString model;
    if (!cbIncludeDecays->isChecked()) {
        model = "primordial";
    } else {
        if (rbDirac->isChecked()) model = "dirac";
        else if (rbBW->isChecked()) model = "bw";
        else model = "ps";
    }

    if (particles.isEmpty() || distributions.isEmpty()) {
        QMessageBox::warning(this, "No selection",
                             "Please select at least one particle and one distribution type.");
        return;
    }

    // Utwórz katalog wyjściowy
    QDir dir;
    if (!dir.exists("output")) {
        if (!dir.mkdir("output")) {
            QMessageBox::critical(this, "Error", "Cannot create output directory.");
            return;
        }
    }

    QString appPath = QCoreApplication::applicationDirPath();
    QString computePrimordial = appPath + "/compute_primordial";
    QString computeDeltaDecays = appPath + "/compute_delta_decays";
    QString computeTotal = appPath + "/compute_total";
    QString plotMt = appPath + "/plot_mt";
    QString plotRapidity = appPath + "/plot_rapidity";

    // Sprawdź, czy pliki wykonywalne istnieją
    if (!QFileInfo::exists(computePrimordial)) {
        QMessageBox::critical(this, "Error", "compute_primordial not found in " + appPath);
        return;
    }
    if (cbIncludeDecays->isChecked() && !QFileInfo::exists(computeDeltaDecays)) {
        QMessageBox::critical(this, "Error", "compute_delta_decays not found in " + appPath);
        return;
    }
    if (!QFileInfo::exists(computeTotal)) {
        QMessageBox::critical(this, "Error", "compute_total not found in " + appPath);
        return;
    }
    if ((distributions.contains("mt") || distributions.contains("rapidity")) && !QFileInfo::exists(plotMt)) {
        QMessageBox::critical(this, "Error", "plot_mt not found in " + appPath);
        return;
    }
    if ((distributions.contains("rapidity") || distributions.contains("ratio")) && !QFileInfo::exists(plotRapidity)) {
        QMessageBox::critical(this, "Error", "plot_rapidity not found in " + appPath);
        return;
    }

    // -----------------------------------------------------------------
    // 1. Uruchom compute_primordial (bez timeoutu)
    // -----------------------------------------------------------------
    QProcess p1;
    p1.start(computePrimordial, QStringList());
    if (!p1.waitForStarted()) {
        QMessageBox::critical(this, "Error", "Could not start compute_primordial: " + p1.errorString());
        return;
    }
    // Czekaj w nieskończoność, aż proces się zakończy
    if (!p1.waitForFinished()) {
        QMessageBox::critical(this, "Error", "compute_primordial did not finish.");
        return;
    }
    if (p1.exitCode() != 0) {
        QMessageBox::critical(this, "Error",
                              "compute_primordial failed.\n" + p1.readAllStandardError());
        return;
    }

    // Sprawdź, czy plik wynikowy został utworzony
    if (!QFile::exists("primordial.root")) {
        QMessageBox::critical(this, "Error", "compute_primordial did not create primordial.root");
        return;
    }

    // -----------------------------------------------------------------
    // 2. Uruchom compute_delta_decays (tylko jeśli włączone rozpady)
    // -----------------------------------------------------------------
    if (cbIncludeDecays->isChecked()) {
        QStringList deltaArgs;
        deltaArgs << "--particles=" + particles.join(",");
        deltaArgs << "--model=" + model;

        QProcess p2;
        p2.start(computeDeltaDecays, deltaArgs);
        if (!p2.waitForStarted()) {
            QMessageBox::critical(this, "Error", "Could not start compute_delta_decays: " + p2.errorString());
            return;
        }
        if (!p2.waitForFinished()) {
            QMessageBox::critical(this, "Error", "compute_delta_decays did not finish.");
            return;
        }
        if (p2.exitCode() != 0) {
            QMessageBox::critical(this, "Error",
                                  "compute_delta_decays failed.\n" + p2.readAllStandardError());
            return;
        }
        // Sprawdź, czy plik wynikowy istnieje (jeśli nie, to znaczy, że nie było żadnych rozpadów – ale to normalne)
    }

    // -----------------------------------------------------------------
    // 3. Uruchom compute_total
    // -----------------------------------------------------------------
    QProcess p3;
    p3.start(computeTotal, QStringList());
    if (!p3.waitForStarted()) {
        QMessageBox::critical(this, "Error", "Could not start compute_total: " + p3.errorString());
        return;
    }
    if (!p3.waitForFinished()) {
        QMessageBox::critical(this, "Error", "compute_total did not finish.");
        return;
    }
    if (p3.exitCode() != 0) {
        QMessageBox::critical(this, "Error",
                              "compute_total failed.\n" + p3.readAllStandardError());
        return;
    }
    if (!QFile::exists("total.root")) {
        QMessageBox::critical(this, "Error", "compute_total did not create total.root");
        return;
    }

    // -----------------------------------------------------------------
    // 4. Uruchom programy rysujące (bez oczekiwania)
    // -----------------------------------------------------------------
    QStringList plotArgs;
    plotArgs << "--distributions=" + distributions.join(",");
    plotArgs << "--particles=" + particles.join(",");

    if (distributions.contains("mt")) {
        QProcess::startDetached(plotMt, plotArgs);
    }
    if (distributions.contains("rapidity")) {
        QProcess::startDetached(plotRapidity, plotArgs);
    }

    QMessageBox::information(this, "Done",
                             "Calculations finished.\n"
                             "Results are in the 'output' folder:\n"
                             " - Data: output/total.root\n"
                             " - Plots: PDF/PNG files");
}