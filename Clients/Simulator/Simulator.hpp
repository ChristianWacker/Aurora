//--- Aurora/Clients/Simulator/Simulator.hpp -----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_SIMULATOR_SIMULATOR_HPP
#define AURORA_CLIENTS_SIMULATOR_SIMULATOR_HPP

#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/TemSimulation.hpp"
#include "AuroraLib/MultisliceParameters.hpp"
#include "AuroraLib/Sample.hpp"
#include "qwt_plot_curve.h"
#include "ui_SimulatorMain.h"

#include <mutex>
#include <thread>

namespace Aurora
{

class Simulator : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT

public:
  Simulator(QMainWindow* parent = nullptr);
  ~Simulator() {}

  void processWave(size_t sliceIndex, size_t numSlices, const CBuffer2D& wave);

signals:
  void waveCalculated();
  void progressUpdated(size_t part, size_t total);

private:
  CtfCoherent mCtf;
  MultisliceParameters mParams;
  SamplePtr mSample;
  Image mProjPotentialImage;
  QPixmap mProjPotentialPixmap;
  CBuffer2D mExitWave;
  CBuffer2D mExitWaveFft;

  QString mConfigPath;
  QString mExitWavePath;
  QString mProjPotentialPath;
  QString mSamplePath;

  QwtPlotCurve mCurveReal;
  QwtPlotCurve mCurveImag;

  void displayWave();

  void drawHistogram(const Histogram& histogram);

  void savePotential(size_t sliceIndex, size_t numSlices,
                     const CBuffer2D& potential);
  void loadSample();

  void updateInterface();
  void updateProgressBar(size_t part, size_t total);

  bool mInterfaceUpdatePending;

  void doSimulate();
  void startSimulation();
  void endSimulation();
  void loadConfigFile();
  void loadSampleClicked();
  void changeOptics();
  void changeSimulation();
  void saveConfigFile();
  void showProjPotentialContextMenu(const QPoint& position);
  void saveProjPotentialImage();
  void saveProjPotentialBuffer();
  void saveProjPotentialBuffer2();
  void extractBox();

  bool mSimulationRunning = false;
  std::mutex mMutex;
  std::thread mWorkThread;
};

} // namespace Aurora

#endif
