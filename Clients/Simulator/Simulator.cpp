//--- Aurora/Clients/Simulator/Simulator.cpp -----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Clients/Simulator/Simulator.hpp"

#include "AuroraLib/Bandlimit.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Histogram.hpp"
#include "AuroraLib/MultisliceOptions.hpp"
#include "AuroraLib/PhaseShift.hpp"
#include "AuroraLib/Potential.hpp"
#include "AuroraLib/Propagator.hpp"
#include "AuroraLib/RBuffer2D.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"
#include "QtSupport/QtSupport.hpp"

#include <fstream>
#include <iostream>
#include <QFileDialog>
#include <qwt_scale_draw.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_grid.h>

namespace Aurora
{

CL::Option<std::string> propagatorName("propagator", "", "Classical-FT");
CL::Option<bool> showProjPotential("showProjPotential", "", true);

class ResolutionScaleDraw : public QwtScaleDraw
{
public:
  QwtText label(double v) const override
  {
    return QwtText(QString::number(10.0 / v, 'f', 2));
  }
};

Simulator::Simulator(QMainWindow* parent) :
  QMainWindow(parent),
  mCtf(1.0), // dummy energy value
  mSample(std::make_shared<Sample>()),
  mProjPotentialPixmap(256, 256),
  mConfigPath(QDir::homePath()),
  mExitWavePath(QDir::homePath()),
  mProjPotentialPath(QDir::homePath()),
  mSamplePath(QDir::homePath()),
  mInterfaceUpdatePending(true)
{
  setupUi(this);

  // abbreviations
  void (QSpinBox::*intChanged)(int) = &QSpinBox::valueChanged;
  void (QDoubleSpinBox::*doubleChanged)(double) = &QDoubleSpinBox::valueChanged;
  auto buttonPushed = &QAbstractButton::clicked;
  auto textChanged = &QComboBox::currentTextChanged;

  connect(energySpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(widthSpinBox, intChanged, this, &Simulator::changeSimulation);
  connect(heightSpinBox, intChanged, this, &Simulator::changeSimulation);
  connect(deltaZSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(defocusSpinBox, doubleChanged, this, &Simulator::changeOptics);
  connect(c3SpinBox, doubleChanged, this, &Simulator::changeOptics);
  connect(c5SpinBox, doubleChanged, this, &Simulator::changeOptics);
  connect(apertureSpinBox, doubleChanged, this, &Simulator::changeOptics);
  connect(loadSampleButton, buttonPushed, this, &Simulator::loadSampleClicked);
  connect(minimumXSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(minimumYSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(minimumZSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(maximumXSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(maximumYSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(maximumZSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(simulateButton, buttonPushed, this, &Simulator::startSimulation);
  connect(phaseShiftComboBox, textChanged, this, &Simulator::changeSimulation);
  connect(potentialComboBox, textChanged, this, &Simulator::changeSimulation);
  connect(propagatorComboBox, textChanged, this, &Simulator::changeSimulation);
  connect(formFactorComboBox, textChanged, this, &Simulator::changeSimulation);
  connect(repetitionsSpinBox, intChanged, this, &Simulator::changeSimulation);
  connect(opticRadioButton, buttonPushed, this, &Simulator::changeOptics);
  connect(exitWaveRadioButton, buttonPushed, this, &Simulator::changeOptics);
  connect(exitWaveFftRadioButton, buttonPushed, this, &Simulator::changeOptics);
  connect(loadConfigButton, buttonPushed, this, &Simulator::loadConfigFile);
  connect(saveConfigButton, buttonPushed, this, &Simulator::saveConfigFile);
  connect(bandlimitComboBox, textChanged, this, &Simulator::changeSimulation);
  connect(bFactorSpinBox, doubleChanged, this, &Simulator::changeSimulation);
  connect(extractBoxButton, buttonPushed, this, &Simulator::extractBox);
  connect(cutoffSpinBox, doubleChanged, this, &Simulator::changeSimulation);

  connect(this, &Simulator::progressUpdated,
          this, &Simulator::updateProgressBar, Qt::QueuedConnection);
  connect(this, &Simulator::waveCalculated,
          this, &Simulator::displayWave, Qt::QueuedConnection);

  // add a context menu to the projected potential
  connect(projPotentialFrame, &QFrame::customContextMenuRequested,
          this, &Simulator::showProjPotentialContextMenu);
  connect(projPotentialLabel, &QLabel::customContextMenuRequested,
          this, &Simulator::showProjPotentialContextMenu);

  auto setToFocus = [=]
  {
    defocusSpinBox->setValue(0.0);
  };
  connect(focusButton, buttonPushed, setToFocus);

  auto setToScherzerDefocus = [=]
  {
    Real scherzerDefocus = mCtf.scherzerDefocus();
    defocusSpinBox->setValue(scherzerDefocus * R(1e-1));
  };
  connect(scherzerButton, buttonPushed, setToScherzerDefocus);

  // setting up the CTF plot
  // make the plot background white including the label and suppress the border
  ctfPlot->setAutoFillBackground(true);
  ctfPlot->setPalette(Qt::white);
  ctfPlot->setCanvasBackground(Qt::white);
  auto canvas = static_cast<QwtPlotCanvas*>(ctfPlot->canvas());
  canvas->setFrameStyle(QFrame::NoFrame);

  // the labels of the axes
  QwtText xTitle("Spatial Frequency [1 / nm]");
  xTitle.setFont(QFont("Times", 12));
  ctfPlot->setAxisTitle(QwtPlot::xBottom, xTitle);

  QwtText x2Title("Resolution [Angstrom]");
  ctfPlot->setAxisScaleDraw(QwtPlot::xTop, new ResolutionScaleDraw);
  x2Title.setFont(QFont("Times", 12));
  ctfPlot->setAxisTitle(QwtPlot::xTop, x2Title);
  ctfPlot->enableAxis(QwtPlot::xTop);

  QwtText yTitle("CTF");
  yTitle.setFont(QFont("Times", 12));
  ctfPlot->setAxisTitle(QwtPlot::yLeft, yTitle);
  ctfPlot->setAxisScale(QwtPlot::yLeft, -1.05, 1.05);

  // add an gray dotted grid
  QwtPlotGrid* grid = new QwtPlotGrid();
  grid->setPen(QPen(Qt::gray, 0, Qt::DotLine));
  grid->attach(ctfPlot);

  mCurveReal.setPen(QPen(Qt::red));
  mCurveImag.setPen(QPen(Qt::blue));

  mCurveReal.attach(ctfPlot);
  mCurveImag.attach(ctfPlot);

  for (auto i : PhaseShift::namesList())
    phaseShiftComboBox->addItem(QString::fromStdString(i));

  for (auto i : Potential::namesList())
    potentialComboBox->addItem(QString::fromStdString(i));

  for (auto i : Propagator::namesList())
    propagatorComboBox->addItem(QString::fromStdString(i));

  for (auto i : FormFactorParametrization::namesList())
    formFactorComboBox->addItem(QString::fromStdString(i));

  for (auto i : bandlimit.valueList())
    bandlimitComboBox->addItem(QString::fromStdString(i));

  mInterfaceUpdatePending = false;

  updateInterface();
  loadSample();
  changeSimulation();
}

void Simulator::processWave(size_t sliceIndex, size_t numSlices,
                            const CBuffer2D& wave)
{
  // This function is executed in the context of the worker thread. Therefore,
  // we cannot call any GUI functions directly. Instead we send signals.

  // Update the progress bar
  emit progressUpdated(sliceIndex, numSlices);

  // we are only interested in the wave after the last slice, the exit wave
  if (numSlices - 1 != sliceIndex)
    return;

  {
    std::unique_lock<std::mutex> lock(mMutex);
    mExitWave = wave.clone();
  }

  emit waveCalculated();
}

/// This function runs in an extra so the GUI won't block.
void runSimulation(Simulator* simulator, const std::string& propagatorName,
                   MultisliceParameters params, const SamplePtr& sample,
                   Real bFactor)
{
  try
  {
    for (auto& atom : sample->atoms())
      atom.setBFactor(bFactor);

//  mParams.setDebugDir(homeDirectory() + "/Desktop/temp/");
//  mParams.setDebug(Debug::phaseShifts);

    TemSimulation tem(params, numRepetitions, sample, propagatorName);

    using namespace std::placeholders;
    tem.onWave = std::bind(&Simulator::processWave, simulator, _1, _2, _3);

    tem.run();
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}

void Simulator::startSimulation()
{
  {
    std::unique_lock<std::mutex> lock(mMutex);

    if (mSimulationRunning)
      return;
    mSimulationRunning = true;

    // deactivate the simulate button
    simulateButton->setEnabled(false);
  }

  updateProgressBar(0, 1);

  std::thread thread(runSimulation, this,
                     propagatorComboBox->currentText().toStdString(),
                     mParams, std::make_shared<Sample>(*mSample),
                     bFactorSpinBox->value());
  thread.detach();
}

void Simulator::showProjPotentialContextMenu(const QPoint& position)
{
  QPoint globalPosition =
    static_cast<QWidget*>(QObject::sender())->mapToGlobal(position);

  QMenu menu;

  QAction saveImage("Save Image...", this);
  menu.addAction(&saveImage);
  connect(&saveImage, &QAction::triggered,
          this, &Simulator::saveProjPotentialImage);

  QAction saveBuffer("Save Buffer...", this);
  menu.addAction(&saveBuffer);
  connect(&saveBuffer, &QAction::triggered,
          this, &Simulator::saveProjPotentialBuffer);

  QAction saveBuffer2("Save Buffer2...", this);
  menu.addAction(&saveBuffer2);
  connect(&saveBuffer2, &QAction::triggered,
          this, &Simulator::saveProjPotentialBuffer2);

  menu.exec(globalPosition);
}

void Simulator::saveProjPotentialImage()
{
  if (!mProjPotentialImage.dataPtr())
    return;

  try
  {
    const QString caption = tr("Save Projected Potential Image");
    const QString filter  = tr("Image (*.tif *.tiff)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "tiff", mProjPotentialPath);

    if (filename.empty())
      return;

    mProjPotentialImage.save(filename);
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void Simulator::saveProjPotentialBuffer()
{
  if (!mSample)
    return;

  try
  {
    const QString caption = tr("Save Projected Potential Buffer");
    const QString filter  = tr("Image (*.tif *.tiff)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "tiff", mProjPotentialPath);

    if (filename.empty())
      return;

    auto rBuffer = mSample->projectedPotential(mParams.cols(),
                                               mParams.rows(),
                                               mParams.minBox(),
                                               mParams.maxBox(),
                                               5.0,
                                               mParams.formFactorName());

    CBuffer2D cBuffer(mParams.cols(), mParams.rows());
    cBuffer.assign(rBuffer);

    cBuffer.save(filename, mParams, "");
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void Simulator::saveProjPotentialBuffer2()
{
  if (!mSample)
    return;

  try
  {
    const QString caption = tr("Save Projected Potential Buffer2");
    const QString filter  = tr("Image (*.tif *.tiff)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "tiff", mProjPotentialPath);

    if (filename.empty())
      return;

    auto potential = Potential::create(mParams.potentialName(), mParams,
                                       mSample);
    CBuffer2D cBuffer(mParams.cols(), mParams.rows(), Complex());

    #pragma omp parallel for
    for (int64_t y = 0; y < mParams.rows(); ++y)
      for (int64_t x = 0; x < mParams.cols(); ++x)
        for (size_t z = 0; z < mParams.numSlices(); ++z)
          cBuffer.pixel(x, y) += potential->potential(z).pixel(x, y);

    cBuffer *= mParams.deltaZ();

    cBuffer.save(filename, mParams, "");
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void Simulator::savePotential(size_t sliceindex, size_t,
                              const CBuffer2D& potential)
{
  Buffer2DInfo info(mParams.cols(), mParams.rows(),
                    mParams.width(), mParams.height(), Real());
  potential.save(QDir::homePath().toStdString() + "/Desktop/potential" +
                 std::to_string(sliceindex) + ".tiff", info, "");
}

void Simulator::displayWave()
{
  std::unique_lock<std::mutex> lock(mMutex);
  auto info = mExitWave.info();

  Real probability = mExitWave.absSqrReduce() / mExitWave.numPixels();
  probabilityLabel->setText(QString::number(probability, 'f', 4));

  QSize sizeWavefunction = calculateSize(exitWaveAbsSqrLabel->width(),
                                         mParams.aspectRatio());

  // show the absolute square of the exit wave function
  using namespace std::placeholders;
  auto extract = std::bind(extractAbsSqr, _1, info.minAbsSqr, info.maxAbsSqr);
  // FIXME
  auto imageAbsSqr = Image::fromCBuffer2D(mExitWave, (Image::ExtractorCColor)extract);
  auto qImageAbsSqr = imageToQImage(imageAbsSqr).scaled(sizeWavefunction,
                                                        Qt::IgnoreAspectRatio,
                                                        Qt::SmoothTransformation);
  exitWaveAbsSqrLabel->setPixmap(QPixmap::fromImage(qImageAbsSqr));

  // show the phase of the exit wave function
  // FIXME
  auto imagePhase = Image::fromCBuffer2D(mExitWave, (Image::ExtractorCColor)extractPhaseColor);
  auto qImagePhase = imageToQImage(imagePhase).scaled(sizeWavefunction,
                                                      Qt::IgnoreAspectRatio,
                                                      Qt::SmoothTransformation);
  exitWavePhaseLabel->setPixmap(QPixmap::fromImage(qImagePhase));

  // Calculate the FFT
  mExitWaveFft = mExitWave.fft(Direction::forward);

  // Create a copy
  auto exitWaveFft = mExitWaveFft.fftShift();

  QSize sizeWavefunctionFft = calculateSize(exitWaveFftLogAbsLabel->width(),
                                            mParams.kXNyquist() / mParams.kYNyquist());

  // show the logarithm of the absolute square of the Fourier transform of the
  // exit wave function
  auto infoFft = exitWaveFft.info();
  auto extractFft = std::bind(extractLog1pAbs, _1, std::log1p(infoFft.minAbs),
                              std::log1p(infoFft.maxAbs));
  // FIXME
  auto imageLogAbsFft = Image::fromCBuffer2D(exitWaveFft, (Image::ExtractorCColor)extractFft);
  auto qImageLogAbsFft = imageToQImage(imageLogAbsFft).scaled(sizeWavefunctionFft,
                                                              Qt::IgnoreAspectRatio,
                                                              Qt::SmoothTransformation);
  exitWaveFftLogAbsLabel->setPixmap(QPixmap::fromImage(qImageLogAbsFft));

  changeOptics();

  mSimulationRunning = false;

  // reactivate the simulate button
  simulateButton->setEnabled(true);
}

void Simulator::updateProgressBar(size_t part, size_t total)
{
  progressBar->setValue(part);
  progressBar->setMaximum(total - 1);
}

void Simulator::loadSample()
{
  try
  {
    if (sampleName.empty() || ("NOTSET" == sampleName))
      return;

    mSample->load(sampleName);

    if (showProjPotential)
    {
      const int64_t cols = projPotentialLabel->width();
      const int64_t rows = projPotentialLabel->height();

      Vector3R minBoxSample = mSample->minQuadraticBoundingBox(R(1.1));
      Vector3R maxBoxSample = mSample->maxQuadraticBoundingBox(R(1.1));

      auto projPotential =
        mSample->projectedPotential(cols, rows, minBoxSample, maxBoxSample, 3);

      auto info = projPotential.info();

      auto extract = [info] (Real x)
      {
        const Real scaledValue = (x - info.min) / (info.max - info.min);
        return Rgba8(saturate(scaledValue * R(255.0)));
      };

      mProjPotentialImage = Image::fromRBuffer2D(projPotential, extract);
      mProjPotentialPixmap =
        QPixmap::fromImage(imageToQImage(mProjPotentialImage));
    }

    numAtomsLabel->setText(QString::number(mSample->atoms().size()));
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void Simulator::loadSampleClicked()
{
  try
  {
    const QString title = tr("Load Sample");
    const QString filter = tr("Samples (*.pdb *.atom)");
    QString filename = QFileDialog::getOpenFileName(this, title, mSamplePath,
                                                    filter);
    if (filename.isEmpty())
      return;

    // extract the path and save it for the next dialog
    mSamplePath = filename.left(filename.lastIndexOf(QDir::separator()));

    sampleName.setValue(filename.toStdString());
    loadSample();
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }

  changeSimulation();
}

void Simulator::changeSimulation()
{
  if (mInterfaceUpdatePending)
    return;

  try
  {
    bandlimit.parse(bandlimitComboBox->currentText().toStdString());
    bFactor.setValue(static_cast<Real>(bFactorSpinBox->value()));
    cutoff.setValue(static_cast<Real>(cutoffSpinBox->value()));
    energy.setValue(static_cast<Real>(energySpinBox->value() * 1e3));
    cols.setValue(widthSpinBox->value());
    rows.setValue(heightSpinBox->value());
    deltaZ.setValue(static_cast<Real>(deltaZSpinBox->value()));
    minimumX.setValue(static_cast<Real>(minimumXSpinBox->value()));
    minimumY.setValue(static_cast<Real>(minimumYSpinBox->value()));
    minimumZ.setValue(static_cast<Real>(minimumZSpinBox->value()));
    maximumX.setValue(static_cast<Real>(maximumXSpinBox->value()));
    maximumY.setValue(static_cast<Real>(maximumYSpinBox->value()));
    maximumZ.setValue(static_cast<Real>(maximumZSpinBox->value()));
    numRepetitions.setValue(repetitionsSpinBox->value());
    formFactorName.setValue(formFactorComboBox->currentText().toStdString());
    phaseShiftName.setValue(phaseShiftComboBox->currentText().toStdString());
    potentialName.setValue(potentialComboBox->currentText().toStdString());
    propagatorName.setValue(propagatorComboBox->currentText().toStdString());

    // set the simulation parameters
    mParams = multisliceParamsFromCommandLine();

    // output some calculations
    wavelengthOutputLabel->setText(QString::number(mParams.wavelength() * R(1e2), 'f', 4));
    dimensionXLabel->setText(QString::number(mParams.width()));
    dimensionYLabel->setText(QString::number(mParams.height()));
    dimensionZLabel->setText(QString::number(mParams.depth()));
    nyquistXLabel->setText(QString::number(mParams.kXNyquist() / Math::twoPi, 'f', 4));
    nyquistYLabel->setText(QString::number(mParams.kYNyquist() / Math::twoPi, 'f', 4));
    pixelAspectRatioLabel->setText(QString::number(mParams.pixelAspectRatio(), 'f', 4));

    if (mSample)
    {
      Vector3R minBoxSample = mSample->minQuadraticBoundingBox(R(1.1));
      Vector3R maxBoxSample = mSample->maxQuadraticBoundingBox(R(1.1));

      Vector3R minBoxSimulation = mParams.minBox();
      Vector3R maxBoxSimulation = mParams.maxBox();

      Vector3R minBoxRescaled = (minBoxSimulation - minBoxSample) /
                                (maxBoxSample - minBoxSample);
      // the size of the simulation box rescaled to normalized coordinates
      Vector3R sizeBoxRescaled = (maxBoxSimulation - minBoxSimulation) /
                                 (maxBoxSample - minBoxSample);

      const int projPotCols = projPotentialLabel->width();
      const int projPotRows = projPotentialLabel->height();

      QPixmap pixmap = mProjPotentialPixmap;
      QPainter painter(&pixmap);
      painter.setPen(Qt::red);
      painter.drawRect(static_cast<int>(minBoxRescaled.x()  * projPotCols),
                       static_cast<int>(minBoxRescaled.y()  * projPotRows),
                       static_cast<int>(sizeBoxRescaled.x() * projPotCols),
                       static_cast<int>(sizeBoxRescaled.y() * projPotRows));

      projPotentialLabel->setPixmap(pixmap);

      numSlicesLabel->setText(QString::number(mParams.numSlices()));
    }

    changeOptics();
  }
  catch (Exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void Simulator::changeOptics()
{
  // prevent double entrance
  if (mInterfaceUpdatePending)
    return;

  mCtf.setAperture(static_cast<Real>(apertureSpinBox->value() * 1e-3));
  mCtf.setC(1, static_cast<Real>(defocusSpinBox->value() * 1e1));
  mCtf.setC(3, static_cast<Real>(c3SpinBox->value() * 1e4));
  mCtf.setC(5, static_cast<Real>(c5SpinBox->value() * 1e4));
  mCtf.setEnergy(static_cast<Real>(energySpinBox->value() * 1e3));

  // Plot the CTF
  // the maximal k value for the plot in 1 / nm
  const Real kkMax = mParams.kXNyquist() * R(1.05) / Math::twoPi * R(10.0);
  ctfPlot->setAxisScale(QwtPlot::xBottom, 0.0, kkMax);
  ctfPlot->setAxisScale(QwtPlot::xTop, 0.0, kkMax);

  QVector<double> dataX;
  QVector<double> dataImag;
  QVector<double> dataReal;

  const int numDataPoints = 1000;
  const Real wavenumber = mCtf.wavenumber();
  for (int i = 0; i < numDataPoints; ++i)
  {
    // electron microscopic spatial frequency in 1 / nm
    Real kk = kkMax * i / numDataPoints;
    // physical spatial frequency in 1 / Angstrom
    Real k = R(0.1) * kk * Math::twoPi;
    dataX.push_back(kk);
    dataReal.push_back(mCtf.value(k / wavenumber, 0.0).real());
    dataImag.push_back(mCtf.value(k / wavenumber, 0.0).imag());
  }

  mCurveReal.setSamples(dataX, dataReal);
  mCurveImag.setSamples(dataX, dataImag);
  ctfPlot->replot();

  // Simulate optics
  if (!mExitWaveFft.data())
    return;

  std::ostringstream description;
  description << "# Simulator, Version: " << version << "\n"
                 "# " << Time::now() << '\n'
              << CL::configString() << '\n';

  if (exitWaveRadioButton->isChecked())
  {
    view1->setBuffer(mExitWave.clone(), mParams, description.str(), false);
  }
  else if (exitWaveFftRadioButton->isChecked())
  {
    view1->setBuffer(mExitWaveFft.clone(), mParams, description.str(), true);
  }
  else if (opticRadioButton->isChecked())
  {
    // Apply the CTF to the exit wave
    auto ctfBuffer = mCtf.buffer(Buffer2DInfoBase(mParams.cols(), mParams.rows(),
                                 mParams.width(), mParams.height()),
                                 R(1.0) / (mParams.cols() * mParams.rows()));

    auto finalWaveFft = mExitWaveFft.clone();
    finalWaveFft *= ctfBuffer;

    view1->setBuffer(finalWaveFft.fft(Direction::backward), mParams,
                     description.str(), false);
  }
  else
    AURORA_UNREACHABLE;
}


void Simulator::updateInterface()
{
  mInterfaceUpdatePending = true;
  bandlimitComboBox->setCurrentText(QString::fromStdString(bandlimit.string()));
  bFactorSpinBox->setValue(bFactor);
  cutoffSpinBox->setValue(cutoff);
  deltaZSpinBox->setValue(deltaZ);
  energySpinBox->setValue(energy * R(1e-3));
  repetitionsSpinBox->setValue(numRepetitions);
  phaseShiftComboBox->setCurrentText(QString::fromStdString(phaseShiftName));
  potentialComboBox->setCurrentText(QString::fromStdString(potentialName));
  propagatorComboBox->setCurrentText(QString::fromStdString(propagatorName));
  formFactorComboBox->setCurrentText(QString::fromStdString(formFactorName));
  widthSpinBox->setValue(cols);
  heightSpinBox->setValue(rows);

  minimumXSpinBox->setValue(minimumX);
  minimumYSpinBox->setValue(minimumY);
  minimumZSpinBox->setValue(minimumZ);
  maximumXSpinBox->setValue(maximumX);
  maximumYSpinBox->setValue(maximumY);
  maximumZSpinBox->setValue(maximumZ);
  mInterfaceUpdatePending = false;
}

void Simulator::loadConfigFile()
{
  try
  {
    const QString title = tr("Load Configuration");
    const QString filter = tr("Configuration (*.txt)");
    const QString filename = QFileDialog::getOpenFileName(this, title,
                                                          mConfigPath, filter);

    if (filename.isEmpty())
      return;

    // extract the path and save it for the next dialog
    mConfigPath = filename.left(filename.lastIndexOf(QDir::separator()));

    CL::parseConfigFile(filename.toStdString());

    updateInterface();

    loadSample();
    changeSimulation();
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void Simulator::saveConfigFile()
{
  try
  {
    const QString caption = tr("Save Configuration");
    const QString filter  = tr("Configuration (*.txt)");

    const std::string filename =
      saveFileDialog(this, caption, filter, "txt", mConfigPath);

    if (filename.empty())
      return;

    std::ofstream configFile(filename);
    configFile << "# Simulator, Version: " << version << "\n"
                  "# " << Time::now() << '\n'
               << CL::configString() << '\n';
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void Simulator::extractBox()
{
  if (!mSample)
    return;

  Vector3R minBoxSample = mSample->minQuadraticBoundingBox(R(1.0));
  Vector3R maxBoxSample = mSample->maxQuadraticBoundingBox(R(1.0));

  mInterfaceUpdatePending = true;
  minimumXSpinBox->setValue(minBoxSample.x());
  minimumYSpinBox->setValue(minBoxSample.y());
  minimumZSpinBox->setValue(minBoxSample.z());
  maximumXSpinBox->setValue(maxBoxSample.x());
  maximumYSpinBox->setValue(maxBoxSample.y());
  maximumZSpinBox->setValue(maxBoxSample.z());
  mInterfaceUpdatePending = false;

  changeSimulation();
}

} // namespace Aurora

int main(int argc, char* argv[])
{
  try
  {
    QLocale::setDefault(QLocale::c());
    QApplication a(argc, argv);
    setlocale(LC_ALL, "C");
    if (!Aurora::CL::parse(argc, argv, "Aurora Simulator"))
      return EXIT_FAILURE;
    Aurora::Simulator simulator;
    simulator.show();
    return a.exec();
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
