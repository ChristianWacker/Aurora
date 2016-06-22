//--- Aurora/Clients/Opsys2/OpSys2.cpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "Clients/OpSys2/OpSys2.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Quadrature.hpp"
#include "AuroraLib/Utils.hpp"

#include "QtSupport/QtSupport.hpp"

#include <iostream>
#include <QApplication>
#include <QFileDialog>
#include <QLineEdit>
#include <qwt_plot_canvas.h>
#include <qwt_plot_grid.h>
#include <vector>

namespace Aurora
{

OpSys2::OpSys2(QMainWindow* parent) :
  QMainWindow(parent),
  mInterfaceUpdatePending(true),
  mCtf(1.0), // dummy energy value
  mLoadPath(QDir::homePath())
{
  setupUi(this);

  // abbreviations
  void (QDoubleSpinBox::*doubleChanged)(double) = &QDoubleSpinBox::valueChanged;
  void (QSpinBox::*intChanged)(int) = &QSpinBox::valueChanged;
  auto stateChanged = &QCheckBox::stateChanged;

  connect(apertureSpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(energySpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(maxFrequencySpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(cCSpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(semiAngleSpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(energySpreadSpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(htRippleSpinBox, doubleChanged, this, &OpSys2::refreshView);
  connect(currentInstabilitySpinBox, doubleChanged, this, &OpSys2::refreshView);

  connect(scherzerOrderSpinBox, intChanged, this, &OpSys2::refreshView);

  auto setToFocus = [=]
  {
    setRealCoefficient("C1", 0.0);
  };
  connect(focusButton, &QPushButton::clicked, setToFocus);

  auto setToScherzerDefocus = [=]
  {
    Real scherzerDefocus =
      mCtf.coherent().scherzerDefocus(scherzerOrderSpinBox->value());
    setRealCoefficient("C1", scherzerDefocus);
  };
  connect(scherzerButton, &QPushButton::clicked, setToScherzerDefocus);

  connect(ctfAngleSpinBox, doubleChanged, this, &OpSys2::refreshView);

  connect(envelopeCheckBox, stateChanged, this, &OpSys2::refreshView);
  connect(ctfRealCheckBox, stateChanged, this, &OpSys2::refreshView);
  connect(ctfImagCheckBox, stateChanged, this, &OpSys2::refreshView);
  connect(ctfWpoaRealCheckBox, stateChanged, this, &OpSys2::refreshView);
  connect(ctfWpoaImagCheckBox, stateChanged, this, &OpSys2::refreshView);

  connect(dataPointsSpinBox, intChanged, this, &OpSys2::refreshView);

  connect(loadExitWaveButton, &QPushButton::clicked,
          this, &OpSys2::showExitWaveLoadDialog);

  // setting up the CTF plot
  // make the plot background white including the label and suppress the border
  ctfPlot->setAutoFillBackground(true);
  ctfPlot->setPalette(Qt::white);
  ctfPlot->setCanvasBackground(Qt::white);
  auto canvas = static_cast<QwtPlotCanvas*>(ctfPlot->canvas());
  canvas->setFrameStyle(QFrame::NoFrame);

  // label the axes
  QwtText xTitleBottom("Spatial Frequency [1 / nm]");
  xTitleBottom.setFont(QFont("Times", 14));
  ctfPlot->setAxisTitle(QwtPlot::xBottom, xTitleBottom);

  QwtText xTitleTop("Scattering Angle [mrad]");
  xTitleTop.setFont(QFont("Times", 14));
  ctfPlot->setAxisTitle(QwtPlot::xTop, xTitleTop);
  ctfPlot->enableAxis(QwtPlot::xTop);

  ctfPlot->setAxisTitle(QwtPlot::yLeft, "CTF");
  ctfPlot->setAxisScale(QwtPlot::yLeft, -1.05, 1.05);

  // add an gray dotted grid
  QwtPlotGrid* grid = new QwtPlotGrid();
  grid->setPen(QPen(Qt::gray, 0, Qt::DotLine));
  grid->attach(ctfPlot);

  mCurveCtfReal.setPen(QPen(Qt::red));
  mCurveCtfImag.setPen(QPen(Qt::blue));
  mCurveEnvelopePos.setPen(QPen(Qt::black));
  mCurveEnvelopeNeg.setPen(QPen(Qt::black));
  mCurveNumCtfReal.setPen(QPen(Qt::green));
  mCurveNumCtfImag.setPen(QPen(Qt::magenta));

  mCurveCtfReal.attach(ctfPlot);
  mCurveCtfImag.attach(ctfPlot);
  mCurveEnvelopePos.attach(ctfPlot);
  mCurveEnvelopeNeg.attach(ctfPlot);
  mCurveNumCtfReal.attach(ctfPlot);
  mCurveNumCtfImag.attach(ctfPlot);

  // conversion factors to Angstrom
  mUnits["m"]  = 1e10;
  mUnits["nm"] = 1e1;
  mUnits["µm"] = 1e4;
  mUnits["mm"] = 1e7;

  // spherical aberrations
  addCoefficient("C1", "nm", false);
  addCoefficient("C3", "µm", false);
  addCoefficient("C5", "µm", false);
  addCoefficient("C7", "mm", false);

  // Astigmatism
  const char* aUnits[] = {"nm", "nm", "nm", "µm", "µm", "µm", "mm", "mm"};
  for (int i = 0; i < 8; ++i)
    addCoefficient('A' + QString::number(i), aUnits[i], true);

  addCoefficient("B2", "nm", true);
  addCoefficient("B4", "µm", true);
  addCoefficient("B6", "mm", true);

  addCoefficient("D4", "µm", true);
  addCoefficient("D6", "mm", true);

  addCoefficient("F6", "µm", true);

  addCoefficient("G7", "mm", true);

  addCoefficient("R5", "µm", true);
  addCoefficient("R7", "µm", true);

  addCoefficient("S3", "nm", true);
  addCoefficient("S5", "µm", true);
  addCoefficient("S7", "mm", true);

  tableWidget->resizeColumnsToContents();
  mInterfaceUpdatePending = false;

  refreshView();
}

void OpSys2::addCoefficient(const QString& name,
                            const QString& defaultUnit, bool withAngles)
{
  int row = tableWidget->rowCount();
  tableWidget->setRowCount(row + 1);

  QLabel* nameLabel = new QLabel(name, tableWidget);
  nameLabel->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
  tableWidget->setCellWidget(row, 0, nameLabel);

  QDoubleSpinBox* absSpinBox = new QDoubleSpinBox(tableWidget);
  void (QDoubleSpinBox::*doubleChanged)(double) = &QDoubleSpinBox::valueChanged;
  connect(absSpinBox, doubleChanged, this, &OpSys2::refreshView);
  tableWidget->setCellWidget(row, 1, absSpinBox);
  absSpinBox->setSingleStep(0.1);
  absSpinBox->setMaximum(+1e10);
  absSpinBox->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
  mConstants[name].abs = absSpinBox;

  // ComboBox for the unit
  QComboBox* unitBox = new QComboBox(tableWidget);
  void (QComboBox::*comboBoxChanged)(const QString&) = &QComboBox::currentIndexChanged;
  tableWidget->setCellWidget(row, 2, unitBox);

  unitBox->insertItem(0, "m");
  unitBox->insertItem(1, "mm");
  unitBox->insertItem(2, "µm");
  unitBox->insertItem(3, "nm");
  unitBox->setCurrentText(defaultUnit);

  mConstants[name].unitBox = unitBox;
  connect(unitBox, comboBoxChanged, this, &OpSys2::refreshView);

  if (withAngles)
  {
    absSpinBox->setMinimum(0.0);

    QDoubleSpinBox* angleSpinBox = new QDoubleSpinBox(tableWidget);
    connect(angleSpinBox, doubleChanged, this, &OpSys2::refreshView);
    tableWidget->setCellWidget(row, 3, angleSpinBox);
    angleSpinBox->setSingleStep(1.0);
    angleSpinBox->setMinimum(-1e4);
    angleSpinBox->setMaximum(+1e4);
    angleSpinBox->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    mConstants[name].angle = angleSpinBox;
  }
  else
  {
    absSpinBox->setMinimum(-1e10);
    mConstants[name].angle = nullptr;
  }

  tableWidget->setRowHeight(row, 24);;
}

Real OpSys2::getRealCoefficient(const QString& name)
{
  Real conversionFactor = mUnits[mConstants[name].unitBox->currentText()];
  return static_cast<Real>(mConstants[name].abs->value()) * conversionFactor;
}

void OpSys2::setRealCoefficient(const QString& name, Real value)
{
  mConstants[name].unitBox->setCurrentText("nm");
  Real conversionFactor = mUnits["nm"];
  mConstants[name].abs->setValue(value / conversionFactor);
}

Complex OpSys2::getComplexCoefficient(const QString& name)
{
  Real conversionFactor = mUnits[mConstants[name].unitBox->currentText()];
  return std::polar(
           static_cast<Real>(mConstants[name].abs->value()) * conversionFactor,
           static_cast<Real>(mConstants[name].angle->value()) * Math::degToRad);
}

void OpSys2::refreshView()
{
  if (mInterfaceUpdatePending)
    return;

  // get the values from the GUI controls
  for (size_t i = 0; i < 8; ++i)
    mCtf.coherent().setA(i, getComplexCoefficient('A' + QString::number(i)));

  mCtf.coherent().setB(2, getComplexCoefficient("B2"));
  mCtf.coherent().setB(4, getComplexCoefficient("B4"));
  mCtf.coherent().setB(6, getComplexCoefficient("B6"));

  mCtf.coherent().setD(4, getComplexCoefficient("D4"));
  mCtf.coherent().setD(6, getComplexCoefficient("D6"));

  mCtf.coherent().setF(6, getComplexCoefficient("F6"));

  mCtf.coherent().setG(7, getComplexCoefficient("G7"));

  mCtf.coherent().setR(5, getComplexCoefficient("R5"));
  mCtf.coherent().setR(7, getComplexCoefficient("R7"));

  mCtf.coherent().setS(3, getComplexCoefficient("S3"));
  mCtf.coherent().setS(5, getComplexCoefficient("S5"));
  mCtf.coherent().setS(7, getComplexCoefficient("S7"));

  mCtf.coherent().setAperture(static_cast<Real>(apertureSpinBox->value() * 1e-3));
  mCtf.coherent().setEnergy(static_cast<Real>(energySpinBox->value() * 1e3));
  mCtf.setIlluminationSemiAngle(static_cast<Real>(semiAngleSpinBox->value() * 1e-3));
  mCtf.setCC(static_cast<Real>(cCSpinBox->value() * 1e4));
  mCtf.setEnergySpread(static_cast<Real>(energySpreadSpinBox->value()));
  mCtf.setHighTensionRipple(static_cast<Real>(htRippleSpinBox->value() * 1e-6));
  mCtf.setCurrentInstability(static_cast<Real>(currentInstabilitySpinBox->value() * 1e-6));

  mCtf.coherent().setC(1, getRealCoefficient("C1"));
  mCtf.coherent().setC(3, getRealCoefficient("C3"));
  mCtf.coherent().setC(5, getRealCoefficient("C5"));
  mCtf.coherent().setC(7, getRealCoefficient("C7"));

  using namespace Math;
  const Real angle = static_cast<Real>(ctfAngleSpinBox->value()) * degToRad;

  // Update data fields
  const Real wavenumber = mCtf.wavenumber();
  const Real wavelength = Math::twoPi / wavenumber;
  wavelengthLabel->setText(QString::number(R(1e2) * wavelength, 'f', 3));
  focusSpreadLabel->setText(QString::number(R(0.1) * mCtf.focusSpread(), 'f', 3));
  const size_t scherzerOrder = scherzerOrderSpinBox->value();
  const Real scherzerDefocus = mCtf.coherent().scherzerDefocus(scherzerOrder);
  scherzerDefocusLabel->setText(QString::number(R(0.1) * scherzerDefocus, 'f', 1));

  // Update 1D plot
  // Maximal spatial frequency for the plot in 1 / nm
  const Real maxSpatialFreqExp = static_cast<Real>(maxFrequencySpinBox->value());

  // Maximal spatial frequency in physical convention in 1 / Angstrom
  const Real maxSpatialFreqPhys = maxSpatialFreqExp * R(0.1) * Math::twoPi;

  // The maximal scattering angle in rad
  const Real maxScatteringAngle = maxSpatialFreqPhys / wavenumber;

  ctfPlot->setAxisScale(QwtPlot::xBottom, 0, maxSpatialFreqExp);
  ctfPlot->setAxisScale(QwtPlot::xTop, 0, maxScatteringAngle * R(1e3));

  QVector<double> dataX;
  QVector<double> dataCtfReal;
  QVector<double> dataCtfImag;
  QVector<double> envelopeDataPos;
  QVector<double> envelopeDataNeg;

  QVector<double> numericalDataCtfReal;
  QVector<double> numericalDataCtfImag;

  // Hermite polynomial for the (numerical) Gauss-Hermite-Quadrature
  HermitePolynomial hermite(dataPointsSpinBox->value());
  const Real aperture = mCtf.coherent().aperture();

  const size_t numDataPoints = 650;
  for (size_t i = 0; i < numDataPoints; ++i)
  {
    // Experimental spatial frequency in 1 / nm
    const Real spatialFrequencyExp = maxSpatialFreqExp * i / numDataPoints;
    // Physical wavenumber in 1 / Angstrom
    const Real spatialFrequencyPhys = maxSpatialFreqPhys * i / numDataPoints;

    const Real theta = spatialFrequencyPhys / wavenumber;

    // aperture
    if (theta > aperture)
    {
      numericalDataCtfReal.push_back(0.0);
      numericalDataCtfImag.push_back(0.0);
      continue;
    }

    const Real thetaX = theta * std::cos(angle);
    const Real thetaY = theta * std::sin(angle);

    const Real envelope = mCtf.spatialEnvelope(theta) *
                          mCtf.temporalEnvelope(theta);

    Complex transfer = mCtf.coherent().pupilAngle(thetaX, thetaY);

    dataX.push_back(spatialFrequencyExp);
    dataCtfReal.push_back(transfer.real() * envelope);
    dataCtfImag.push_back(transfer.imag() * envelope);
    envelopeDataPos.push_back(envelope);
    envelopeDataNeg.push_back(-envelope);

    // skip the expensive numerical calculations, if they are not needed
    if (!ctfRealCheckBox->checkState() && !ctfImagCheckBox->checkState())
      continue;

    using namespace std::placeholders;
    auto integrand1 = [&] (Real betaX, Real betaY)
    {
      auto f = std::bind(&CtfCoherent::pupilAngle, mCtf.coherent(),
                         thetaX + betaX, thetaY + betaY, _1);
      return integrateGaussHermite<Complex>(f,
                                            hermite, 0.0,
                                            powerOf<2>(mCtf.focusSpread()));
    };

    auto integrand2 = [&] (Real betaY)
    {
      return integrateGaussHermite<Complex>(std::bind(integrand1, _1, betaY),
                                            hermite, 0.0,
                                            powerOf<2>(mCtf.sigmaBeta()));
    };

    Complex integral = integrateGaussHermite<Complex>(integrand2, hermite, 0.0,
                                                      powerOf<2>(mCtf.sigmaBeta()));
    numericalDataCtfReal.push_back(integral.real());
    numericalDataCtfImag.push_back(integral.imag());
  }

  // Calculate the first zero of the imaginary part of the CTF
  // (point-resolution).
  QVector<double>* data;
  // Use the numerical data if available, since it includes all effects.
  if (!ctfRealCheckBox->checkState() && !ctfImagCheckBox->checkState())
    data = &dataCtfImag;
  else
    data = &numericalDataCtfImag;

  int signFlipIndex = 0;
  for (; signFlipIndex < data->size() - 1; ++signFlipIndex)
  {
    if ((*data)[signFlipIndex] * (*data)[signFlipIndex + 1] < 0.0)
      break;
  }

  if (signFlipIndex == numDataPoints - 1)
    pointResolutionLabel->setText("?");
  else
    pointResolutionLabel->setText(QString::number(dataX[signFlipIndex], 'f', 2));

  // set the data for the 1D plot
  mCurveCtfReal.setSamples(dataX, dataCtfReal);
  mCurveCtfImag.setSamples(dataX, dataCtfImag);
  mCurveEnvelopePos.setSamples(dataX, envelopeDataPos);
  mCurveEnvelopeNeg.setSamples(dataX, envelopeDataNeg);
  mCurveNumCtfReal.setSamples(dataX, numericalDataCtfReal);
  mCurveNumCtfImag.setSamples(dataX, numericalDataCtfImag);

  mCurveCtfReal.setVisible(ctfWpoaRealCheckBox->checkState());
  mCurveCtfImag.setVisible(ctfWpoaImagCheckBox->checkState());
  mCurveEnvelopePos.setVisible(envelopeCheckBox->checkState());
  mCurveEnvelopeNeg.setVisible(envelopeCheckBox->checkState());
  mCurveNumCtfReal.setVisible(ctfRealCheckBox->checkState());
  mCurveNumCtfImag.setVisible(ctfImagCheckBox->checkState());

  ctfPlot->replot();

  // Update the 2D CTF plots
  const int cols = ctfRealLabel->width();
  const int rows = ctfRealLabel->height();

  auto ctfBuffer = mCtf.buffer(Buffer2DInfoBase(cols, rows,
                               Math::pi * cols / maxSpatialFreqPhys,
                               Math::pi * rows / maxSpatialFreqPhys)).fftShift();

  QImage imageCtfReal(cols, rows, QImage::Format_RGB32);
  QImage imageCtfImag(cols, rows, QImage::Format_RGB32);

  #pragma omp parallel for
  for (int y = 0; y < rows; ++y)
  {
    QRgb* dstReal = reinterpret_cast<QRgb*>(imageCtfReal.scanLine(y));
    QRgb* dstImag = reinterpret_cast<QRgb*>(imageCtfImag.scanLine(y));

    for (int64_t x = 0; x < cols; ++x)
    {
      // Map the range [-1, 1] to [0, 255]. Invert y coordinate, to correct for
      // the different definitions of the origin.
      uint8_t realValue = static_cast<uint8_t>(R(127.5) + R(127.4) * ctfBuffer.pixel(x, rows - y - 1).real());
      uint8_t imagValue = static_cast<uint8_t>(R(127.5) + R(127.4) * ctfBuffer.pixel(x, rows - y - 1).imag());

      dstReal[x] = qRgb(realValue, realValue, realValue);
      dstImag[x] = qRgb(imagValue, imagValue, imagValue);
    }
  }

  // add red direction line
  QPixmap ctfPixmapReal = QPixmap::fromImage(imageCtfReal);
  QPixmap ctfPixmapImag = QPixmap::fromImage(imageCtfImag);

  QPainter painterReal(&ctfPixmapReal);
  QPainter painterImag(&ctfPixmapImag);

  painterReal.setPen(Qt::red);
  painterImag.setPen(Qt::blue);

  const int centerX = cols / 2;
  const int centerY = rows / 2;
  const int radius  = cols / 2;

  const int lineEndX = static_cast<int>(centerX + radius * std::cos(angle));
  const int lineEndY = static_cast<int>(centerY - radius * std::sin(angle));

  painterReal.drawLine(centerX, centerY, lineEndX, lineEndY);
  painterImag.drawLine(centerX, centerY, lineEndX, lineEndY);

  ctfRealLabel->setPixmap(ctfPixmapReal);
  ctfImagLabel->setPixmap(ctfPixmapImag);

  // apply the CTF to the exit wave
  if ((0 != mExitWave.cols()) && (0 != mExitWave.rows()))
  {
    CBuffer2D waveFourier = mExitWave.fft(Direction::forward);
    CBuffer2D ctfBuffer = mCtf.buffer(mExitWaveInfo,
                                      1 / mExitWaveInfo.numPixels<Real>());
    waveFourier *= ctfBuffer;

    view->setBuffer(waveFourier.fft(Direction::backward), mExitWaveInfo,
                    mDescription);
  }
}

void OpSys2::setBuffer(const CBuffer2D& buffer, const Buffer2DInfo& info,
                       const std::string& description)
{
  mExitWave = buffer.clone();
  mExitWaveInfo = info;
  mDescription = description;

  auto pixelInfo = mExitWave.info();
  using namespace std::placeholders;
  auto extract = std::bind(extractAbsSqr, _1, pixelInfo.minAbsSqr, pixelInfo.maxAbsSqr);
  // FIXME:
  Image exitWaveImage = Image::fromCBuffer2D(mExitWave, (Image::ExtractorCColor)extract);
  QSize sizeExitWave = calculateSize(exitWaveLabel->width(),
                                     info.aspectRatio());
  QImage exitWaveQImage = imageToQImage(exitWaveImage).scaled(sizeExitWave,
                                                              Qt::IgnoreAspectRatio,
                                                              Qt::SmoothTransformation);
  exitWaveLabel->setPixmap(QPixmap::fromImage(exitWaveQImage));

  mInterfaceUpdatePending = true;
  energySpinBox->setValue(mExitWaveInfo.energy() * R(1e-3));
  Real maxFrequency = std::hypot(mExitWaveInfo.kXNyquist(),
                                 mExitWaveInfo.kYNyquist());
  maxFrequencySpinBox->setValue(maxFrequency * 10 / Math::twoPi);
  mInterfaceUpdatePending = false;

  refreshView();
}

void OpSys2::loadExitWave(const QString& filename)
{
  try
  {
    if (filename.isEmpty())
      return;

    mLoadPath = filename.left(filename.lastIndexOf(QDir::separator()));

    Buffer2DInfo info;
    std::string description;
    CBuffer2D buffer(filename.toStdString(), info, description);

    setBuffer(buffer, info, description);
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void OpSys2::showExitWaveLoadDialog()
{
  const QString title = tr("Load ComplexBuffer");
  const QString filter = tr("ComplexBuffer (*.tif *.tiff)");
  const QString filename =
    QFileDialog::getOpenFileName(this, title, mLoadPath, filter);

  loadExitWave(filename);
}

void OpSys2::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/plain"))
    event->acceptProposedAction();
}

void OpSys2::dropEvent(QDropEvent* event)
{
  QUrl url = QUrl::fromEncoded(event->mimeData()->text().trimmed().toUtf8());
  if (!url.isLocalFile())
    return;
  QString filename = url.toLocalFile();

  event->accept();
  loadExitWave(filename);
}

} // namespace Aurora

int main(int argc, char* argv[])
{
  try
  {
    QLocale::setDefault(QLocale::c());
    QApplication a(argc, argv);
    setlocale(LC_ALL, "C");
    Aurora::OpSys2 v;
    v.show();
    return a.exec();
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
