//--- Aurora/QtSupport/QView.cpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "QtSupport/QView.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Histogram.hpp"
#include "AuroraLib/Utils.hpp"
#include "QtSupport/QtSupport.hpp"

#include <fstream>
#include <iostream>
#include <QApplication>
#include <QFileDialog>
#include <QFrame>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QMenu>
#include <QPushButton>
#include <QUrl>
#include <QVBoxLayout>
#include <qwt_plot_canvas.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_renderer.h>

namespace Aurora
{

QView::QView(QWidget* parent, bool intern) :
  QWidget(parent),
  mIntern(intern),
  mInterfaceUpdatePending(false),
  mAspectRatio(1.0),
  mFourier(false),
  mLoadPath(QDir::homePath()),
  mSavePath(QDir::homePath()),
  mHistogramSavePath(QDir::homePath()),
  mProbability(0.0),
  mHistogramText(nullptr)
{
  if (!intern)
    setAcceptDrops(true);

  QVBoxLayout* verticalLayout0 = new QVBoxLayout(this);
  verticalLayout0->setMargin(0);

  QFrame* frame = new QFrame();
  frame->setFixedSize(532, 532);
  frame->setFrameShape(QFrame::Shape::Panel);
  frame->setFrameShadow(QFrame::Shadow::Raised);
  verticalLayout0->addWidget(frame);

  mImageLabel = new QLabel(frame);

  frame->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(frame, &QFrame::customContextMenuRequested,
          this, &QView::showBufferContextMenu);
  mImageLabel->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(mImageLabel, &QLabel::customContextMenuRequested,
          this, &QView::showBufferContextMenu);

  QHBoxLayout* horizontalLayout2 = new QHBoxLayout();
  horizontalLayout2->setMargin(0);
  verticalLayout0->addLayout(horizontalLayout2);
  verticalLayout0->setAlignment(Qt::AlignTop);

  mHistogramPlot = new QwtPlot();
  verticalLayout0->addWidget(mHistogramPlot);

  mHistogramPlot->setFixedSize(532, 180);
  // make the plot background white and suppress the border
  mHistogramPlot->setAutoFillBackground(true);
  mHistogramPlot->setPalette(Qt::white);
  mHistogramPlot->setCanvasBackground(Qt::white);
  auto canvas = static_cast<QwtPlotCanvas*>(mHistogramPlot->canvas());
  canvas->setFrameStyle(QFrame::NoFrame);

  // add an gray dotted grid
  QwtPlotGrid* grid = new QwtPlotGrid();
  grid->setPen(QPen(Qt::gray, 0, Qt::DotLine));
  grid->attach(mHistogramPlot);

  mHistogramSigma = new QwtPlotCurve();
  mHistogramSigma->setBrush(QBrush(QColor(255, 200, 200)));
  mHistogramSigma->setBaseline(Qt::Horizontal);
  mHistogramSigma->attach(mHistogramPlot);

  mHistogramCurve = new QwtPlotCurve();
  mHistogramCurve->setPen(QPen(Qt::blue));
  mHistogramCurve->attach(mHistogramPlot);

  mHistogramMode = new QwtPlotCurve();
  mHistogramMode->attach(mHistogramPlot);

  mHistogramMean = new QwtPlotCurve();
  mHistogramMean->setPen(QPen(Qt::red));
  mHistogramMean->attach(mHistogramPlot);

  mHistogramPlot->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(mHistogramPlot, &QwtPlot::customContextMenuRequested,
          this, &QView::showHistogramContextMenu);

  QVBoxLayout* verticalLayout1 = new QVBoxLayout();
  horizontalLayout2->addLayout(verticalLayout1);
  verticalLayout1->setSpacing(5);
  verticalLayout1->setAlignment(Qt::AlignTop);

  QVBoxLayout* verticalLayout2 = new QVBoxLayout();
  horizontalLayout2->addLayout(verticalLayout2);
  verticalLayout2->setSpacing(0);
  verticalLayout2->setAlignment(Qt::AlignTop);

  QVBoxLayout* verticalLayout3 = new QVBoxLayout();
  horizontalLayout2->addLayout(verticalLayout3);
  verticalLayout3->setSpacing(0);
  verticalLayout3->setAlignment(Qt::AlignTop);

  QVBoxLayout* verticalLayout4 = new QVBoxLayout();
  horizontalLayout2->addLayout(verticalLayout4);
  verticalLayout4->setSpacing(0);
  verticalLayout4->setAlignment(Qt::AlignTop);

  QGroupBox* preprocessRadioGroup = new QGroupBox("Preprocess (1)");
  verticalLayout2->addWidget(preprocessRadioGroup);
  preprocessRadioGroup->setMaximumSize(110, QWIDGETSIZE_MAX);

  QVBoxLayout* layoutPreprocess = new QVBoxLayout();
  preprocessRadioGroup->setLayout(layoutPreprocess);
  layoutPreprocess->setSpacing(0);
  layoutPreprocess->setAlignment(Qt::AlignTop);

  // abbreviations
  void (QDoubleSpinBox::*doubleChanged)(double) = &QDoubleSpinBox::valueChanged;
  auto buttonPushed = &QAbstractButton::clicked;
  auto stateChanged = &QCheckBox::stateChanged;

  mFftCheckBox = new QCheckBox("FFT (2)");
  verticalLayout4->addWidget(mFftCheckBox);
  connect(mFftCheckBox, stateChanged, this, &QView::updateView);

  mLogCheckBox = new QCheckBox("Log (4)");
  verticalLayout4->addWidget(mLogCheckBox);
  connect(mLogCheckBox, stateChanged, this, &QView::updateView);

  mLog1pCheckBox = new QCheckBox("Log1p (5)");
  verticalLayout4->addWidget(mLog1pCheckBox);
  connect(mLog1pCheckBox, stateChanged, this, &QView::updateView);

  mFilterCheckBox = new QCheckBox("Filter");
  verticalLayout4->addWidget(mFilterCheckBox);
  mFilterCheckBox->setChecked(true);
  connect(mFilterCheckBox, stateChanged, this, &QView::updateView);

  QGroupBox* representationRadioGroup = new QGroupBox("Representation (3)");
  verticalLayout3->addWidget(representationRadioGroup);

  QHBoxLayout* layoutRepresentation = new QHBoxLayout();
  representationRadioGroup->setLayout(layoutRepresentation);
  layoutRepresentation->setSpacing(0);
  layoutRepresentation->setAlignment(Qt::AlignTop);

  QVBoxLayout* layoutRepresentationLeft = new QVBoxLayout();
  layoutRepresentation->addLayout(layoutRepresentationLeft);
  layoutRepresentationLeft->setSpacing(0);
  layoutRepresentationLeft->setAlignment(Qt::AlignTop);

  QVBoxLayout* layoutRepresentationRight = new QVBoxLayout();
  layoutRepresentation->addLayout(layoutRepresentationRight);
  layoutRepresentationRight->setSpacing(0);
  layoutRepresentationRight->setAlignment(Qt::AlignTop);

  if (!intern)
  {
    QPushButton* loadBufferButton = new QPushButton("Load ComplexBuffer");
    verticalLayout1->addWidget(loadBufferButton);
    connect(loadBufferButton, buttonPushed, this, &QView::showBufferLoadDialog);
  }

  mAutoContrastCheckbox = new QCheckBox("Auto Contrast");
  verticalLayout1->addWidget(mAutoContrastCheckbox);
  connect(mAutoContrastCheckbox, stateChanged, this, &QView::updateView);
  mAutoContrastCheckbox->setChecked(true);

  mMinSpinBox = new QDoubleSpinBox();
  mMinSpinBox->setSingleStep(0.1);
  mMinSpinBox->setDecimals(20);
  mMinSpinBox->setRange(-1e18, 1e18);
  mMinSpinBox->setAlignment(Qt::AlignRight);
  verticalLayout1->addWidget(mMinSpinBox);
  connect(mMinSpinBox, doubleChanged, this, &QView::updateView);

  mMaxSpinBox = new QDoubleSpinBox();
  mMaxSpinBox->setSingleStep(0.1);
  mMaxSpinBox->setDecimals(20);
  mMaxSpinBox->setRange(-1e18, 1e18);
  mMaxSpinBox->setAlignment(Qt::AlignRight);
  verticalLayout1->addWidget(mMaxSpinBox);
  connect(mMaxSpinBox, doubleChanged, this, &QView::updateView);

  mNoneRadioButton = new QRadioButton("None");
  layoutPreprocess->addWidget(mNoneRadioButton);
  connect(mNoneRadioButton, buttonPushed, this, &QView::updateView);
  mNoneRadioButton->setChecked(true);

  mAbsRadioButton = new QRadioButton("Abs");
  layoutPreprocess->addWidget(mAbsRadioButton);
  connect(mAbsRadioButton, buttonPushed, this, &QView::updateView);

  mAbsSqrRadioButton = new QRadioButton("AbsSqr");
  layoutPreprocess->addWidget(mAbsSqrRadioButton);
  connect(mAbsSqrRadioButton, buttonPushed, this, &QView::updateView);

  mRealButton = new QRadioButton("Real");
  layoutRepresentationLeft->addWidget(mRealButton);
  connect(mRealButton, buttonPushed, this, &QView::updateView);

  mImagButton = new QRadioButton("Imag");
  layoutRepresentationLeft->addWidget(mImagButton);
  connect(mImagButton, buttonPushed, this, &QView::updateView);

  mAbsButton = new QRadioButton("Abs");
  layoutRepresentationLeft->addWidget(mAbsButton);
  connect(mAbsButton, buttonPushed, this, &QView::updateView);
  mAbsButton->setChecked(true);

  mAbsSqrButton = new QRadioButton("AbsSqr");
  layoutRepresentationRight->addWidget(mAbsSqrButton);
  connect(mAbsSqrButton, buttonPushed, this, &QView::updateView);

  mPhaseButton = new QRadioButton("Phase");
  layoutRepresentationRight->addWidget(mPhaseButton);
  connect(mPhaseButton, buttonPushed, this, &QView::updateView);

  mPhaseColorButton = new QRadioButton("Phase (Color)");
  layoutRepresentationRight->addWidget(mPhaseColorButton);
  connect(mPhaseColorButton, buttonPushed, this, &QView::updateView);
}

void QView::showBufferContextMenu(const QPoint& position)
{
  QPoint globalPosition =
    static_cast<QWidget*>(QObject::sender())->mapToGlobal(position);

  QMenu menu;
  QAction saveBuffer("Save ComplexBuffer...", this);
  QAction saveAsText("Save processed buffer as text...", this);
  QAction loadImage("Load Image...", this);
  QAction saveImage("Save Image...", this);

  if (mIntern)
  {
    menu.addAction(&saveBuffer);
    connect(&saveBuffer, &QAction::triggered, this, &QView::saveBuffer);
  }
  else
  {
    menu.addAction(&loadImage);
    connect(&loadImage, &QAction::triggered, this, &QView::loadImage);
  }

  menu.addAction(&saveAsText);
  connect(&saveAsText, &QAction::triggered, this, &QView::saveBufferAsText);

  menu.addAction(&saveImage);
  connect(&saveImage, &QAction::triggered, this, &QView::saveImage);

  menu.exec(globalPosition);
}

void QView::showHistogramContextMenu(const QPoint& position)
{
  if (!mHistogram)
    return;

  QPoint globalPosition = mHistogramPlot->mapToGlobal(position);

  QMenu menu;
  QAction saveHistogramImage("Save Image...", this);
  QAction saveHistogramData("Save Data...", this);
  menu.addAction(&saveHistogramImage);
  menu.addAction(&saveHistogramData);
  connect(&saveHistogramImage, &QAction::triggered,
          this, &QView::saveHistogramImage);
  connect(&saveHistogramData, &QAction::triggered,
          this, &QView::saveHistogramData);
  menu.exec(globalPosition);
}

void QView::saveHistogramImage()
{
  if (!mHistogram)
    return;

  try
  {
    const QString caption = tr("Save histogram image");
    const QString filter  = tr("Image formats (*.*)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "pdf", mHistogramSavePath);

    if (filename.empty())
      return;

    QwtPlotRenderer renderer;
    renderer.renderDocument(mHistogramPlot, QString::fromStdString(filename),
                            QSizeF(200, 100), 300);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void QView::saveHistogramData()
{
  if (!mHistogram)
    return;

  try
  {
    const QString caption = tr("Save histogram data");
    const QString filter  = tr("Textfiles (*.txt)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "txt", mHistogramSavePath);

    if (filename.empty())
      return;

    std::ofstream outFile(filename);

    outFile << "# Mean: " << mHistogram->mean() << '\n';
    outFile << "# Mode: " << mHistogram->mode() << '\n';
    outFile << "# Sigma: " << sqrt(mHistogram->variance()) << '\n';
    for (size_t i = 0, numBins = mHistogram->numBins(); i < numBins; ++i)
      outFile << mHistogram->binToValue(i) << ' ' << (*mHistogram)[i] << '\n';
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void QView::drawHistogram()
{
  if ((!std::isfinite(mHistogram->minValue())) ||
      (!std::isfinite(mHistogram->maxValue())))
    return;

  // data
  QVector<double> histogramDataX;
  QVector<double> histogramDataY;
  for (size_t i = 0, numBins = mHistogram->numBins(); i < numBins; ++i)
  {
    histogramDataX.push_back(mHistogram->binToValue(i));
    histogramDataY.push_back((*mHistogram)[i]);
  }
  mHistogramCurve->setSamples(histogramDataX, histogramDataY);

  // mode
  const size_t mode = mHistogram->mode();
  const Real modeX = mHistogram->binToValue(mode);
  const Real modeY = (*mHistogram)[mode];
  QVector<double> modeDataX = {modeX, modeX};
  QVector<double> modeDataY = {0, modeY};
  mHistogramMode->setSamples(modeDataX, modeDataY);

  // mean
  const Real mean = mHistogram->mean();
  const Real meanY = (*mHistogram)[mHistogram->valueToBin(mean)];
  QVector<double> meanDataX = {mean, mean};
  QVector<double> meanDataY = {0, meanY};
  mHistogramMean->setSamples(meanDataX, meanDataY);

  // sigma environment
  const Real variance = mHistogram->variance();
  const size_t sigmaBinFirst = mHistogram->valueToBin(mean - std::sqrt(variance));
  const size_t sigmaBinLast  = mHistogram->valueToBin(mean + std::sqrt(variance));
  QVector<double> sigmaDataX;
  QVector<double> sigmaDataY;
  for (size_t i = sigmaBinFirst; i <= sigmaBinLast; ++i)
  {
    sigmaDataX.push_back(mHistogram->binToValue(i));
    sigmaDataY.push_back((*mHistogram)[i]);
  }
  mHistogramSigma->setSamples(sigmaDataX, sigmaDataY);

  mHistogramPlot->setAxisScale(QwtPlot::xBottom,
                               mHistogram->minValue(), mHistogram->maxValue());
  mHistogramPlot->setAxisScale(QwtPlot::yLeft, 0, modeY);

  QwtText text("Mode: " + QString::number(modeX, 'f', 8) + "\n"
               "Mean: " + QString::number(mean, 'f', 8) + "\n"
               "Std Dev: " + QString::number(sqrt(variance), 'f', 8) + "\n"
               "Probability: " + QString::number(mProbability, 'f', 8));
  text.setRenderFlags(Qt::AlignTop | Qt::AlignRight);

  if (!mHistogramText)
    mHistogramText = new QwtPlotTextLabel();

  mHistogramText->setText(text);
  mHistogramText->attach(mHistogramPlot);
  mHistogramPlot->replot();
}

void QView::loadBuffer(const QString& filename)
{
  try
  {
    if (filename.isEmpty())
      return;

    // extract the path and save it for the next dialog
    mLoadPath = filename.left(filename.lastIndexOf(QDir::separator()));

    std::string description;
    Buffer2DInfo info;
    CBuffer2D buffer;
    buffer.load(filename.toStdString(), info, description);

    setBuffer(std::move(buffer), info, description);
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }

  emit(bufferChanged());
}

void QView::showBufferLoadDialog()
{
  try
  {
    const QString title = tr("Load ComplexBuffer");
    const QString filter = tr("ComplexBuffer (*.tif *.tiff)");
    const QString filename = QFileDialog::getOpenFileName(this, title,
                                                          mLoadPath, filter);
    loadBuffer(filename);
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void QView::loadImage()
{
  try
  {
    const QString title = tr("Load Image");
    const QString filter = tr("Images (*.tif *.tiff)");
    const QString filename = QFileDialog::getOpenFileName(this, title,
                                                          mLoadPath, filter);

    if (filename.isEmpty())
      return;

    // extract the path and save it for the next dialog
    mLoadPath = filename.left(filename.lastIndexOf(QDir::separator()));

    Image image;
    image.load(filename.toStdString());
    Buffer2DInfo bufferInfo(image.cols(), image.rows(), R(0.0), R(0.0), R(0.0));
    setBuffer(CBuffer2D::fromImage(image), bufferInfo, "");
  }
  catch (const Exception& e)
  {
    std::cerr << e.string() << std::endl;
  }
}

void QView::saveBuffer()
{
  if (!mBuffer.data())
    return;

  try
  {
    const QString caption = tr("Save Buffer");
    const QString filter  = tr("ComplexBuffer (*.tif *.tiff)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "tiff", mSavePath);

    if (!filename.empty())
      mBuffer.save(filename, mBufferInfo, mDescription);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void QView::saveBufferAsText()
{
  if (!mProcessedBuffer.data())
    return;

  try
  {
    const QString caption = tr("Save Float");
    const QString filter  = tr("Text file (*.txt)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "txt", mSavePath);

    if (filename.empty())
      return;

    const int64_t cols = mBuffer.cols();
    const int64_t rows = mBuffer.rows();

    std::ofstream outFile(filename);
    outFile << "# cols: " << cols << ", rows: " << rows << "\n"
               "# x        y         value\n";

    for (int64_t y = 0; y < rows; ++y)
    {
      outFile << '\n';
      for (int64_t x = 0; x < cols; ++x)
      {
        outFile << x << ' ' << y << ' '
                << mProcessedBuffer.pixel(x, y).real() << '\n';
      }
    }
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void QView::saveImage()
{
  if (!mImage.dataPtr())
    return;

  try
  {
    const QString caption = tr("Save Image");
    const QString filter  = tr("Images (*.tif *.tiff)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "tiff", mSavePath);

    if (!filename.empty())
      mImage.save(filename);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void QView::updateView()
{
  if (mInterfaceUpdatePending)
    return;

  if (!mBuffer.data())
    return;

  CBuffer2D buffer;

  // first processing step
  if (mAbsRadioButton->isChecked())
    buffer = mBuffer.bufferAbs();
  else if (mAbsSqrRadioButton->isChecked())
    buffer = mBuffer.bufferAbsSqr();
  else if (mNoneRadioButton->isChecked())
    buffer = mBuffer.clone();
  else
    AURORA_UNREACHABLE;

  // second processing step: FFT
  if (mFftCheckBox->isChecked())
  {
    Direction direction = mFourier ? Direction::backward : Direction::forward;
    buffer = buffer.fft(direction);
  }

  /// In Fourier space we want the lowest frequencies to be in the middle
  if (mFourier != mFftCheckBox->isChecked())
    buffer = buffer.fftShift();

  // third processing step: reduce from complex values to real values
  std::function<Real (const Complex&)> extract;

  if (mRealButton->isChecked())
    extract = Aurora::real;
  else if (mImagButton->isChecked())
    extract = Aurora::imag;
  else if (mAbsButton->isChecked())
    extract = Aurora::abs;
  else if (mAbsSqrButton->isChecked())
    extract = Aurora::absSqr;
  else if (mPhaseButton->isChecked() || mPhaseColorButton->isChecked())
    extract = Aurora::arg;
  else
    AURORA_UNREACHABLE;

  const int64_t cols = mBuffer.cols();
  const int64_t rows = mBuffer.rows();

  // save the new values in the real component of the pixels
  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      buffer.pixel(x, y).real(extract(buffer.pixel(x, y)));

  // fourth processing step: logarithm log|x|
  if (mLogCheckBox->isChecked())
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        buffer.pixel(x, y).real(std::log(std::abs(buffer.pixel(x, y))));
  }

  // fifth processing step: logarithm log(1 + |x|)
  if (mLog1pCheckBox->isChecked())
  {
    #pragma omp parallel for
    for (int64_t y = 0; y < rows; ++y)
      for (int64_t x = 0; x < cols; ++x)
        buffer.pixel(x, y).real(std::log1p(std::abs(buffer.pixel(x, y))));
  }

  mProcessedBuffer = buffer.clone();

  auto info = buffer.info();

  // catch exceptional conditions
  if ((info.minReal == info.maxReal) || !std::isfinite(info.minReal) ||
      !std::isfinite(info.maxReal))
  {
    info.minReal = 0;
    info.maxReal = 1;
  }

  // create a histogram
  mHistogram = std::make_unique<Histogram>(info.minReal, info.maxReal, 256);
  for (int64_t y = 0; y < rows; ++y)
    for (int64_t x = 0; x < cols; ++x)
      mHistogram->insert(buffer.pixel(x, y).real());

  if (mAutoContrastCheckbox->isChecked())
  {
    mInterfaceUpdatePending = true;
    mMinSpinBox->setValue(info.minReal);
    mMaxSpinBox->setValue(info.maxReal);
    mInterfaceUpdatePending = false;
  }

  const Real minValue = static_cast<Real>(mMinSpinBox->value());
  const Real maxValue = static_cast<Real>(mMaxSpinBox->value());

  Real aspectRatio = mBufferInfo.aspectRatio();
  // for an image in Fourier space the aspect ratio must be corrected
  if (mFourier != mFftCheckBox->isChecked())
    aspectRatio = cols / (rows * aspectRatio);

  QSize size = calculateSize(512, aspectRatio);
  mImageLabel->resize(size);
  mImageLabel->move(266 - size.width() / 2, 266 - size.height() / 2);

  // For the two phase representations the min and max values are ignored.
  // The range from -pi to pi is mapped either on all gray or on all color
  // values.
  if (mPhaseButton->isChecked())
  {
    auto extractImage = [] (const Complex& z)
    {
      const Real scaledValue = R(0.5) + z.real() / Math::twoPi;
      return Gray8(saturate(scaledValue * R(255.0)));
    };
    // FIXME: MSVC - remove cast
    mImage = Image::fromCBuffer2D(buffer, (Image::ExtractorCGray)extractImage);
  }
  else if (mPhaseColorButton->isChecked())
  {
    auto extractImage = [] (const Complex& z)
    {
      Vector3R hsv(z.real() + Math::pi, R(1.0), R(1.0));
      Vector3R rgb = hsv2Rgb(hsv);

      return Rgba8(saturate(rgb.x() * R(255.0)),
                   saturate(rgb.y() * R(255.0)),
                   saturate(rgb.z() * R(255.0)));
    };
    // FIXME: MSVC - remove cast
    mImage = Image::fromCBuffer2D(buffer, (Image::ExtractorCColor)extractImage);
  }
  else
  {
    auto extractImage = [minValue, maxValue] (const Complex& z)
    {
      const Real scaledValue = (z.real() - minValue) / (maxValue - minValue);
      return Gray8(saturate(scaledValue * R(255.0)));
    };
    // FIXME: MSVC - remove cast
    mImage = Image::fromCBuffer2D(buffer, (Image::ExtractorCGray)extractImage);
  }

  Qt::TransformationMode trafoMode = mFilterCheckBox->isChecked() ?
                                     Qt::SmoothTransformation :
                                     Qt::FastTransformation;
  auto qImage = imageToQImage(mImage).scaled(size, Qt::IgnoreAspectRatio,
                                             trafoMode);
  mImageLabel->setPixmap(QPixmap::fromImage(qImage));
  mImageLabel->setScaledContents(true);

  mProbability = mBuffer.absSqrReduce() / mBuffer.numPixels();
  drawHistogram();
}

void QView::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/plain"))
    event->acceptProposedAction();
}

void QView::dropEvent(QDropEvent* event)
{
  QUrl url = QUrl::fromEncoded(event->mimeData()->text().trimmed().toUtf8());
  if (!url.isLocalFile())
    return;
  QString filename = url.toLocalFile();

  event->accept();
  loadBuffer(filename);
}

#if 0
// fire colors
Vector3 fireColor(Real x)
{
  x = clamp(x, R(0.0), R(1.0));

  Vector3 red(1.0, 0.0, 0.0);
  Vector3 orange(1.0, 0.64, 0.0);
  Vector3 yellow(1.0, 1.0, 0.0);
  Vector3 white(1.0, 1.0, 1.0);

  if (x < 0.4)
  {
    x = x / 0.4;
    return x * red;
  }

  if (x < 0.5)
  {
    x = (x - 0.4) / 0.1;
    return x * orange + (1 - x) * red;
  }

  if (x < 0.6)
  {
    x = (x - 0.5) / 0.1;
    return x * yellow + (1 - x) * orange;
  }

  x = (x - 0.6) / 0.4;
  return x * white + (1 - x) * yellow;
}
#endif

} // namespace Aurora
