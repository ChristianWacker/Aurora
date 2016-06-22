//--- Aurora/Clients/View/View.cpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "View.hpp"

#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/FourierRingCorrelation.hpp"
#include "AuroraLib/Utils.hpp"
#include "QtSupport/QtSupport.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <QApplication>
#include <QDir>
#include <QFileDialog>
#include <qwt_plot_canvas.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_renderer.h>
#include <qwt_text.h>
#include <vector>

namespace Aurora
{

CL::Option<std::string> inFilename1("i1", "Filename for view1");
CL::Option<std::string> inFilename2("i2", "Filename for view2");

View::View(QMainWindow* parent) :
  QMainWindow(parent),
  mFscSavePath(QDir::homePath())
{
  setupUi(this);

  connect(view1, &QView::bufferChanged, this, &View::updateFrc);
  connect(view2, &QView::bufferChanged, this, &View::updateFrc);

    void (QSpinBox::*intChanged)(int) = &QSpinBox::valueChanged;
  connect(frcPointsSpinBox, intChanged, this, &View::updateFrc);

  // make the plot background white including the label and suppress the border
  plot->setAutoFillBackground(true);
  plot->setPalette(Qt::white);
  plot->setCanvasBackground(Qt::white);
  auto canvas = static_cast<QwtPlotCanvas*>(plot->canvas());
  canvas->setFrameStyle(QFrame::NoFrame);

  // give the axes titles:
  QwtText xTitle("Spatial Frequency [1 / Angstrom]");
  xTitle.setFont(QFont("Times", 14));
  plot->setAxisTitle(QwtPlot::xBottom, xTitle);

  QwtText yTitle("Correlation");
  yTitle.setFont(QFont("Times", 14));
  plot->setAxisTitle(QwtPlot::yLeft, yTitle);

  // add an gray dotted grid
  QwtPlotGrid* grid = new QwtPlotGrid();
  grid->setPen(QPen(Qt::gray, 0, Qt::DotLine));
  grid->attach(plot);

  mCurveBlack.attach(plot);
  mCurveRed1.attach(plot);
  mCurveRed2.attach(plot);
  mCurveBlue1.attach(plot);
  mCurveBlue2.attach(plot);

  mCurveBlack.setPen(QPen(Qt::black));
  mCurveRed1.setPen(QPen(Qt::red));
  mCurveRed2.setPen(QPen(Qt::red));
  mCurveBlue1.setPen(QPen(Qt::blue));
  mCurveBlue2.setPen(QPen(Qt::blue));

  plot->setAxisScale(QwtPlot::yLeft, -1.05, 1.05);

  plot->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(plot, &QwtPlot::customContextMenuRequested,
          this, &View::showContextMenu);

  updateFrc();

  view1->loadBuffer(QString::fromStdString(inFilename1));
  view2->loadBuffer(QString::fromStdString(inFilename2));
}

void View::saveFscData()
{
  if (!mFsc)
    return;

  try
  {
    const QString caption = tr("Save Fourier Shell Correlation Data");
    const QString filter  = tr("text file (*.txt)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "txt", mFscSavePath);

    if (filename.empty())
      return;

    std::ofstream file(filename);
    for (int64_t i = 0, numBins = mFsc->numBins(); i < numBins; ++i)
      file << i * mMaxFrequency / numBins << ' ' << mFsc->bin(i) << '\n';
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void View::saveFscPlot()
{
  const QString caption = tr("Save Fourier Shell Correlation Plot");
  const QString filter  = tr("Image files (*.*)");
  const std::string filename =
    saveFileDialog(this, caption, filter, "pdf", mFscSavePath);

  if (filename.empty())
    return;

  QwtPlotRenderer renderer;
  renderer.renderDocument(plot, QString::fromStdString(filename),
                          QSizeF(200, 200), 300);
}

void View::updateFrc()
{
  const auto& buffer1 = view1->buffer();
  const auto& buffer2 = view2->buffer();

  // for the Fourier shell correlation the two buffers must be compatible
  if (!buffer1.data() || !buffer2.data() || !buffer1.compatible(buffer2))
  {
    // paint a red cross
    mCurveBlack.setVisible(false);
    mCurveRed1.setSamples(QVector<double>{ 0.0,  1.0},
                          QVector<double>{-1.0,  1.0});
    mCurveRed2.setSamples(QVector<double>{ 0.0,  1.0},
                          QVector<double>{ 1.0, -1.0});
    mCurveBlue1.setVisible(false);
    mCurveBlue2.setVisible(false);

    plot->setAxisScale(QwtPlot::xBottom, 0, 1.0);

    plot->replot();

    return;
  }

  QVector<double> dataX;
  QVector<double> dataY;

  auto complexBuffer1Abs = buffer1.bufferAbs();
  auto complexBuffer2Abs = buffer2.bufferAbs();

  Buffer2DInfo& bufferInfo = view1->bufferInfo();

  // the next frequencies are in 1 / Angstrom (EM notation)
  const Real kXNyquist = bufferInfo.kXNyquist() / Math::twoPi;
  const Real kYNyquist = bufferInfo.kYNyquist() / Math::twoPi;
  mMaxFrequency = std::hypot(kXNyquist, kYNyquist);

  // size of a bin in inverse Angstroms
  const Real binSize = mMaxFrequency / frcPointsSpinBox->value();

  using namespace std::chrono;
  const auto startTime = high_resolution_clock::now();

  mFsc = std::make_unique<FourierRingCorrelation>(complexBuffer1Abs,
                                                   complexBuffer2Abs,
                                                   bufferInfo.deltaKX(),
                                                   bufferInfo.deltaKY(),
                                                   binSize * Math::twoPi);

  for (int64_t i = 0, numBins = mFsc->numBins(); i < numBins; ++i)
  {
    dataX.push_back(i * mMaxFrequency / mFsc->numBins());
    dataY.push_back(mFsc->bin(i));
  }

  const auto stopTime = high_resolution_clock::now();
  std::cout << "Time needed for FSC: "
            << duration(stopTime - startTime) << '\n';

  // make graph visible and assign data to it:
  mCurveBlack.setVisible(true);
  mCurveBlack.setSamples(dataX, dataY);

  // y values for a vertical line
  QVector<double> verticalLine = {-1.0, 1.0};

  // mark the Nyquist frequencies
  mCurveRed1.setSamples(QVector<double>{kXNyquist, kXNyquist}, verticalLine);
  mCurveRed2.setSamples(QVector<double>{kYNyquist, kYNyquist}, verticalLine);

  // mark the Bandlimit
  mCurveBlue1.setVisible(true);
  mCurveBlue2.setVisible(true);
  const Real bandlimitX = R(2.0) / R(3.0) * kXNyquist;
  mCurveBlue1.setSamples(QVector<double>{bandlimitX, bandlimitX}, verticalLine);
  const Real bandlimitY = R(2.0) / R(3.0) * kYNyquist;
  mCurveBlue2.setSamples(QVector<double>{bandlimitY, bandlimitY}, verticalLine);

  // set axes ranges, so we can see all the data:
  plot->setAxisScale(QwtPlot::xBottom, 0, mMaxFrequency);

  plot->replot();
}

void View::showContextMenu(const QPoint& position)
{
  QPoint globalPosition = plot->mapToGlobal(position);

  QMenu menu;
  QAction savePlot("Save Plot...", this);
  QAction savePlotData("Save Data...", this);
  menu.addAction(&savePlot);
  menu.addAction(&savePlotData);
  connect(&savePlot, &QAction::triggered, this, &View::saveFscPlot);
  connect(&savePlotData, &QAction::triggered, this, &View::saveFscData);
  menu.exec(globalPosition);
}

} // namespace Aurora

int main(int argc, char* argv[])
{
  try
  {
    QLocale::setDefault(QLocale::c());
    QApplication a(argc, argv);
    setlocale(LC_ALL, "C");
    if (!Aurora::CL::parse(argc, argv, "Aurora View"))
      return 1;
    Aurora::View v;
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
