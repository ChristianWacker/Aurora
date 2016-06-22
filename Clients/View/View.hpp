//--- Aurora/Clients/View/View.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_VIEW_VIEW_HPP
#define AURORA_CLIENTS_VIEW_VIEW_HPP

#include "AuroraLib/FourierRingCorrelation.hpp"
#include "AuroraLib/Image.hpp"
#include "qwt_plot_curve.h"
#include "ui_MainView.h"

namespace Aurora
{

class View : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT

public:
  View(QMainWindow* parent = nullptr);
  ~View() { }

private:
  QString mFscSavePath;

  QwtPlotCurve mCurveBlack;
  QwtPlotCurve mCurveRed1;
  QwtPlotCurve mCurveRed2;
  QwtPlotCurve mCurveBlue1;
  QwtPlotCurve mCurveBlue2;

  std::unique_ptr<FourierRingCorrelation> mFsc;
  Real mMaxFrequency;

private slots:
  /// Saves the data from the Fourier Shell Correlation (FSC)
  void saveFscData();
  /// Saves the Fourier Shell Correlation (FSC) Plot
  void saveFscPlot();
  void showContextMenu(const QPoint& position);
  /// Update the Fourier Ring Correlation
  void updateFrc();
};

} // namespace Aurora

#endif
