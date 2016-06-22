//--- Aurora/Clients/OpSys2/OpSys2.hpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_OP_SYS2_OP_SYS2_HPP
#define AURORA_CLIENTS_OP_SYS2_OP_SYS2_HPP

#include "AuroraLib/Ctf.hpp"
#include "AuroraLib/CtfIncoherent.hpp"
#include "AuroraLib/Image.hpp"
#include "ui_OpSys2Main.h"

#include <qwt_plot_curve.h>

#include <QComboBox>

namespace Aurora
{

class OpSys2 : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT

public:
  OpSys2(QMainWindow* parent = nullptr);
  ~OpSys2() { }

private:
  bool mInterfaceUpdatePending;

  CtfIncoherent mCtf;

  QwtPlotCurve mCurveCtfReal;
  QwtPlotCurve mCurveCtfImag;
  QwtPlotCurve mCurveEnvelopePos;
  QwtPlotCurve mCurveEnvelopeNeg;
  QwtPlotCurve mCurveNumCtfReal;
  QwtPlotCurve mCurveNumCtfImag;

  QString mLoadPath;
  CBuffer2D mExitWave;
  Buffer2DInfo mExitWaveInfo;
  std::string mDescription;

  struct Constant
  {
    QDoubleSpinBox* abs;
    QDoubleSpinBox* angle;
    QComboBox* unitBox;
  };
  std::map<QString, Constant> mConstants;

  std::map<QString, Real> mUnits;

  void setScherzerDefocus();
  void addCoefficient(const QString& name, const QString& defaultUnit,
                      bool withAngles);
  Complex getComplexCoefficient(const QString& name);
  Real getRealCoefficient(const QString& name);
  void setRealCoefficient(const QString& name, Real value);

  void refreshView();

  void setBuffer(const CBuffer2D& buffer, const Buffer2DInfo& info,
                    const std::string& description);
  void loadExitWave(const QString& filename);
  void showExitWaveLoadDialog();

  void dragEnterEvent(QDragEnterEvent *event) override;
  void dropEvent(QDropEvent* event) override;
};

} // namespace Aurora

#endif
