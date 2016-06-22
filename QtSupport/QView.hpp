//--- Aurora/QtSupport/QView.hpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_QT_SUPPORT_Q_VIEW_HPP
#define AURORA_QT_SUPPORT_Q_VIEW_HPP

#include "AuroraLib/Buffer2DInfo.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Histogram.hpp"
#include "AuroraLib/Image.hpp"

#include "QtSupport/QtSupport.hpp"

#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QDropEvent>
#include <QLabel>
#include <QRadioButton>
#include <QMimeData>
#include <QWidget>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_textlabel.h>

namespace Aurora
{

class AURORA_QT_SUPPORT_API QView : public QWidget
{
  Q_OBJECT

public:
  QView(QWidget* parent, bool intern);

  const CBuffer2D& buffer() const
  {
    return mBuffer;
  }

  /// Sets a new ComplexBuffer that should be shown by the widget.
  /// @param fourier
  ///  This must be true, if the buffer is in Fourier space. It is assumed that
  ///  the highest frequencies are in the center of the buffer.
  void setBuffer(CBuffer2D&& buffer,
                 const Buffer2DInfo& bufferInfo,
                 const std::string& description,
                 bool fourier = false)
  {
    mBuffer = std::move(buffer);
    mBufferInfo = bufferInfo;
    mDescription = description;
    mFourier = fourier;

    updateView();

    emit(bufferChanged());
  }

  std::string description()
  {
    return mDescription;
  }

  Buffer2DInfo& bufferInfo()
  {
    return mBufferInfo;
  }

  void loadBuffer(const QString& filename);

  void dragEnterEvent(QDragEnterEvent *event) override;
  void dropEvent(QDropEvent* event) override;

signals:
  void bufferChanged();

private:
  bool mIntern;
  bool mInterfaceUpdatePending;
  Real mAspectRatio;
  bool mFourier;
  QString mLoadPath;
  QString mSavePath;
  QString mHistogramSavePath;
  Real mProbability;

  CBuffer2D mBuffer;
  Buffer2DInfo mBufferInfo;
  CBuffer2D mProcessedBuffer;
  std::string mDescription;
  Image mImage;
  std::unique_ptr<Histogram> mHistogram;
  QLabel* mImageLabel;
  QCheckBox* mAutoContrastCheckbox;
  QRadioButton* mRealButton;
  QRadioButton* mImagButton;
  QRadioButton* mAbsButton;
  QRadioButton* mAbsSqrButton;
  QRadioButton* mPhaseButton;
  QRadioButton* mPhaseColorButton;
  QRadioButton* mNoneRadioButton;
  QRadioButton* mAbsRadioButton;
  QRadioButton* mAbsSqrRadioButton;
  QCheckBox* mFftCheckBox;
  QCheckBox* mLogCheckBox;
  QCheckBox* mLog1pCheckBox;
  QCheckBox* mFilterCheckBox;
  QwtPlot* mHistogramPlot;
  QwtPlotCurve* mHistogramCurve;
  QwtPlotCurve* mHistogramMode;
  QwtPlotCurve* mHistogramMean;
  QwtPlotCurve* mHistogramSigma;
  QwtPlotTextLabel* mHistogramText;
  QDoubleSpinBox* mMinSpinBox;
  QDoubleSpinBox* mMaxSpinBox;

  void drawHistogram();

  void updateView();

private slots:
  void loadImage();
  void showBufferLoadDialog();
  void saveBufferAsText();
  void saveImage();
  void saveBuffer();

  void showBufferContextMenu(const QPoint& position);
  void showHistogramContextMenu(const QPoint& position);
  void saveHistogramImage();
  void saveHistogramData();
};

class AURORA_QT_SUPPORT_API QViewFile : public QView
{
  Q_OBJECT

public:
  QViewFile(QWidget* parent = nullptr) : QView(parent, false) { }
};

class AURORA_QT_SUPPORT_API QViewIntern : public QView
{
  Q_OBJECT

public:
  QViewIntern(QWidget* parent = nullptr) : QView(parent, true) { }
};

} // namespace Aurora

#endif
