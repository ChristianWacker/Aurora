//--- Aurora/Clients/ModelBuilder/ModelBuilder.hpp -----------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENT_MODEL_BUILDER_MODEL_BUILDER_HPP
#define AURORA_CLIENT_MODEL_BUILDER_MODEL_BUILDER_HPP

#include "AuroraLib/Image.hpp"
#include "AuroraLib/Sample.hpp"
#include "ui_ModelBuilderMain.h"

#include <complex>
#include <vector>

namespace Aurora
{

class ModelBuilder : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT

public:
  ModelBuilder(QMainWindow* parent = nullptr);
  ~ModelBuilder() { }

private:
  QString mSaveImagePath;
  QString mSaveSamplePath;

  Image mProjPotentialImage;
  SamplePtr mSample;

private slots:
  void update();
  void createSphere();
  void saveImage();
  void saveSample();
  void view();
};

} // namespace Aurora

#endif
