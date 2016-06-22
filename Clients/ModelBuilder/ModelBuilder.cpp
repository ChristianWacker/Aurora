//--- Aurora/Clients/ModelBuilder/ModelBuilder.cpp -----------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "ModelBuilder.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Crystal.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/RBuffer2D.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/Utils.hpp"
#include "QtSupport/QtSupport.hpp"

#include <cstdio>
#include <Eigen/Geometry>
#include <iostream>
#include <QFileDialog>
#include <thread>

namespace Aurora
{

CL::Option<bool> showProjPotential("showProjPotential", "", true);

ModelBuilder::ModelBuilder(QMainWindow* parent) :
  QMainWindow(parent),
  mSaveImagePath(QDir::homePath()),
  mSaveSamplePath(QDir::homePath())
{
  setupUi(this);

  for (auto i : Crystal::namesList())
    crystalNameComboBox->addItem(QString::fromStdString(i));

  auto buttonPushed = &QAbstractButton::clicked;
  connect(saveSampleButton, buttonPushed, this, &ModelBuilder::saveSample);
  connect(saveImageButton, buttonPushed, this, &ModelBuilder::saveImage);
  connect(viewButton, buttonPushed, this, &ModelBuilder::view);
  connect(updateButton, buttonPushed, this, &ModelBuilder::update);
  connect(sphereButton, buttonPushed, this, &ModelBuilder::createSphere);

  update();
}

void ModelBuilder::update()
{
  std::string crystalName = crystalNameComboBox->currentText().toStdString();
  auto crystal = Crystal::create(crystalName);

  // get the rotations matrix from the GUI
  Matrix3x3R rotation = Matrix3x3R::Identity();
  if (((directionHSpinBox->value() != 0) ||
       (directionKSpinBox->value() != 0) ||
       (directionLSpinBox->value() != 0)) &&
      ((rightHSpinBox->value() != 0) ||
       (rightKSpinBox->value() != 0) ||
       (rightLSpinBox->value() != 0)))
  {
    Vector3R direction =
      directionHSpinBox->value() * crystal->reciprocalLatticeVec1() +
      directionKSpinBox->value() * crystal->reciprocalLatticeVec2() +
      directionLSpinBox->value() * crystal->reciprocalLatticeVec3();
    direction.normalize();

    const Vector3R right =
      rightHSpinBox->value() * crystal->reciprocalLatticeVec1() +
      rightKSpinBox->value() * crystal->reciprocalLatticeVec2() +
      rightLSpinBox->value() * crystal->reciprocalLatticeVec3();

    // construct an orthogonal vector to direction
    Vector3R ortho1 = right - direction.dot(right) * direction;
    ortho1.normalize();

    const Vector3R ortho2 = direction.cross(ortho1);

    // this matrix rotates direction in the z-direction and ortho1 in the
    // x-direction
    rotation.row(0) = ortho1;
    rotation.row(1) = ortho2;
    rotation.row(2) = direction;
  }

  // get the translation
  const Vector3R translation(static_cast<Real>(translationXSpinBox->value()),
                             static_cast<Real>(translationYSpinBox->value()),
                             static_cast<Real>(translationZSpinBox->value()));

  const Vector3R minBox(static_cast<Real>(minBoxXSpinBox->value()),
                        static_cast<Real>(minBoxYSpinBox->value()),
                        static_cast<Real>(minBoxZSpinBox->value()));

  const Vector3R maxBox(static_cast<Real>(maxBoxXSpinBox->value()),
                        static_cast<Real>(maxBoxYSpinBox->value()),
                        static_cast<Real>(maxBoxZSpinBox->value()));

  mSample = crystal->generateSample(minBox, maxBox, rotation, translation);

  mSample->preprocess();

  Vector3R minBoxSample = mSample->minQuadraticBoundingBox(R(1.1));
  Vector3R maxBoxSample = mSample->maxQuadraticBoundingBox(R(1.1));

  if (showProjPotential)
  {
    const int64_t cols = projectedPotentialLabel->width();
    const int64_t rows = projectedPotentialLabel->height();

    auto projPotential =
      mSample->projectedPotential(cols, rows, minBoxSample, maxBoxSample, 3);

    auto info = projPotential.info();

    auto extract = [info] (Real x)
    {
      const Real value = (x - info.min) / (info.max - info.min);
      return Rgba8(saturate(value * R(255.0)));
    };

    mProjPotentialImage = Image::fromRBuffer2D(projPotential, extract);
    auto pixmap = QPixmap::fromImage(imageToQImage(mProjPotentialImage));
    projectedPotentialLabel->setPixmap(pixmap);
  }

  numAtomsLabel->setText(QString::number(mSample->numAtoms()));
}

void ModelBuilder::createSphere()
{
  Vector3R center(0.0, 0.0, 0.0);

  // calculate the center
  for (auto atom : mSample->atoms())
    center += atom.position();

  center /= mSample->numAtoms();

  Real radius = static_cast<Real>(radiusSpinBox->value());
  AtomVector newAtoms;

  for (auto atom : mSample->atoms())
    if ((atom.position() - center).norm() < radius)
      newAtoms.push_back(atom);

  mSample->atoms() = newAtoms;
}

void ModelBuilder::saveImage()
{
  try
  {
    const QString caption = tr("Save Image");
    const QString filter = tr("Images (*.tif *.tiff)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "tiff", mSaveImagePath);

    mProjPotentialImage.save(filename);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void ModelBuilder::saveSample()
{
  try
  {
    const QString caption = tr("Save Sample");
    const QString filter = tr("Samples (*.pdb *.atom)");
    const std::string filename =
      saveFileDialog(this, caption, filter, "pdb", mSaveSamplePath);

    if (!filename.empty())
      mSample->save(filename);
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
}

void spawnVmd(const std::string& pdbFilename)
{
  std::string vmdCommand = QDir::homePath().toStdString() + "/vmd/bin/vmd";

  if (-1 == system((vmdCommand + " " + pdbFilename).c_str()))
    std::cerr << "Calling VMD failed." << std::endl;

  if (0 != remove(pdbFilename.c_str()))
    std::cerr << "\"" << pdbFilename << "\" could not be remove." << std::endl;
}

void ModelBuilder::view()
{
  std::string tempFilename(tmpnam(nullptr));
  tempFilename += ".pdb";
  std::cout << tempFilename << std::endl;

  mSample->save(tempFilename);

  std::thread thread(spawnVmd, tempFilename);
  thread.detach();
}

} //namespace Aurora

int main(int argc, char* argv[])
{
  try
  {
    QLocale::setDefault(QLocale::c());
    QApplication a(argc, argv);
    setlocale(LC_ALL, "C");
    if (!Aurora::CL::parse(argc, argv))
      return EXIT_FAILURE;
    Aurora::ModelBuilder modelBuilder;
    modelBuilder.show();
    return a.exec();
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
