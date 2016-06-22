//--- Aurora/AuroraLib/SampleCodecPdb.cpp --------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/SampleCodecPdb.hpp"

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/MemoryBlock.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/StringRef.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

namespace Aurora
{

CL::Option<bool> flipZ("flipZ", "", false);

void SampleCodecPdb::load(const std::string& filename, Sample& sample)
{
  AtomVector& atoms = sample.atoms();
  atoms.clear();

  MemoryBlock buffer = MemoryBlock::openFile(filename);

  if ((verbosity >= 1) && (flipZ == 1))
    std::cout << "Flipping z Coordinate\n";

  StringRef s(buffer.data(), buffer.size());

  while (!s.empty())
  {
    size_t pos = s.find('\n');

    StringRef line = s.subStr(0, pos);
    // remove the line
    s = s.frontTruncated(pos + 1);

    if (line.empty())
      continue; // next line

    StringRef tag = line.subStr(0, 6);
    if (("ATOM  " != tag) && ("HETATM" != tag))
      continue; // next line

    Real x = line.subStr(30, 8).toReal();
    Real y = line.subStr(38, 8).toReal();
    Real z = line.subStr(46, 8).toReal();

    const std::string elementSymbol = line.subStr(76, 2).trimmed().toUpper();
    const size_t atomicNumber = Element::atomicNumber(elementSymbol);

    StringRef chargeString = line.subStr(78, 2).trimmed();
    int32_t charge = 0;

    if (2 == chargeString.length())
    {
      charge = chargeString[0] - '0';
      if ('-' == chargeString[1])
        charge *= -1;
    }

    if (flipZ)
      z = -z;

    atoms.push_back(Atom(Vector3R(x, y, z), atomicNumber, charge));
  }

  sample.preprocess();
}

void SampleCodecPdb::save(const Sample& sample, const std::string& filename)
{
  std::ofstream outFile(filename);
  outFile.fill(' ');
  outFile.setf(std::ios::fixed);

  int64_t serialNumber = 0;

  for (auto atom : sample.atoms())
  {
    ++serialNumber;

    const std::string elementSymbol = Element::symbol(atom.atomicNumber());
    std::string atomName = elementSymbol;

    if (1 == atomName.length())
      atomName = " " + atomName;

    const Vector3R atomPosition = atom.position();

    const Real occupancy  = R(1.0);
    const Real tempFactor = R(0.0);

    const std::string charge = "  ";

    std::string serialString;
    if (serialNumber > 99999)
      serialString = "XXXXX";
    else
      serialString = std::to_string(serialNumber);

    outFile << "ATOM  " << std::right << std::setw(5) << serialString << ' '
            << std::left << std::setw(4) << atomName << "         1    "
            << std::setprecision(3) << std::right
            << std::setw(8) << atomPosition.x()
            << std::setw(8) << atomPosition.y()
            << std::setw(8) << atomPosition.z()
            << std::setprecision(2)
            << std::setw(6) << occupancy
            << std::setw(6) << tempFactor
            << std::left <<  "          "
            << std::setw(2) << std::right << elementSymbol
            << std::setw(2) << charge << '\n';
  }
}

class RegisterSampleCodecPdb
{
public:
  RegisterSampleCodecPdb()
  {
    auto codec = std::make_shared<SampleCodecPdb>();

    Sample::registerLoader("pdb", codec);
    Sample::registerSaver("pdb", codec);
  }
} registerSampleCodecPdb;

} // namespace Aurora
