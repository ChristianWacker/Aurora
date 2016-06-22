//--- Aurora/AuroraLib/SampleCodecAtom.cpp -------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/SampleCodecAtom.hpp"

#include "AuroraLib/BinaryIStream.hpp"
#include "AuroraLib/BinaryOStream.hpp"
#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Sample.hpp"
#include "AuroraLib/Utils.hpp"

#include <iomanip>
#include <iostream>

namespace Aurora
{

namespace
{
  struct AtomFileHeader
  {
    uint32_t version;
    uint64_t numberOfAtoms;
  };

  const uint32_t atomFileVersion = 0;
}

void SampleCodecAtom::load(const std::string& filename, Sample& sample)
{
  AtomVector& atoms = sample.atoms();
  atoms.clear();

  BinaryIStream stream = BinaryIStream::openFile(filename);

  if (stream.size() <= 16)
    AURORA_THROW(EInOut, "Cannot load \"" + filename + "\"");

  if (!stream.checkSignature("ATOM"))
    AURORA_THROW(EInOut, "Cannot load \"" + filename + "\"");

  const uint32_t version = stream.read<uint32_t>();
  // Check the version number
  if (atomFileVersion != version)
    AURORA_THROW(EInOut, "Cannot load \"" + filename + "\"");

  uint64_t numAtoms = stream.read<int64_t>();
  std::cout << "File contains " << numAtoms << '\n';
  while(numAtoms--)
  {
    Vector3R atomPosition;
    atomPosition.x() = stream.read<float>();
    atomPosition.y() = stream.read<float>();
    atomPosition.z() = stream.read<float>();

    const size_t atomicNumber = stream.read<uint16_t>();
    const int charge = stream.read<int16_t>();

    atoms.push_back(Atom(atomPosition, atomicNumber, charge));
  }

  sample.preprocess();
}

void SampleCodecAtom::save(const Sample& sample, const std::string& filename)
{
  BinaryOStream stream = BinaryOStream::createFile(filename);

  /// Write a simple header
  stream.write<char>("ATOM", 4);
  stream.write<uint32_t>(atomFileVersion);
  stream.write<uint64_t>(sample.numAtoms());

  for (auto atom : sample.atoms())
  {
    const Eigen::Vector3f atomPosition = atom.position().cast<float>();

    stream.write<float>(atomPosition.x());
    stream.write<float>(atomPosition.y());
    stream.write<float>(atomPosition.z());

    stream.write<uint16_t>(atom.atomicNumber());
    stream.write<int16_t>(atom.charge());
  }
}

class RegisterSampleCodecAtom
{
public:
  RegisterSampleCodecAtom()
  {
    auto codec = std::make_shared<SampleCodecAtom>();

    Sample::registerLoader("atom", codec);
    Sample::registerSaver("atom", codec);
  }
} registerSampleCodecAtom;

} // namespace Aurora
