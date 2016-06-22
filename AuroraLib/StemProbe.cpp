//--- Aurora/AuroraLib/StemProbe.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/StemProbe.hpp"

#include "AuroraLib/Ctf.hpp"

namespace Aurora
{

StemProbe::StemProbe(const std::shared_ptr<CtfCoherent>& ctf,
                     const MultisliceParameters& params,
                     CBuffer2D& targetBuffer) :
  mCtf(ctf), mParams(params), mFourierBuffer(mParams.cols(), mParams.rows()),
  mFft(Fft::create(mFourierBuffer, targetBuffer, Direction::backward, false)),
  mScalingFactor(1.0)
{
  // The probe should normed. A test calculation is the simplest way to
  // calculate the norm.
  generate(0.0, 0.0);
  mScalingFactor = 1 / std::sqrt(targetBuffer.absSqrReduce());
}

void StemProbe::generate(Real centerX, Real centerY)
{
  // We are using the beam shift parameters of the CTF to move the probe.
  mCtf->setA(0, Complex(-centerX, -centerY));
  mFourierBuffer.assign(mCtf->buffer(mParams, mScalingFactor));
  mFft->transform();
}

} // namespace Aurora
