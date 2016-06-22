//--- Aurora/AuroraLib/FftKiss.hpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------
//
/// @file
/// @brief This FFT implementation is based on Kiss FFT ("Keep It Simple,
///  Stupid") by Mark Borgerding.
///  Website: http://sourceforge.net/projects/kissfft/.
//
//------------------------------------------------------------------------------

// Copyright (c) 2003-2010, Mark Borgerding
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//  * Neither the author nor the names of any contributors may be used to
//    endorse or promote products derived from this software without specific
//    prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef AURORA_AURORA_LIB_FFT_KISS_HPP
#define AURORA_AURORA_LIB_FFT_KISS_HPP

#include "AuroraLib/Complex.hpp"
#include "AuroraLib/Fft.hpp"

#include <vector>

namespace Aurora
{

namespace KissFft
{

class OneDim
{
public:
  OneDim(size_t numElements, Direction direction);
  void transform(const Complex* in, Complex* out, size_t inStride);

private:
  void butterfly2(Complex* out, size_t stride, size_t m);
  void butterfly3(Complex* out, size_t stride, size_t m);
  void butterfly4(Complex* out, size_t stride, size_t m);
  void butterfly5(Complex* out, size_t stride, size_t m);
  void butterflyGeneric(Complex* out, size_t stride, size_t m, size_t radix);
  void doTransform(Complex* out, const Complex* f, size_t stride,
                   const size_t* factors, size_t inStride);

  const size_t mNumElements;
  const Direction mDirection;
  std::vector<Complex> mTwiddles;
  std::vector<size_t> mFactors;
};

/// This class provides a FFT in arbitrary number of dimensions.
class MultiDim
{
public:
  MultiDim(const std::vector<size_t>& dimensions, Direction direction);
  void transform(const Complex* in, Complex* out);

private:
  size_t calcNumElements(const std::vector<size_t>& dimensions);

  std::vector<size_t> mDimensions;
  // buffer capable of holding the entire input
  std::vector<Complex> mTempBuffer;
  std::vector<std::unique_ptr<OneDim>> mFfts;
};

} // namespace KissFft

/// This is actually a wrapper class around @ref KissFft::MultiDim. It provides
/// the @ref Fft interface.
class AURORA_API FftKiss : public Fft
{
public:
  FftKiss(const CBuffer2D& bufferIn, CBuffer2D& bufferOut, Direction direction);

  void transform() override;

private:
  std::unique_ptr<KissFft::MultiDim> mMultiDim;
};

} // namespace Aurora

#endif
