//--- Aurora/AuroraLib/FftKiss.cpp ---------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
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

#include "AuroraLib/FftKiss.hpp"

#include "AuroraLib/Complex.hpp"
#include "AuroraLib/CBuffer2D.hpp"
#include "AuroraLib/Math.hpp"

namespace Aurora
{

namespace KissFft
{

OneDim::OneDim(size_t numElements, Direction direction) :
  mNumElements(numElements),
  mDirection(direction),
  mTwiddles(numElements)
{
  for (size_t i = 0; i < numElements; ++i)
  {
    Real phase = -Math::twoPi * i / numElements;
    if (Direction::forward == direction)
      phase *= -1;

    mTwiddles[i] = fromArg(phase);
  }

  /// mFactors is populated by p1, m1, p2, m2, ...
  /// where
  /// p[i] * m[i] = m[i - 1]
  /// m0 = n
  size_t p = 4;
  size_t n = numElements;
  const Real floorSqrt = std::floor(std::sqrt(static_cast<Real>(n)));

  // factor out powers of 4, powers of 2, then any remaining primes
  do
  {
    while (n % p)
    {
      switch (p)
      {
      case 2:
        p = 3;
        break;
      case 4:
        p = 2;
        break;
      default:
        p += 2;
        break;
      }

      if (p > floorSqrt)
        p = n; // no more factors, skip to end
    }

    n /= p;
    mFactors.push_back(p);
    mFactors.push_back(n);
  }
  while (n > 1);
}

void OneDim::butterfly2(Complex* out, size_t stride, size_t m)
{
  Complex* out2 = out + m;
  const Complex* tw1 = mTwiddles.data();

  do
  {
    Complex t = *out2 * *tw1;
    tw1 += stride;
    *out2 = *out - t;
    *out += t;
    ++out2;
    ++out;
  }
  while (--m);
}

void OneDim::butterfly3(Complex* out, size_t stride, size_t m)
{
  size_t k = m;
  const size_t m2 = 2 * m;
  const Complex* tw1 = mTwiddles.data();
  const Complex* tw2 = mTwiddles.data();
  Complex epi3 = tw1[stride * m];

  do
  {
    Complex scratch1 = out[m] * *tw1;
    Complex scratch2 = out[m2] * *tw2;

    Complex scratch3 = scratch1 + scratch2;
    Complex scratch0 = scratch1 - scratch2;

    tw1 += stride;
    tw2 += stride * 2;

    out[m] = out[0] - R(0.5) * scratch3;

    scratch0 *= epi3.imag();

    out[0] += scratch3;

    out[m2].real(out[m].real() + scratch0.imag());
    out[m2].imag(out[m].imag() - scratch0.real());

    out[m ].real(out[m].real() - scratch0.imag());
    out[m ].imag(out[m].imag() + scratch0.real());

    ++out;
  }
  while (--k);
}

void OneDim::butterfly4(Complex* out, size_t stride, size_t m)
{
  const Complex* tw1 = mTwiddles.data();
  const Complex* tw2 = mTwiddles.data();
  const Complex* tw3 = mTwiddles.data();
  size_t k = m;
  const size_t m2 = 2 * m;
  const size_t m3 = 3 * m;

  do
  {
    Complex scratch0 = out[m] * *tw1;
    Complex scratch1 = out[m2] * *tw2;
    Complex scratch2 = out[m3] * *tw3;

    Complex scratch5 = *out - scratch1;
    *out += scratch1;
    Complex scratch3 = scratch0 + scratch2;
    Complex scratch4 = scratch0 - scratch2;
    out[m2] = *out - scratch3;
    tw1 += stride;
    tw2 += stride * 2;
    tw3 += stride * 3;
    *out += scratch3;

    if (Direction::forward == mDirection)
    {
      out[m ].real(scratch5.real() - scratch4.imag());
      out[m ].imag(scratch5.imag() + scratch4.real());
      out[m3].real(scratch5.real() + scratch4.imag());
      out[m3].imag(scratch5.imag() - scratch4.real());
    }
    else
    {
      out[m ].real(scratch5.real() + scratch4.imag());
      out[m ].imag(scratch5.imag() - scratch4.real());
      out[m3].real(scratch5.real() - scratch4.imag());
      out[m3].imag(scratch5.imag() + scratch4.real());
    }

    ++out;
  }
  while (--k);
}

void OneDim::butterfly5(Complex* out, size_t stride, size_t m)
{
  Complex* out0 = out;
  Complex* out1 = out + m;
  Complex* out2 = out + 2 * m;
  Complex* out3 = out + 3 * m;
  Complex* out4 = out + 4 * m;
  const Complex* tw(mTwiddles.data());
  Complex ya = mTwiddles[stride * m];
  Complex yb = mTwiddles[stride * 2 * m];

  for (size_t u = 0; u < m; ++u)
  {
    Complex scratch0  = *out0;

    Complex scratch1  = *out1 * tw[    u * stride];
    Complex scratch2  = *out2 * tw[2 * u * stride];
    Complex scratch3  = *out3 * tw[3 * u * stride];
    Complex scratch4  = *out4 * tw[4 * u * stride];

    Complex scratch7  = scratch1 + scratch4;
    Complex scratch10 = scratch1 - scratch4;
    Complex scratch8  = scratch2 + scratch3;
    Complex scratch9  = scratch2 - scratch3;

    *out0++ += scratch7 + scratch8;

    Complex scratch5 = scratch0 + scratch7 * ya.real() + scratch8 * yb.real();
    Complex scratch6( scratch10.imag() * ya.imag() + scratch9.imag() * yb.imag(),
                     -scratch10.real() * ya.imag() - scratch9.real() * yb.imag());

    *out1++ = scratch5 - scratch6;
    *out4++ = scratch5 + scratch6;

    Complex scratch11 = scratch0 + scratch7 * yb.real() + scratch8 * ya.real();
    Complex scratch12(-scratch10.imag() * yb.imag() + scratch9.imag() * ya.imag(),
                       scratch10.real() * yb.imag() - scratch9.real() * ya.imag());

    *out2++ = scratch11 + scratch12;
    *out3++ = scratch11 - scratch12;
  }
}

// Performs the butterfly for one stage of a mixed radix FFT.
void OneDim::butterflyGeneric(Complex* out, size_t stride, size_t m,
                              size_t radix)
{
  std::vector<Complex> scratch(radix);

  for (size_t u = 0; u < m; ++u)
  {
    size_t k = u;
    for (size_t q1 = 0; q1 < radix; ++q1)
    {
      scratch[q1] = out[k];
      k += m;
    }

    k = u;
    for (size_t q1 = 0; q1 < radix; ++q1)
    {
      size_t twidx = 0;
      out[k] = scratch[0];
      for (size_t q = 1; q < radix; ++q)
      {
        twidx += stride * k;
        if (twidx >= mNumElements)
          twidx -= mNumElements;

        out[k] += scratch[q] * mTwiddles[twidx];
      }
      k += m;
    }
  }
}

void OneDim::doTransform(Complex* out, const Complex* f, size_t stride,
                         const size_t* factors, size_t inStride)
{
  Complex* outBegin = out;
  const size_t radix = *factors++;
  const size_t m = *factors++; // stage's FFT length/radix
  const Complex* outEnd = out + radix * m;

  if (1 == m)
  {
    do
    {
      *out = *f;
      f += stride * inStride;
    }
    while (++out != outEnd);
  }
  else
  {
    do
    {
      // recursive call:
      // DFT of size m * p performed by doing p instances of smaller DFTs of
      // size m, each one takes a decimated version of the input.
      doTransform(out, f, stride * radix, factors, inStride);
      f += stride * inStride;
    }
    while ((out += m) != outEnd);
  }

  out = outBegin;

  // recombine the p smaller DFTs
  switch (radix)
  {
  case 2:
    butterfly2(out, stride, m);
    break;
  case 3:
    butterfly3(out, stride, m);
    break;
  case 4:
    butterfly4(out, stride, m);
    break;
  case 5:
    butterfly5(out, stride, m);
    break;
  default:
    butterflyGeneric(out, stride, m, radix);
  }
}

void OneDim::transform(const Complex* in, Complex* out, size_t inStride)
{
  assert(in != out && "Kiss FFT does not support in-place transformations");
  doTransform(out, in, 1, mFactors.data(), inStride);
}

MultiDim::MultiDim(const std::vector<size_t>& dimensions, Direction direction) :
  mDimensions(dimensions),
  mTempBuffer(calcNumElements(dimensions)),
  mFfts(dimensions.size())
{
  for (size_t i = 0, numDimensions = dimensions.size(); i < numDimensions; ++i)
    mFfts[i] = std::make_unique<OneDim>(dimensions[i], direction);
}

size_t MultiDim::calcNumElements(const std::vector<size_t>& dimensions)
{
  size_t numElements = 1;

  for (size_t dimension : dimensions)
    numElements *= dimension;

  return numElements;
}

//  This works by tackling one dimension at a time.
//
//  In effect,
//  Each stage starts out by reshaping the matrix into a DixSi 2d matrix.
//  A Di-sized FFT is taken of each column, transposing the matrix as it goes.
//
// Here's a 3-d example:
// Take a 2x3x4 matrix, laid out in memory as a contiguous buffer
//  [ [ [ a b c d ] [ e f g h ] [ i j k l ] ]
//    [ [ m n o p ] [ q r s t ] [ u v w x ] ] ]
//
// Stage 0 ( D=2): treat the buffer as a 2x12 matrix
//    [ [a b ... k l]
//      [m n ... w x] ]
//
//    FFT each column with size 2.
//    Transpose the matrix at the same time using fft_stride.
//
//    [ [ a+m a-m ]
//      [ b+n b-n]
//      ...
//      [ k+w k-w ]
//      [ l+x l-x ] ]
//
//    Note fft([x y]) == [x+y x-y]
//
// Stage 1 (D = 3) treats the buffer (the output of stage D=2) as an 3x8 matrix,
//    [ [ a+m a-m b+n b-n c+o c-o d+p d-p ]
//      [ e+q e-q f+r f-r g+s g-s h+t h-t ]
//      [ i+u i-u j+v j-v k+w k-w l+x l-x ] ]
//
//    And perform FFTs (size=3) on each of the columns as above, transposing
//    the matrix as it goes.  The output of stage 1 is
//        (Legend: ap = [ a+m e+q i+u ]
//                 am = [ a-m e-q i-u ] )
//
//    [ [ sum(ap) fft(ap)[0] fft(ap)[1] ]
//      [ sum(am) fft(am)[0] fft(am)[1] ]
//      [ sum(bp) fft(bp)[0] fft(bp)[1] ]
//      [ sum(bm) fft(bm)[0] fft(bm)[1] ]
//      [ sum(cp) fft(cp)[0] fft(cp)[1] ]
//      [ sum(cm) fft(cm)[0] fft(cm)[1] ]
//      [ sum(dp) fft(dp)[0] fft(dp)[1] ]
//      [ sum(dm) fft(dm)[0] fft(dm)[1] ]  ]
//
// Stage 2 (D=4) treats this buffer as a 4*6 matrix,
//    [ [ sum(ap) fft(ap)[0] fft(ap)[1] sum(am) fft(am)[0] fft(am)[1] ]
//      [ sum(bp) fft(bp)[0] fft(bp)[1] sum(bm) fft(bm)[0] fft(bm)[1] ]
//      [ sum(cp) fft(cp)[0] fft(cp)[1] sum(cm) fft(cm)[0] fft(cm)[1] ]
//      [ sum(dp) fft(dp)[0] fft(dp)[1] sum(dm) fft(dm)[0] fft(dm)[1] ]  ]
//
//    Then FFTs each column, transposing as it goes.
//
//    The resulting matrix is the 3d FFT of the 2x3x4 input matrix.
//
//    Note as a sanity check that the first element of the final
//    stage's output (DC term) is
//    sum( [ sum(ap) sum(bp) sum(cp) sum(dp) ] )
//    , i.e. the summation of all 24 input elements.
void MultiDim::transform(const Complex* in, Complex* out)
{
  assert(in != out && "Kiss FFT does not support in-place transformations");

  const Complex* bufferIn = in;
  Complex* bufferOut;

  // arrange it so the last bufferOut == out
  if (mDimensions.size() & 1)
  {
    bufferOut = out;
    if (in == out)
    {
      std::copy(in, in + mTempBuffer.size(), mTempBuffer.data());
      bufferIn = mTempBuffer.data();
    }
  }
  else
  {
    bufferOut = mTempBuffer.data();
  }

  for (size_t k = 0, numDimensions = mDimensions.size(); k < numDimensions; ++k)
  {
    const int64_t stride = mTempBuffer.size() / mDimensions[k];

    #pragma omp parallel for
    for (int64_t i = 0; i < stride; ++i)
      mFfts[k]->transform(bufferIn + i, bufferOut + i * mDimensions[k], stride);

    // toggle back and forth between the two buffers
    if (mTempBuffer.data() == bufferOut)
    {
      bufferIn = mTempBuffer.data();
      bufferOut = out;
    }
    else
    {
      bufferIn = out;
      bufferOut = mTempBuffer.data();
    }
  }
}

} // namespace KissFft

FftKiss::FftKiss(const CBuffer2D& bufferIn, CBuffer2D& bufferOut,
                 Direction direction) :
  Fft(bufferIn, bufferOut, direction)
{
  assert(bufferIn.compatible(bufferOut));

  std::vector<size_t> dimensions;
  dimensions.push_back(mInBuffer.rows());
  dimensions.push_back(mInBuffer.cols());

  mMultiDim = std::make_unique<KissFft::MultiDim>(dimensions, direction);
}

void FftKiss::transform()
{
  // call the FFT algorithm
  mMultiDim->transform(mInBuffer.data(), mOutBuffer.data());
}

} // namespace Aurora
