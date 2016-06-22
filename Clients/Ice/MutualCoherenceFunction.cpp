//--- Aurora/Clients/Ice/MutualCoherenceFunction.cpp ---------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "MutualCoherenceFunction.hpp"

#include "AuroraLib/FftFftw.hpp"

namespace Aurora
{

/*static*/ MutualCoherenceFunction MutualCoherenceFunction::zero(int64_t cols,
                                                                 int64_t rows)
{
  MutualCoherenceFunction result(cols, rows);
  std::fill(result.data(), result.data() + result.numElements(), 0);
  return result;
}

/*static*/ MutualCoherenceFunction MutualCoherenceFunction::source(
                                          const Buffer2DInfo& info, Real radius)
{
  const int64_t cols = info.cols();
  const int64_t rows = info.rows();
  (void) radius;
  //const Real radiusSqr = radius * radius;

  MutualCoherenceFunction result = zero(info.cols(), info.rows());

  #pragma omp parallel for
  for (int64_t y = 0; y < rows; ++y)
  {
    //const Real yPos = y * info.deltaY();
    for (int64_t x = 0; x < cols; ++x)
    {
      //const Real xPos = x * info.deltaX();
    //  const Real rSqr = xPos * xPos + yPos * yPos;
  //    if (rSqr < radiusSqr)
        result.element(x, y, x, y) = Complex(1.0);
    }
  }

  return result;
}

/*static*/ MutualCoherenceFunction MutualCoherenceFunction::pure(
                                                        const CBuffer2D& buffer)
{
  const int64_t cols = buffer.cols();
  const int64_t rows = buffer.rows();

  MutualCoherenceFunction result(cols, rows);

  #pragma omp parallel for
  for (int64_t y2 = 0; y2 < rows; ++y2)
  {
    for (int64_t x2 = 0; x2 < cols; ++x2)
    {
      for (int64_t y1 = 0; y1 < rows; ++y1)
      {
        for (int64_t x1 = 0; x1 < cols; ++x1)
        {
          result.element(x1, y1, x2, y2) =
            buffer.pixel(x1, y1) * std::conj(buffer.pixel(x2, y2));
        }
      }
    }
  }

  return result;
}

void MutualCoherenceFunction::init()
{
  Fftw::init();
  // The access to FFTW planner routines must be serialized
  std::unique_lock<std::mutex> lock(Fftw::getMutex());

  // todo: destroy plans

  // our data is stored in column-major format. Hence, the last variable, y2
  // varies most slowly.
  fftw_iodim64 firstDimensions[2];
  firstDimensions[0].n  = mRows;
  firstDimensions[0].is = mCols;
  firstDimensions[0].os = mCols;
  firstDimensions[1].n  = mCols;
  firstDimensions[1].is = 1;
  firstDimensions[1].os = 1;

  fftw_iodim64 secondDimensions[2];
  secondDimensions[0].n  = mRows;
  secondDimensions[0].is = mRows * mCols * mCols;
  secondDimensions[0].os = mRows * mCols * mCols;
  secondDimensions[1].n  = mCols;
  secondDimensions[1].is = mRows * mCols;
  secondDimensions[1].os = mRows * mCols;

  mFftPlan1 = fftw_plan_guru64_dft(2, firstDimensions, 2, secondDimensions,
    reinterpret_cast<fftw_complex*>(mData.get()),
    reinterpret_cast<fftw_complex*>(mData.get()),
    FFTW_BACKWARD, FFTW_ESTIMATE);

  mFftPlan2 = fftw_plan_guru64_dft(2, secondDimensions, 2, firstDimensions,
    reinterpret_cast<fftw_complex*>(mData.get()),
    reinterpret_cast<fftw_complex*>(mData.get()),
    FFTW_FORWARD, FFTW_ESTIMATE);

  mFftInvPlan1 = fftw_plan_guru64_dft(2, firstDimensions, 2, secondDimensions,
    reinterpret_cast<fftw_complex*>(mData.get()),
    reinterpret_cast<fftw_complex*>(mData.get()),
    FFTW_FORWARD, FFTW_ESTIMATE);

  mFftInvPlan2 = fftw_plan_guru64_dft(2, secondDimensions, 2, firstDimensions,
    reinterpret_cast<fftw_complex*>(mData.get()),
    reinterpret_cast<fftw_complex*>(mData.get()),
    FFTW_BACKWARD, FFTW_ESTIMATE);
}

void MutualCoherenceFunction::doFft(Direction direction)
{
  if (Direction::forward == direction)
  {
    fftw_execute(mFftPlan1);
    fftw_execute(mFftPlan2);
  }
  else
  {
    fftw_execute(mFftInvPlan1);
    fftw_execute(mFftInvPlan2);
  }
}

void MutualCoherenceFunction::transform(Direction direction)
{
  // Fourier transformation in x1 and y1
  doFft(direction);
}

void MutualCoherenceFunction::apply(const CBuffer2D& buffer)
{
  assert((buffer.cols()== mCols) && (buffer.rows() == mRows) &&
         "Buffer is incompatible.");

  #pragma omp parallel for
  for (int64_t y2 = 0; y2 < mRows; ++y2)
  {
    for (int64_t x2 = 0; x2 < mCols; ++x2)
    {
      for (int64_t y1 = 0; y1 < mRows; ++y1)
      {
        for (int64_t x1 = 0; x1 < mCols; ++x1)
        {
          element(x1, y1, x2, y2) *=
            buffer.pixel(x1, y1) * std::conj(buffer.pixel(x2, y2));
        }
      }
    }
  }
}

CBuffer2D MutualCoherenceFunction::slice(int64_t y) const
{
  CBuffer2D result(mCols, mCols);
  if (-1 == y)
    y = mRows / 2;

  #pragma omp parallel for
  for (int64_t x1 = 0; x1 < mCols; ++x1)
    for (int64_t x2 = 0; x2 < mCols; ++x2)
      result.pixel(x1, x2) = element(x1, y, x2, y);

  return result;
}

} // namespace Aurora
