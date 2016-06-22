//--- Aurora/Clients/Ice/MutualCoherenceFunction.hpp ---------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_ICE_MUTUAL_COHERENCE_FUNCTION_HPP
#define AURORA_CLIENTS_ICE_MUTUAL_COHERENCE_FUNCTION_HPP

#include "AuroraLib/FftFftw.hpp"
#include "AuroraLib/PhaseShift.hpp"

namespace Aurora
{

/// This class represents a mutual coherence function (MCF). It is a
/// four-dimensional object.
/// This class is movable, but not copyable.
class MutualCoherenceFunction
{
public:
  MutualCoherenceFunction() :
    mCols(0), mRows(0), mData(nullptr, Memory::alignedDelete)
  {}

  MutualCoherenceFunction(int64_t cols, int64_t rows) :
    mCols(cols), mRows(rows),
    mData(Memory::alignedNew<Complex>(numElements()), Memory::alignedDelete)
  {
    assert((cols > 0) && (rows > 0));
    init();
  }

  MutualCoherenceFunction(const MutualCoherenceFunction&) = delete;
  MutualCoherenceFunction(MutualCoherenceFunction&&) = default;
  MutualCoherenceFunction& operator=(const MutualCoherenceFunction&) = delete;
  MutualCoherenceFunction& operator=(MutualCoherenceFunction&&) = default;

  /// Return a MCF with all entries set to zero.
  static MutualCoherenceFunction zero(int64_t cols, int64_t rows);

  /// Return a MCF for an incoherent source.
  static MutualCoherenceFunction source(const Buffer2DInfo& info, Real radius);

  static MutualCoherenceFunction pure(const CBuffer2D& buffer);

  void init();

  void transform(Direction direction);

  // I/O
  //void load(const std::string& filename);
  //void save(const std::string& filename);

  /// Return a pointer to the internal data structure
  Complex* data()
  {
    return mData.get();
  }

  /// Return a pointer to the internal data structure
  const Complex* data() const
  {
    return mData.get();
  }

  Complex element(int64_t x1, int64_t y1, int64_t x2, int64_t y2) const
  {
    return mData[(((y2 * mCols) + x2) * mRows + y1) * mCols + x1];
  }

  Complex& element(int64_t x1, int64_t y1, int64_t x2, int64_t y2)
  {
    return mData[(((y2 * mCols) + x2) * mRows + y1) * mCols + x1];
  }

  /// Return the elements on the diagonal (x1 = x2, y1 = y2). In real space
  /// this is the image intensity.
  RBuffer2D diagonal() const
  {
    RBuffer2D result(mCols, mRows);

    #pragma omp parallel for
    for (int64_t y = 0; y < mRows; ++y)
      for (int64_t x = 0; x < mCols; ++x)
        result.pixel(x, y) = element(x, y, x, y).real();

    return result;
  }

  void apply(const CBuffer2D& buffer);

  CBuffer2D slice(int64_t y = -1) const;

  size_t numElements() const
  {
    return mCols * mCols * mRows * mRows;
  }

  size_t numBytes() const
  {
    return numElements() * sizeof(Complex);
  }

  int64_t cols() const
  {
    return mCols;
  }

  int64_t rows() const
  {
    return mRows;
  }

private:
  int64_t  mCols, mRows;

  typedef std::unique_ptr<Complex[], void(*)(void*)> Data;
  Data     mData;

  fftw_plan mFftPlan1;
  fftw_plan mFftPlan2;
  fftw_plan mFftInvPlan1;
  fftw_plan mFftInvPlan2;

  void     doFft(Direction direction);
};

} // namespace Aurora

#endif
