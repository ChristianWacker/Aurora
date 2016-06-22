//--- Aurora/Clients/Matrix/PassBandlimit.hpp ----------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_CLIENTS_MATRIX_PASS_BANDLIMIT_HPP
#define AURORA_CLIENTS_MATRIX_PASS_BANDLIMIT_HPP

#include "AuroraLib/PhaseShift.hpp"
#include "Pass.hpp"

namespace Aurora
{

class PassBandlimit : public Pass
{
public:
  using Pass::Pass;
  void run() override;
};

void calcBandlimitMatrix(const MultisliceParameters& params,
                         MatrixNC& matrix);

} // namespace Aurora

#endif
