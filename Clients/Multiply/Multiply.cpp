//--- Aurora/Clients/Multiply/Multiply.cpp -------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/MatrixSupport.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <iostream>

namespace Aurora
{

enum class Mode
{
  abs, absSqr, current, real, scatteredCurrent
};

CL::Option<std::string> inName1("i1", "", homeDirectory());
CL::Option<std::string> inName2("i2", "", homeDirectory());
CL::Option<size_t> power1("power1", "", 1);
CL::Option<size_t> power2("power2", "", 1);
CL::Option<size_t> slices1("slices1", "", 10);
CL::Option<std::string> mode("mode", "{S, T}", "T");
CL::Option<std::string> outName("o", "", homeDirectory());

class Multiply
{
public:
  Multiply(int argc, char** argv);
};

Multiply::Multiply(int argc, char** argv)
{
  std::cerr << "Reduce Client\n";

  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
    std::cout << buildString() << '\n';

  std::string u1Name = inName1 + "/A-" + toString(slices1, 4) + ".mat";
  std::string u2Name = inName2 + "/A-0000.mat";
  std::string m1Name = inName1 + "/" + mode + "-power" + toString(power1, 4);
  std::string m2Name = inName2 + "/" + mode + "-power" + toString(power2, 4);
  std::string finalName =
    outName + "/" + mode + "-power" + toString(power1 + power2, 4);

  if (mode == "F")
  {
    MatrixNC result;
    loadMatrix(m2Name + ".mat", result);

    MatrixNC m;
    loadEigenvectors(u2Name, m);
    result *= m.adjoint();

    loadEigenvectors(u1Name, m);
    result *= m;

    loadMatrix(m1Name + ".mat", m);
    result *= m;

    saveMatrix(result, finalName + ".mat");

    return;
  }

  if (mode == "S")
  {
    swap(u1Name, u2Name);
    swap(m1Name, m2Name);
  }

  {
    MatrixNC u2;
    loadEigenvectors(u2Name, u2);

    {
      MatrixNC m2;
      loadMatrix(m2Name + "-11.mat", m2);
      MatrixNC a;
      a.noalias() = m2 * u2.adjoint();
      saveMatrix(a, tempDirectory() + "/A-11.mat");
    }

    {
      MatrixNC m2;
      loadMatrix(m2Name + "-12.mat", m2);
      MatrixNC a;
      a.noalias() = m2 * u2.adjoint();
      saveMatrix(a, tempDirectory() + "/A-12.mat");
    }

    {
      MatrixNC m2;
      loadMatrix(m2Name + "-21.mat", m2);
      MatrixNC a;
      a.noalias() = m2 * u2.adjoint();
      saveMatrix(a, tempDirectory() + "/A-21.mat");
    }

    {
      MatrixNC m2;
      loadMatrix(m2Name + "-22.mat", m2);
      MatrixNC a;
      a.noalias() = m2 * u2.adjoint();
      saveMatrix(a, tempDirectory() + "/A-22.mat");
    }
  }

  {
    MatrixNC u1;
    loadEigenvectors(u1Name, u1);

    {
      MatrixNC m1;
      loadMatrix(m1Name + "-11.mat", m1);
      MatrixNC b;
      b.noalias() = u1 * m1;
      saveMatrix(b, tempDirectory() + "/B-11.mat");
    }

    {
      MatrixNC m1;
      loadMatrix(m1Name + "-12.mat", m1);
      MatrixNC b;
      b.noalias() = u1 * m1;
      saveMatrix(b, tempDirectory() + "/B-12.mat");
    }

    {
      MatrixNC m1;
      loadMatrix(m1Name + "-21.mat", m1);
      MatrixNC b;
      b.noalias() = u1 * m1;
      saveMatrix(b, tempDirectory() + "/B-21.mat");
    }

    {
      MatrixNC m1;
      loadMatrix(m1Name + "-22.mat", m1);
      MatrixNC b;
      b.noalias() = u1 * m1;
      saveMatrix(b, tempDirectory() + "/B-22.mat");
    }
  }

  multiplyLarge(tempDirectory() + "/A", tempDirectory() + "/B", finalName);

  removeMatricesLarge(tempDirectory() + "/A");
  removeMatricesLarge(tempDirectory() + "/B");
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::Multiply multiply(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
