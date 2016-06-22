//--- Aurora/Clients/FormFactor/FormFactor.cpp ---------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CommandLine.hpp"
#include "AuroraLib/Element.hpp"
#include "AuroraLib/Energy.hpp"
#include "AuroraLib/Formatter.hpp"
#include "AuroraLib/Quadrature.hpp"
#include "AuroraLib/Utils.hpp"
#include "AuroraLib/Version.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace Aurora
{

enum class Mode
{
  formFactor, crossSection, potential, hybridPotential, projPotential,
  phaseShift, hybridPhaseShift
};

CL::Option<Mode> mode("mode", "", Mode::formFactor,
  enumVal(Mode::formFactor), enumVal(Mode::crossSection),
  enumVal(Mode::potential), enumVal(Mode::hybridPotential),
  enumVal(Mode::projPotential), enumVal(Mode::phaseShift),
  enumVal(Mode::hybridPhaseShift));
CL::Option<std::string> elementName("element", "Chemical Element", "C");
CL::Option<Energy> energy("energy", "Kinetic energy in eV", R(20e3));
CL::Option<Real> rho("rho", "in Angstrom", 1);
CL::Option<Real> kappa("kappa", "in 1 / Angstrom", 1);
CL::Option<Real> maxValue("max", "Maximal value", 1);
CL::Option<size_t> numSamples("numSamples", "Number of samples", 512);
CL::Option<std::string> formFactor("formFactor", "", "Peng");
CL::Option<Real> bFactor("bFactor", "Angstrom", 0);

class FormFactor
{
public:
  FormFactor(int argc, char** argv);
};

FormFactor::FormFactor(int argc, char** argv)
{
  if (verbosity >= 1)
    std::cerr << "FormFactor Version " << version << '\n';

  // add the form factor names to the option description
  formFactor.setDescription(FormFactorParametrization::namesString());

  if (!CL::parse(argc, argv))
    return;

  auto parametrization = FormFactorParametrization::create(formFactor);

  const Element& element = parametrization->element(elementName);

  std::cout << std::left;

  if (Mode::formFactor == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal k = " << maxValue << " / Angstrom\n";

    std::cout << "# Electronic form factor\n"
                 "# bFactor = " << bFactor << "\n"
                 "# k                  form factor\n";

    const Real deltaK = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real k = i * deltaK;
      std::cout << std::setw(20) << k << ' '
                << std::setw(20) << element.electronicFormFactor(k * k, bFactor)
                << '\n';
    }
  }
  else if (Mode::crossSection == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal angle = " << maxValue << '\n';

    const Real wavenumber = energy.wavenumber();

    std::cout << "# Differential elastic cross section\n"
                 "# Kinetic energy = " << energy << " eV\n"
                 "# Wavelength = " << energy.wavelength() * 100 << " pm\n"
                 "# bFactor = " << bFactor << "\n"
                 "# theta              cross section (sigma)\n";

    const Real deltaTheta = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real theta = i * deltaTheta;
      const Real kSqr = 4 * powerOf<2>(wavenumber * std::sin(R(0.5) * theta));
      const Real formfactor = element.electronicFormFactor(kSqr, bFactor);
      std::cout << std::setw(20) << theta << ' '
                << std::setw(20) << powerOf<2>(formfactor * R(1e-10)) << '\n';
    }
  }
  else if (Mode::potential == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal radius = " << maxValue << " Angstrom\n";

    std::cout << "# Potential\n"
                 "# bFactor = " << bFactor << "\n"
                 "# radius (r)         potential\n";

    const Real deltaR = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real r = i * deltaR;
      std::cout << std::setw(20) << r << ' '
                << std::setw(20) << element.potential(r * r, bFactor) << '\n';
    }
  }
  else if (Mode::hybridPotential == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal z = " << maxValue << " Angstrom\n";

    std::cout << "# Hybrid potential\n"
                 "# Kappa = " << kappa << " / Angstrom\n"
                 "# bFactor = " << bFactor << "\n"
                 "# z                  hybrid potential\n";

    const Real kappaSqr = kappa * kappa;
    const Real deltaZ = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real z = i * deltaZ;
      std::cout << std::setw(20) << z << ' '
                << std::setw(20) << element.hybridPotential(kappaSqr, z, bFactor)
                << '\n';
    }
  }
  else if (Mode::phaseShift == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal z = " << maxValue << " Angstrom\n";

    std::cout << "# Phase shift\n"
                 "# Rho = " << rho << " Angstrom\n"
                 "# bFactor = " << bFactor << "\n"
                 "# z                  phase shift\n";

    const Real deltaZ = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real z = i * deltaZ;
      std::cout << std::setw(20) << z << ' '
                << std::setw(20) << element.phaseShift(rho, z, bFactor) << '\n';
    }
  }
  else if (Mode::hybridPhaseShift == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal z = " << maxValue << " Angstrom\n";

    std::cout << "# Hybrid phase shift\n"
                 "# Kappa = " << kappa << " / Angstrom\n"
                 "# bFactor = " << bFactor << "\n"
                 "# z                  hybrid phase shift\n";

    const Real kappaSqr = kappa * kappa;
    const Real deltaZ = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real z = i * deltaZ;
      std::cout << std::setw(20) << z << ' '
                << std::setw(20) << element.hybridPhaseShift(kappaSqr, z, bFactor)
                << '\n';
    }
  }
  else if (Mode::projPotential == mode)
  {
    if (verbosity >= 1)
      std::cerr << "Maximal rho = " << maxValue << " Angstrom\n";

    std::cout << "# Projected potential\n"
                 "# bFactor = " << bFactor << "\n"
                 "# rho                projected potential\n";

    const Real deltaRho = maxValue / numSamples;
    for (size_t i = 0; i < numSamples; ++i)
    {
      const Real rho = i * deltaRho;
      std::cout << std::setw(20) << rho << ' '
                << std::setw(20) << element.projectedPotential(rho, bFactor)
                << '\n';
    }
  }
  else
  {
    AURORA_THROW(EInvalidParameter, "Unknown mode");
  }
}

} // namespace Aurora

int main(int argc, char** argv)
{
  try
  {
    Aurora::FormFactor formFactor(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
