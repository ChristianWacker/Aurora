//--- Aurora/Clients/CrossSection/CrossSection.cpp -----------------------------
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
#include <iostream>

namespace Aurora
{

enum class Mode
{
  elastic, inelastic, ionicEnergies, r2
};

CL::Option<Mode> mode("mode", "", Mode::elastic,
  enumVal(Mode::elastic), enumVal(Mode::inelastic),
  enumVal(Mode::ionicEnergies), enumVal(Mode::r2));
CL::Option<Energy> energy("energy", "Kinetic energy in eV", R(20e3));
CL::Option<std::string> formFactor("formFactor", "", "Peng");
CL::Option<Real> minIntegration("minIntegration", "", 0);
CL::Option<Real> maxIntegration("maxIntegration", "", Math::pi);

class CrossSection
{
public:
  CrossSection(int argc, char** argv);
};

// todo: literature.
// https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
const std::array<Real, 104> ionicEnergies =
  {{R(13.59844), R(24.58741), R( 5.39172), R( 9.3227 ), R( 8.29803),
    R(11.26030), R(14.53414), R(13.61806), R(17.42282), R(21.5646 ),
    R( 5.13908), R( 7.64624), R( 5.98577), R( 8.15169), R(10.48669),
    R(10.36001), R(12.96764), R(15.75962), R( 4.34066), R( 6.11316),
    R( 6.5615 ), R( 6.8281 ), R( 6.7462 ), R( 6.7665 ), R( 7.43402),
    R( 7.9024 ), R( 7.8810 ), R( 7.6398 ), R( 7.72638), R( 9.3942 ),
    R( 5.99930), R( 7.8994 ), R( 9.7886 ), R( 9.75238), R(11.81381),
    R(13.99961), R( 4.17713), R( 5.6949 ), R( 6.2171 ), R( 6.63390),
    R( 6.75885), R( 7.09243), R( 7.28   ), R( 7.36050), R( 7.45890),
    R( 8.3369 ), R( 7.5762 ), R( 8.9938 ), R( 5.78636), R( 7.3439 ),
    R( 8.6084 ), R( 9.0096 ), R(10.45126), R(12.1298 ), R( 3.89390),
    R( 5.21170), R( 5.5769),  R( 5.5387 ), R( 5.473  ), R( 5.5250 ),
    R( 5.582  ), R( 5.6436 ), R( 5.6704 ), R( 6.1501 ), R( 5.8638 ),
    R( 5.9389 ), R( 6.0215 ), R( 6.1077 ), R( 6.18431), R( 6.25416),
    R( 5.4259 ), R( 6.82507), R( 7.5496 ), R( 7.8640 ), R( 7.8335 ),
    R( 8.4382 ), R( 8.9670 ), R( 8.9587 ), R( 9.2255 ), R(10.43750),
    R( 6.1082 ), R( 7.41666), R( 7.2856 ), R( 8.417  ), R( 9.31751),
    R(10.74850), R( 4.0727 ), R( 5.2784 ), R( 5.17   ), R( 6.3067 ),
    R( 5.89   ), R( 6.19405), R( 6.2657 ), R( 6.0262 ), R( 5.9738 ),
    R( 5.9915 ), R( 6.1979 ), R( 6.2817 ), R( 6.42   ), R( 6.50   ),
    R( 6.58   ), R( 6.65   ), R( 4.9    ), R( 6.0    )}};

CrossSection::CrossSection(int argc, char** argv)
{
  if (verbosity >= 1)
    std::cerr << "CrossSection Client\n";

  // add the form factor names to the option description
  formFactor.setDescription(FormFactorParametrization::namesString());

  if (!CL::parse(argc, argv))
    return;

  if (verbosity >= 1)
    std::cerr << buildString() << '\n';

  auto parametrization = FormFactorParametrization::create(formFactor);

  if (Mode::elastic == mode)
  {
    const Real wavenumber = energy.wavenumber();
    const Real gamma = energy.lorentzFactor();

    std::cout << "# Total elastic cross section\n"
                 "# Kinetic energy = " << energy << " eV\n"
                 "# Wavelength = " << energy.wavelength() * 100 << " pm\n"
                 "# Relativistic factor (gamma): " << gamma << "\n"
                 "#\n"
                 "# atomic number (Z)            cross section (sigma)\n";

    for (size_t i = 1; i <= Element::mNumElements; ++i)
    {
      const Element& element = parametrization->element(i, 0);

      auto integrand = [&] (Real theta)
      {
        const Real kSqr = 4 * powerOf<2>(wavenumber * std::sin(R(0.5) * theta));
        return powerOf<2>(gamma * element.electronicFormFactor(kSqr))
          * std::sin(theta);
      };

      const Real totalCrossSection =
        Math::twoPi * integrateOpen(integrand, minIntegration, maxIntegration);
      std::cout << i << ' ' << totalCrossSection << '\n';
    }
  }
  else if (Mode::inelastic == mode)
  {
    const Real wavenumber = energy.wavenumber();

    std::cout << "# Total inelastic cross section\n"
                 "# Kinetic energy = " << energy << " eV\n"
                 "# Wavelength = " << energy.wavelength() * 100 << " pm\n"
                 "#\n"
                 "# atomic number (Z)            cross section (sigma)\n";

    for (size_t i = 1,
         e = std::min(Element::mNumElements, ionicEnergies.size());
         i <= e; ++i)
    {
      const Element& element = parametrization->element(i, 0);

      auto integrand = [&] (Real theta)
      {
        const Real thetaE = ionicEnergies[i - 1] / 4 / energy;
        const Real kSqr =
          powerOf<2>(wavenumber) * (powerOf<2>(theta) + powerOf<2>(thetaE));

        Real result = i - powerOf<2>(element.xRayFormFactor(kSqr)) / i;
        result /= powerOf<2>(kSqr);
        result *= powerOf<2>(2 / Phys::a0) * std::sin(theta);
        return result;
      };

      const Real totalCrossSection =
        Math::twoPi * integrateOpen(integrand, minIntegration, maxIntegration);
      std::cout << i << ' ' << totalCrossSection << '\n';
    }
  }
  else if (Mode::ionicEnergies == mode)
  {
    std::cout << "# Ionic energies\n"
                 "#\n"
                 "# atomic number (Z)      (E) [eV]\n";

    for (size_t i = 1; i <= ionicEnergies.size(); ++i)
      std::cout << i << ' ' << ionicEnergies[i - 1] << '\n';
  }
  else if (Mode::r2 == mode)
  {
    std::cout << "# Mean quadratic distance\n"
                 "#\n"
                 "# atomic number (Z)      R^2\n";

    for (size_t i = 1; i <= Element::mNumElements; ++i)
    {
      const Element& element = parametrization->element(i, 0);
      std::cout << i << ' ' << element.radiusSqr() << '\n';
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
    Aurora::CrossSection crossSection(argc, argv);
  }
  catch (std::exception& e)
  {
    std::cerr << Aurora::Color::red << Aurora::Terminal::bold
              << "An exception occured:\n"
              << Aurora::Terminal::reset << e.what() << '\n';
  }
}
