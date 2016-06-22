//--- Aurora/AuroraLib/ElementDoyle.cpp ----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/ElementGaussian.hpp"
#include "AuroraLib/ElementHydrogen.hpp"
#include "AuroraLib/Exception.hpp"

namespace Aurora
{

class AURORA_API FormFactorParametrizationDoyle :
  public FormFactorParametrization
{
public:
  FormFactorParametrizationDoyle()
  {
    addElement(std::make_unique<ElementHydrogen>());
    addElement(std::make_unique<ElementHydrogenIon>());

    // [Doyle1968], table 4a
    addElement(ConstElementPtr(new ElementGaussian<4>( 2, 0, // HE
      {{R(  0.0906), R( 0.1814), R(0.1095), R(0.0362)}},
      {{R( 18.1834), R( 6.2109), R(1.8026), R(0.2844)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>( 6, 0, // C
      {{R(  0.7307), R( 1.1951), R(0.4563), R(0.1247)}},
      {{R( 36.9951), R(11.2966), R(2.8139), R(0.3456)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>( 7, 0, // N
      {{R(  0.5717), R( 1.0425), R(0.4647), R(0.1311)}},
      {{R( 28.8465), R( 9.0542), R(2.4213), R(0.3167)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>( 8, 0, // O
      {{R(  0.4548), R( 0.9173), R(0.4719), R(0.1384)}},
      {{R( 23.7803), R( 7.6220), R(2.1440), R(0.2959)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>(11, 0, // NA
      {{R(  2.2406), R( 1.3326), R(0.9070), R(0.2863)}},
      {{R(108.0039), R(24.5047), R(3.3914), R(0.4346)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>(12, 0, // MG
      {{R(  2.2682), R( 1.8025), R(0.8394), R(0.2892)}},
      {{R( 73.6704), R(20.1749), R(3.0181), R(0.4046)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>(14, 0, // SI
      {{R(  2.1293), R( 2.5333), R(0.8349), R(0.3216)}},
      {{R( 57.7748), R(16.4756), R(2.8796), R(0.3860)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>(15, 0, // P
      {{R(  1.8882), R( 2.4685), R(0.8046), R(0.3204)}},
      {{R( 44.8756), R(13.5383), R(2.6424), R(0.3608)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>(16, 0, // S
      {{R(  1.6591), R( 2.3863), R(0.7899), R(0.3208)}},
      {{R( 36.6500), R(11.4881), R(2.4686), R(0.3403)}})));
    addElement(ConstElementPtr(new ElementGaussian<4>(92, 0, // U
      {{R(  6.7668), R( 6.7287), R(4.0135), R(1.5607)}},
      {{R( 85.9510), R(15.6415), R(2.9364), R(0.3348)}})));
  }
};

class AURORA_API RegisterDoyle
{
public:
  RegisterDoyle()
  {
    FormFactorParametrization::registerClass<FormFactorParametrizationDoyle>
      ("Doyle");
  }
} registerDoyle;

} // namespace Aurora
