//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AirFluidProperties.h"

registerMooseObject("acaesApp", AirFluidProperties);

InputParameters
AirFluidProperties::validParams()
{
  InputParameters params = IdealGasFluidProperties::validParams();
  params.addClassDescription("Fluid properties for air");
  return params;
}

AirFluidProperties::AirFluidProperties(const InputParameters & parameters)
  : IdealGasFluidProperties(parameters),
    _p_critical(3.786e6),
    _T_critical(132.5306),
    _rho_critical(342.685),
    _p_triple(5.264181e3),
    _T_triple(59.75)
{
}

std::string
AirFluidProperties::fluidName() const
{
  return "air";
}

Real
AirFluidProperties::molarMass() const
{
  return _molar_mass;
}

Real
AirFluidProperties::criticalTemperature() const
{
  return _T_critical;
}

Real
AirFluidProperties::criticalPressure() const
{
  return _p_critical;
}

Real
AirFluidProperties::criticalDensity() const
{
  return _rho_critical;
}

std::vector<Real>
AirFluidProperties::henryCoefficients() const
{
  return {-9.625741, 4.659272, 11.642974};
}
