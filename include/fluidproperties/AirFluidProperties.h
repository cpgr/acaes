//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IdealGasFluidProperties.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

/**
 * Air fluid properties (where air is represented as an ideal gas)
 * Default parameters are for air at atmospheric pressure and temperature
 */
class AirFluidProperties : public IdealGasFluidProperties
{
public:
  static InputParameters validParams();

  AirFluidProperties(const InputParameters & parameters);

  virtual std::string fluidName() const override;
  virtual Real molarMass() const override;
  virtual Real criticalTemperature() const override;
  virtual Real criticalPressure() const override;
  virtual Real criticalDensity() const override;
  /// 78:22 weighted Henry coefficients of N2 and O2
  virtual std::vector<Real> henryCoefficients() const override;

protected:
  /// Critical pressure (Pa)
  const Real _p_critical;
  /// Critical temperature (K)
  const Real _T_critical;
  /// Critical density (kg/m^3)
  const Real _rho_critical;
  /// Triple point pressure (Pa)
  const Real _p_triple;
  /// Triple point temperature (K)
  const Real _T_triple;
};

#pragma GCC diagnostic pop
