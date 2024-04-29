//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowEnthalpySink.h"

class SinglePhaseFluidProperties;
class PorousFlowDictator;
class Function;

/**
 * Applies a flux sink of heat energy to a boundary with specified mass flux and inlet temperature.
 */
class ACAESIsentropicCompBC : public PorousFlowEnthalpySink
{
public:
  static InputParameters validParams();

  ACAESIsentropicCompBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Derivative of residual with respect to the jvar variable
  Real jac(unsigned int jvar) const;

  /// Isentropic efficiency of compression (default is 0.85)
  const Real _iec;
  /// Ambient pressure
  const Real _p_ambient;
  /// Ambient temperature
  const Real _T_ambient;
  /// Enthalpy at ambient pressure and temperature
  const Real _h_ambient;
  /// Entropy at ambient pressure and temperature
  const Real _s_ambient;
};
