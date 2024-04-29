//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACAESIsentropicCompBC.h"
#include "Function.h"
#include "PorousFlowDictator.h"
#include "SinglePhaseFluidProperties.h"
// #include "MooseVariable.h"
// #include "PorousFlowSinkBC.h"

// #include "libmesh/quadrature.h"

// #include <iostream>

registerMooseObject("acaesApp", ACAESIsentropicCompBC);

InputParameters
ACAESIsentropicCompBC::validParams()
{
  InputParameters params = PorousFlowEnthalpySink::validParams();
  params.suppressParameter<Real>("T_in");
  params.set<Real>("T_in", 300);
  params.addParam<Real>("iec", 0.85, "Isentropic compression efficiency. Detault is 0.85");
  params.addClassDescription(
      "Applies a source equal to the product of the mass flux and the "
      "fluid enthalpy. The enthalpy is computed at temperature T_in and pressure equal to the "
      "porepressure in the porous medium, if fluid_phase is given, otherwise at the supplied "
      "porepressure. Hence this adds heat energy to the porous medium at rate corresponding to a "
      "fluid being injected at (porepressure, T_in) at rate (-flux_function).");
  return params;
}

ACAESIsentropicCompBC::ACAESIsentropicCompBC(const InputParameters & parameters)
  : PorousFlowEnthalpySink(parameters),
    _iec(getParam<Real>("iec")),
    _p_ambient(1.01325e5),
    _T_ambient(298.15),
    _h_ambient(_fp.h_from_p_T(_p_ambient, _T_ambient)),
    _s_ambient(_fp.s_from_p_T(_p_ambient, _T_ambient))
{
}

Real
ACAESIsentropicCompBC::computeQpResidual()
{
  Real p;
  if (_pressure)
    p = (*_pressure)[_qp];
  else
    p = _pp[_i][_ph];

  const Real gamma = _fp.gamma_from_p_T(_p_ambient, _T_ambient);
  const Real T =
      _T_ambient * (1.0 / _iec * (std::pow(p / _p_ambient, (gamma - 1.0) / gamma) - 1.0) + 1.0);

  const Real h = _fp.h_from_p_T(p, T);

  return _test[_i][_qp] * _m_func.value(_t, _q_point[_qp]) * h;
}

Real
ACAESIsentropicCompBC::computeQpJacobian()
{
  return jac(_var.number());
}

Real
ACAESIsentropicCompBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  return jac(jvar);
}

Real
ACAESIsentropicCompBC::jac(unsigned int jvar) const
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);

  if (_i != _j)
    return 0.0;

  Real jac;
  if (_pressure)
  {
    jac = 0.;
  }
  else
  {
    Real h, dh_dpp, dh_dT;
    const Real gamma = _fp.gamma_from_p_T(_p_ambient, _T_ambient);
    const Real T =
        _T_ambient *
        (1.0 / _iec * (std::pow(_pp[_i][_ph] / _p_ambient, (gamma - 1.0) / gamma) - 1.0) + 1.0);

    const Real dT_dp = _T_ambient / (_iec * _p_ambient) * ((gamma - 1.0) / gamma) *
                       std::pow(_pp[_i][_ph] / _p_ambient, -1.0 / gamma);

    _fp.h_from_p_T(_pp[_i][_ph], T, h, dh_dpp, dh_dT);

    jac = dh_dpp * _dpp_dvar[_i][_ph][pvar] + dh_dT * dT_dp * _dpp_dvar[_i][_ph][pvar];
  }
  return _test[_i][_qp] * _m_func.value(_t, _q_point[_qp]) * jac;
}
