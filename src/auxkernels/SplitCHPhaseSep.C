//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
// 
// Author: Raghunandan Pratoori

# include "SplitCHPhaseSep.h"

registerMooseObject("pantherApp", SplitCHPhaseSep);

InputParameters
SplitCHPhaseSep::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Split formulation Cahn-Hilliard Kernel for the chemical potential variable and nonlocal term");
  params.addParam<MaterialPropertyName>("mob_name", "mobtemp", "The mobility used with the kernel");
  params.addCoupledVar("args", "Vector of variable arguments of the mobility");
  params.deprecateCoupledVar("args", "coupled_variables", "02/27/2024");
  params.addCoupledVar(
      "w", "Coupled chemical potential (if not specified kernel variable will be used)");
  params.addCoupledVar("c", "Concentration variable");
  params.addCoupledVar("bar_c", "Element average concentration variable");
  params.addParam<MaterialPropertyName>("sigma_name", "The sigma used with the kernel");
  params.addParam<MaterialPropertyName>("kappa_name", "The kappa used with the kernel");
  return params;
}

SplitCHPhaseSep::SplitCHPhaseSep(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _mob_name(getParam<MaterialPropertyName>("mob_name")),
    _mob(getMaterialProperty<Real>("mob_name")),
    _is_coupled(isCoupled("w")),
    _w_var(_is_coupled ? coupled("w") : _var.number()),
    _grad_w(_is_coupled ? coupledGradient("w") : _grad_u),
    _dmobdarg(_n_args),
    _c(coupledValue("c")),
    _bar_c(coupledValue("bar_c")),
    _sigma(getMaterialProperty<Real>("sigma_name")),
    _kappa(getMaterialProperty<Real>("kappa_name"))
{
  // Iterate over all coupled variables
  for (unsigned int i = 0; i < _n_args; ++i)
    _dmobdarg[i] = &getMaterialPropertyDerivative<Real>(_mob_name, i);

  // computeAverageConcentration();
}

// void
// SplitCHPhaseSep::computeAverageConcentration()
// {
//   _bar_c = 0.0;
//   for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
//     _bar_c += _c[qp];
//   _bar_c /= _qrule->n_points();
// }

Real
SplitCHPhaseSep::computeQpResidual()
{
  Real residual = _mob[_qp] * _grad_w[_qp] * _grad_test[_i][_qp]; //First term
  residual += _mob[_qp] * _sigma[_qp] * (_c[_qp] - _bar_c[_qp]) * _test[_i][_qp];//; //Second term
  return residual;
}

Real
SplitCHPhaseSep::computeQpJacobian()
{
  return (_is_coupled && _w_var != _var.number()) ? 0.0 : computeQpWJacobian();
}

Real
SplitCHPhaseSep::computeQpWJacobian()
{
  Real jacobian = _mob[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  jacobian += _mob[_qp] * _sigma[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  return jacobian;
}

Real
SplitCHPhaseSep::computeQpOffDiagJacobian(unsigned int jvar)
{
  // c Off-Diagonal Jacobian
  if (_w_var == jvar)
    return computeQpWJacobian();

  // get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);
  Real off_diag_jacobian = (*_dmobdarg[cvar])[_qp] * _phi[_j][_qp] * _grad_w[_qp] * _grad_test[_i][_qp];
  off_diag_jacobian += (*_dmobdarg[cvar])[_qp] * _sigma[_qp] * _phi[_j][_qp] * _test[_i][_qp];
  return off_diag_jacobian;
}
