//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ScaledSumIC.h"

registerMooseObject("pantherApp", ScaledSumIC);

InputParameters
ScaledSumIC::validParams()
{
  InputParameters params = InitialCondition::validParams();

  params.addClassDescription("Sets the initial condition as the sum of other variables");

  params.addRequiredCoupledVar("values", "Vector of values to sum");

  params.addParam<std::vector<Real>>("prefactor", {}, "Prefactor to multiply the sum term with.");

  return params;
}

ScaledSumIC::ScaledSumIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _n_values(coupledComponents("values")),
    _prefactor(_n_values, 1.0)
{
  // we need at least one variable in the sum
  if (_n_values == 0)
    mooseError("Please supply at least one variable to sum in ScaledSumIC ", name());

  // get prefactor values if not 1.0
  std::vector<Real> p = this->template getParam<std::vector<Real>>("prefactor");

  // if prefactor is used we need the same number of prefactors as sum materials
  if (_n_values == p.size())
    _prefactor = p;
  else if (p.size() != 0)
    mooseError("Supply the same number of sum materials and prefactors.");

  for (unsigned int i = 0; i < _n_values; i++)
    _values.push_back(&coupledValue("values", i));
}

Real
ScaledSumIC::value(const Point & /*p*/)
{
  Real sum = 0;
  for (unsigned int i = 0; i < _n_values; i++)
    sum += (*(_values[i]))[_qp] * _prefactor[i];

  return sum;
}
