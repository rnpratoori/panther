// Author: Raghunandan Pratoori

#include "RandomConstraintIC.h"

#include "libmesh/point.h"
#include "Distribution.h"

registerMooseObject("pantherApp", RandomConstraintIC);

InputParameters
RandomConstraintIC::validParams()
{
  InputParameters params = RandomICBase::validParams();
  params += DistributionInterface::validParams();
  params.addParam<Real>(
      "min", 0.0, "Lower bound of uniformly distributed randomly generated values");
  params.addParam<Real>(
      "max", 1.0, "Upper bound of uniformly distributed randomly generated values");
  params.addParam<DistributionName>(
      "distribution", "Name of distribution defining distribution of randomly generated values");

  params.addClassDescription("Initialize a variable with randomly generated numbers following "
                             "either a uniform distribution or a user-defined distribution");
  params.addRequiredCoupledVar("coupled", "Coupled variable");
  return params;
}

RandomConstraintIC::RandomConstraintIC(const InputParameters & parameters)
  : RandomICBase(parameters),
    DistributionInterface(this),
    _min(getParam<Real>("min")),
    _max(getParam<Real>("max")),
    _distribution(nullptr),
    _coupled_val(coupledValue("coupled"))
{
  if (_min >= _max)
    paramError("min", "Min >= Max for RandomConstraintIC!");

  if (parameters.isParamSetByUser("distribution"))
  {
    _distribution = &getDistributionByName(getParam<DistributionName>("distribution"));
    if (parameters.isParamSetByUser("min") || parameters.isParamSetByUser("max"))
      paramError("distribution", "Cannot use together with 'min' or 'max' parameter");
  }
}

Real
RandomConstraintIC::value(const Point & /*p*/)
{
    // 1 - coupled_val[_qp] gives the phase fraction to be
    // occupied by the other phases
    // Phase fraction of the current variable is scaled to
    // the avaiable range
    if (_coupled_val[_qp] < 1)
    {
        if (_distribution)
            return _distribution->quantile(generateRandom() * (1 - _coupled_val[_qp]));
        else
            return generateRandom() * (1 - _coupled_val[_qp]) * (_max - _min) + (1 - _coupled_val[_qp]) * _min;
    }
    else
        return 0.0;  
}
