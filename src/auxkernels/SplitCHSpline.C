#include "SplitCHSpline.h"

registerMooseObject("pantherApp", SplitCHSpline);

template <>
InputParameters
validParams<SplitCHSpline>()
{
  InputParameters params = validParams<SplitCHCRes>();
  params.addRequiredParam<MaterialPropertyName>("f_name", "Base name of spline free energy material");
  return params;
}

SplitCHSpline::SplitCHSpline(const InputParameters & parameters)
  : SplitCHCRes(parameters),
    _dfdc(getMaterialProperty<Real>(parameters.get<MaterialPropertyName>("f_name") + "_df_dc")),
    _d2fdc2(getMaterialProperty<Real>(parameters.get<MaterialPropertyName>("f_name") + "_d2f_dc2"))
{
}

Real
SplitCHSpline::computeDFDC(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
      return _dfdc[_qp];
    case Jacobian:
      return _d2fdc2[_qp] * _phi[_j][_qp];
  }

  mooseError("Unknown PFFunctionType");
}

Real
SplitCHSpline::computeQpOffDiagJacobian(unsigned int jvar)
{
  // If this is the w-variable, let the base class handle it
  if (jvar == _w_var)
    return SplitCHCRes::computeQpOffDiagJacobian(jvar);

  // For now, assume no coupling with other variables (no args)
  return 0.0;
}
