#include "SplineMaterial.h"

registerMooseObject("pantherApp", SplineMaterial);

InputParameters
SplineMaterial::validParams()
{
  InputParameters params = ParsedMaterialTempl::validParams();
  params.addRequiredParam<FunctionName>("function", "The name of the spline function");
  params.addRequiredCoupledVar("c", "The variable on which the spline function depends");

  return params;
}

SplineMaterial::SplineMaterial(const InputParameters & parameters)
  : ParsedMaterialTempl(parameters),
    _func(static_cast<const SplineFunction &>(getFunction("function"))),
    _f(declareProperty<Real>("f")),
    _df_dc(declareProperty<Real>("df_dc")),
    _d2f_dc2(declareProperty<Real>("d2f_dc2")),
    _c(coupledValue("c"))
{
}

void
SplineMaterial::computeQpProperties()
{
    Real c_val = _c[_qp]; // concentration value at this quad point
    Point  p(c_val);

    _f[_qp]       = _func.value(0.0, p);
    _df_dc[_qp]   = _func.derivative(c_val);
    _d2f_dc2[_qp] = _func.secondDerivative(c_val);

}
