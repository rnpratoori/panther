#include "DerivativeSpline2Material.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "libmesh/utility.h"

registerMooseObject("pantherApp", DerivativeSpline2Material);

InputParameters
DerivativeSpline2Material::validParams()
{
    InputParameters params = Material::validParams();
    params.addRequiredParam<FileName>("triangle_file", "Path to file containing triangle vertex coordinates");
    params.addRequiredParam<FileName>("coefficient_file", "Path to file containing triangle coefficients");
    params.addRequiredParam<std::string>("property_name", "Name of the free energy property");
    params.addRequiredCoupledVar("coupled_variables", "Variables to take derivatives with respect to");
    params.addParam<unsigned int>("derivative_order", 2, "Maximum order of derivatives to calculate");
    params.addClassDescription("Material that evaluates a piecewise linear spline over triangular regions with automatic derivatives");
    return params;
}

DerivativeSpline2Material::DerivativeSpline2Material(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _F(declareProperty<Real>(getParam<std::string>("property_name"))),
    _derivative_order(getParam<unsigned int>("derivative_order")),
    _tri_file(getParam<FileName>("triangle_file")),
    _coeff_file(getParam<FileName>("coefficient_file"))
{
  // Get coupled variables
  unsigned int n_vars = coupledComponents("coupled_variables");
  _vars.resize(n_vars);
  _dF.resize(n_vars);
  
  // Declare derivative properties
  for (unsigned int i = 0; i < n_vars; ++i)
  {
    _vars[i] = &coupledValue("coupled_variables", i);
    _dF[i] = &declarePropertyDerivative<Real>(
        getParam<std::string>("property_name"),
        getVar("coupled_variables", i)->name());
  }
  
  // Second derivatives if requested
  if (_derivative_order > 1)
  {
    _d2F.resize(n_vars * n_vars);
    for (unsigned int i = 0; i < n_vars; ++i)
      for (unsigned int j = 0; j < n_vars; ++j)
        _d2F[i * n_vars + j] = &declarePropertyDerivative<Real>(
            getParam<std::string>("property_name"),
            getVar("coupled_variables", i)->name(),
            getVar("coupled_variables", j)->name());
  }

  // Load triangle vertex data
  std::ifstream tfile(_tri_file);
  if (!tfile)
    mooseError("Failed to open triangle file: " + _tri_file);

  Real x1, y1, x2, y2, x3, y3;
  while (tfile >> x1 >> y1 >> x2 >> y2 >> x3 >> y3)
    _triangles.push_back({x1, y1, x2, y2, x3, y3});
  tfile.close();

  // Load triangle coefficients
  std::ifstream cfile(_coeff_file);
  if (!cfile)
    mooseError("Failed to open coefficient file: " + _coeff_file);

  Real a, b, c;
  while (cfile >> a >> b >> c)
    _coefficients.push_back({a, b, c});
  cfile.close();

  if (_triangles.size() != _coefficients.size())
    mooseError("Mismatch: number of triangles and coefficient sets must be equal");
}

void
DerivativeSpline2Material::computeQpProperties()
{
    const Real x = _q_point[_qp](0);
    const Real y = _q_point[_qp](1);
    
    _F[_qp] = evaluateFunction(x, y);
    
    // Calculate derivatives here if needed
}

Real
DerivativeSpline2Material::evaluateFunction(Real x, Real y) const
{
for (std::size_t i = 0; i < _triangles.size(); ++i)
{
if (pointInTriangle(x, y, _triangles[i]))
{
const auto & coeffs = _coefficients[i];
return coeffs[0] * x + coeffs[1] * y + coeffs[2];
}
}
return 0.0; // Or throw an error if x, y is outside all triangles
}

bool
DerivativeSpline2Material::pointInTriangle(Real x, Real y, const std::array<Real, 6> & tri) const
{
Real x1 = tri[0], y1 = tri[1];
Real x2 = tri[2], y2 = tri[3];
Real x3 = tri[4], y3 = tri[5];

Real det = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
Real a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / det;
Real b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / det;
Real c = 1 - a - b;

return (a >= 0) && (b >= 0) && (c >= 0);
}