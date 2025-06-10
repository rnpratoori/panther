#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"

class DerivativeSpline2Material : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  DerivativeSpline2Material(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  Real evaluateFunction(Real x, Real y) const;
  bool pointInTriangle(Real x, Real y, const std::array<Real, 6> & tri) const;

  /// The free energy property
  MaterialProperty<Real> & _F;
  
  /// First derivatives wrt coupled variables
  std::vector<MaterialProperty<Real> *> _dF;
  
  /// Second derivatives wrt coupled variables
  std::vector<MaterialProperty<Real> *> _d2F;
  
  /// Coupled variable values
  std::vector<const VariableValue *> _vars;
  
  /// Maximum order of derivatives to calculate
  const unsigned int _derivative_order;

  /// Triangle vertices data
  std::vector<std::array<Real, 6>> _triangles;
  
  /// Triangle coefficients data
  std::vector<std::array<Real, 3>> _coefficients;

  /// Input files
  const FileName & _tri_file;
  const FileName & _coeff_file;
};