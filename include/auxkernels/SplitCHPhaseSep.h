// Author: Raghunandan Pratoori

#pragma once

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

/**
 * SplitCHWresBase implements the residual for the chemical
 * potential in the split form of the Cahn-Hilliard
 * equation in a general way that can be templated to a scalar or
 * tensor mobility.
 */
class SplitCHPhaseSep : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  static InputParameters validParams();

  SplitCHPhaseSep(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpWJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// calculate the average of concentration
  void computeAverageConcentration();

  const MaterialPropertyName _mob_name;
  const MaterialProperty<Real> & _mob;

  /// is the kernel used in a coupled form?
  const bool _is_coupled;

  /// int label for the chemical potential
  unsigned int _w_var;

  /// Variable value for the chemical potential
  const VariableGradient & _grad_w;

  /// derivatives of the mobility
  std::vector<const MaterialProperty<Real> *> _dmobdarg;

  /// Variables for nonlocal contribution
  const VariableValue & _c;
  const VariableValue & m;
  const MaterialProperty<Real> & _sigma;
  const MaterialProperty<Real> & _kappa;
};
