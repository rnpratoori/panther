#pragma once

#include "SplitCHCRes.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

/**
 * SplitCHSpline: Cahn-Hilliard split kernel using spline-based free energy (supports multicomponent).
 */
class SplitCHSpline : public DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHCRes>>
{
public:
  static InputParameters validParams();
  SplitCHSpline(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real computeDFDC(PFFunctionType type) override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const MaterialProperty<Real> & _dfdc;
  const MaterialProperty<Real> & _d2fdc2;
  std::vector<const MaterialProperty<Real> *> _d2fdcdarg;
};
