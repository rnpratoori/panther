// Author: Raghunandan Pratoori

#pragma once

#include "AuxKernel.h"

/**
 * Coupled auxiliary value 
 * Used to calculate the polymer volume fraction equivalent to experiments
 * pvf = (c + 1)/2
 */
class CHEAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  CHEAux(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeValue() override;

  const VariableValue & _coupled_val;
};
