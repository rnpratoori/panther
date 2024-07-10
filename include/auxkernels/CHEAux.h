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
