// Author: Raghunandan Pratoori 

#pragma once

#include "RandomICBase.h"
#include "DistributionInterface.h"

// Forward Declarations
class InputParameters;
class Distribution;
namespace libMesh
{
class Point;
}

template <typename T>
InputParameters validParams();

/**
 * RandomConstraintIC just returns a Random value on
 * selected elements. The elements are selected based on the
 * values of the coupled varible.
 */
class RandomConstraintIC : public RandomICBase, public DistributionInterface
{
public:
  /**
   * Constructor
   * @param parameters The parameters object holding data for the class to use.
   */
  RandomConstraintIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

  static InputParameters validParams();

protected:
  /// The lower bound of the random number range
  const Real _min;

  /// The upper bound of the random number range
  const Real _max;

  /// Distribution object optionally used to define distribution of random numbers
  Distribution const * _distribution;

  /// Coupled varible
  const VariableValue & _coupled_val;
};
