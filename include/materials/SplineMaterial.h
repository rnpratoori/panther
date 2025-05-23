#pragma once

#include "ParsedMaterial.h"
#include "SplineFunction.h"

class SplineMaterial : public ParsedMaterialTempl<false>
{
public:
    static InputParameters validParams();
    SplineMaterial(const InputParameters & parameters);

protected:
    virtual void computeQpProperties() override;

    const SplineFunction & _func;

    MaterialProperty<Real> & _f;
    MaterialProperty<Real> & _df_dc;
    MaterialProperty<Real> & _d2f_dc2;
    const VariableValue & _c;
};
