#pragma once

#include "Function.h"
#include <vector>

class Spline2Function : public Function
{
public:
  static InputParameters validParams();

  Spline2Function(const InputParameters & parameters);

  virtual Real value(const Point & p) const override;

private:
  struct Triangle {
    Point v0;
    Point v1;
    Point v2;
    Real a, b, c; // coefficients for z = a*x + b*y + c
  };

  std::vector<Triangle> _triangles;

  void loadTrianglesFromFile(const std::string & filename);
  bool inTriangle(const Point & p, const Triangle & tri) const;
};