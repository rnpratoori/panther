#include "Spline2Function.h"
#include "MooseMesh.h"
#include "InputParameters.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

registerMooseObject("YourApp", Spline2Function);

InputParameters
Spline2Function::validParams()
{
  InputParameters params = validParams<Function>();
  params.addRequiredParam<FileName>("input_file", "Text file with triangle vertices and coefficients");
  return params;
}

Spline2Function::Spline2Function(const InputParameters & parameters)
  : Function(parameters)
{
  loadTrianglesFromFile(getParam<FileName>("input_file"));
}

void
Spline2Function::loadTrianglesFromFile(const std::string & filename)
{
  std::ifstream in(filename);
  if (!in)
    mooseError("Cannot open triangle data file: " + filename);

  Real x0, y0, x1, y1, x2, y2, a, b, c;
  while (in >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> a >> b >> c)
  {
    Triangle tri;
    tri.v0 = Point(x0, y0);
    tri.v1 = Point(x1, y1);
    tri.v2 = Point(x2, y2);
    tri.a = a;
    tri.b = b;
    tri.c = c;
    _triangles.push_back(tri);
  }
}

Real
Spline2Function::value(const Point & p) const
{
  for (const auto & tri : _triangles)
  {
    if (inTriangle(p, tri))
      return tri.a * p(0) + tri.b * p(1) + tri.c;
  }
  mooseError("Point not inside any triangle: (" + std::to_string(p(0)) + ", " + std::to_string(p(1)) + ")");
  return 0.0;
}

bool
Spline2Function::inTriangle(const Point & p, const Triangle & tri) const
{
  auto sign = [](const Point & p1, const Point & p2, const Point & p3)
  {
    return (p1(0) - p3(0)) * (p2(1) - p3(1)) - (p2(0) - p3(0)) * (p1(1) - p3(1));
  };

  bool b1 = sign(p, tri.v0, tri.v1) < 0.0;
  bool b2 = sign(p, tri.v1, tri.v2) < 0.0;
  bool b3 = sign(p, tri.v2, tri.v0) < 0.0;

  return ((b1 == b2) && (b2 == b3));
}