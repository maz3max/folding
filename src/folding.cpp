#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Aff_transformation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::FT FT;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

Aff_transformation_3 createRotationTransformation(const Vector_3 &normal_a, const Vector_3 &normal_b)
{
  if (normal_a == normal_b)
  {
    return Aff_transformation_3(CGAL::IDENTITY);
  }
  else if (normal_a == -normal_b)
  {
    return Aff_transformation_3(-1, 0, 0, 0, -1, 0, 0, 0, -1);
  }

  // get normalized rotation axis
  Vector_3 rotation_axis = CGAL::cross_product(normal_a, normal_b);
  rotation_axis /= CGAL::approximate_sqrt(rotation_axis.squared_length());

  // get sine and cosine of the angle (the actual angle is not needed)
  const FT cos_x = CGAL::scalar_product(normal_a, normal_b);
  const FT sin_x = CGAL::approximate_sqrt(CGAL::cross_product(normal_a, normal_b).squared_length());

  // build rotation matrix using Rodrigues' rotation formula
  const FT m00 = cos_x + (rotation_axis.x() * rotation_axis.x()) * (1 - cos_x);
  const FT m01 = (rotation_axis.x() * rotation_axis.y()) * (1 - cos_x) - rotation_axis.z() * sin_x;
  const FT m02 = (rotation_axis.x() * rotation_axis.z()) * (1 - cos_x) + rotation_axis.y() * sin_x;
  const FT m10 = (rotation_axis.y() * rotation_axis.x()) * (1 - cos_x) + rotation_axis.z() * sin_x;
  const FT m11 = cos_x + (rotation_axis.y() * rotation_axis.y()) * (1 - cos_x);
  const FT m12 = (rotation_axis.y() * rotation_axis.z()) * (1 - cos_x) - rotation_axis.x() * sin_x;
  const FT m20 = (rotation_axis.z() * rotation_axis.x()) * (1 - cos_x) - rotation_axis.y() * sin_x;
  const FT m21 = (rotation_axis.z() * rotation_axis.y()) * (1 - cos_x) + rotation_axis.x() * sin_x;
  const FT m22 = cos_x + (rotation_axis.z() * rotation_axis.z()) * (1 - cos_x);
  return Aff_transformation_3(m00, m01, m02, m10, m11, m12, m20, m21, m22);
}

Aff_transformation_3 getFallingDownTransformation(const Point_3 &facetPoint, const Vector_3 &facetNormal, const Plane_3 &plane)
{
  // Project the facetPoint onto the plane
  Point_3 projectedPoint = plane.projection(facetPoint);

  // Calculate translation vector
  Vector_3 translationVector = projectedPoint - facetPoint;
  Aff_transformation_3 translate(CGAL::TRANSLATION, translationVector);

  // Calculate rotation matrix
  Vector_3 planeNormal = plane.orthogonal_vector();
  Aff_transformation_3 rotate = createRotationTransformation(facetNormal, planeNormal);

  // Combine translation and rotation
  return translate * rotate;
}
