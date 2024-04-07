#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Aff_transformation_3.h>
#include <Eigen/Geometry>
#include <cmath>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// Function to create a rotation transformation using Eigen quaternions
Aff_transformation_3 createRotationTransformation(const Vector_3 &axis, double angle)
{
  Eigen::Vector3d axis_eigen(axis.x(), axis.y(), axis.z());
  Eigen::Quaterniond quaternion(Eigen::AngleAxisd(angle, axis_eigen.normalized()));
  Eigen::Matrix3d rotationMatrix = quaternion.toRotationMatrix();

  // Convert Eigen matrix to CGAL Aff_transformation_3
  return Aff_transformation_3(rotationMatrix(0, 0), rotationMatrix(0, 1), rotationMatrix(0, 2), 0,
                              rotationMatrix(1, 0), rotationMatrix(1, 1), rotationMatrix(1, 2), 0,
                              rotationMatrix(2, 0), rotationMatrix(2, 1), rotationMatrix(2, 2), 0);
}

Aff_transformation_3 getFallingDownTransformation(const Point_3 &facetPoint, const Vector_3 &facetNormal, const Plane_3 &plane)
{
  // Project the facetPoint onto the plane
  Point_3 projectedPoint = plane.projection(facetPoint);

  // Calculate translation vector
  Vector_3 translationVector = projectedPoint - facetPoint;
  Aff_transformation_3 translate(CGAL::TRANSLATION, translationVector);

  // Check if the normals are identical
  Vector_3 planeNormal = plane.orthogonal_vector();
  if (planeNormal == facetNormal)
  {
    return translate;
  }

  // Calculate rotation
  Vector_3 rotationAxis = CGAL::cross_product(facetNormal, planeNormal);
  rotationAxis = rotationAxis / CGAL::approximate_sqrt(rotationAxis.squared_length());
  double angle = std::acos(CGAL::scalar_product(facetNormal, planeNormal));

  Aff_transformation_3 rotate = createRotationTransformation(rotationAxis, angle);

  // Combine translation and rotation
  return translate * rotate;
}

/*
Idea:

A spanning tree (dual to a cut tree) is generated using steepest edge cut first.
After that, the polyhedron is unfolded along the spanning tree:

1. Find the facet that is facing down most (normal vector is closest to (0, 0, -1))
2. Use getFallingDownTransformation to get the transformation that makes the facet fall down
3. Along the spanning tree, each facet inherits the transformation from its parent
4. Additionally, each facet is rotated around the common edge with its parent

Lastly, intersections between facets are checked. The process is repeated until no intersections are found.
*/

int main(int argc, char *argv[])
{
  Polyhedron P;
  std::ifstream in1((argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off"));
  in1 >> P;
  CGAL::draw(P);

  return EXIT_SUCCESS;
}
