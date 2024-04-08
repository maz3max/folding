#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Aff_transformation_3.h>
#include <algorithm>
#include <iostream>

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

Polyhedron unfoldPolyhedron(const Polyhedron &P,
                            const std::vector<std::pair<Polyhedron::Facet_handle, Polyhedron::Facet_handle>> &spanningTree)
{
  Polyhedron unfoldedPolyhedron = P;

  // this will need to be replaced by a random normal
  Vector_3 downNormal(0, 0, -1);

  // Find the facet that is facing down most
  Polyhedron::Facet_handle downFacet = unfoldedPolyhedron.facets_begin();
  FT downFacetZ = 0;
  for (Polyhedron::Facet_handle facet = unfoldedPolyhedron.facets_begin(); facet != unfoldedPolyhedron.facets_end(); ++facet)
  {
    // Calculate and store the normal of the facet
    Polyhedron::Halfedge_const_handle h = facet->halfedge();
    Vector_3 normal = CGAL::cross_product(
          h->next()->vertex()->point() - h->vertex()->point(),
          h->next()->next()->vertex()->point() - h->next()->vertex()->point());
    normal /= CGAL::approximate_sqrt(normal.squared_length());
    facet->plane() = Plane_3(h->vertex()->point(), normal);

    FT z = CGAL::scalar_product(normal, downNormal);
    std::cout << "Normal: " << normal << "z: " << z << std::endl;
    if (z > downFacetZ)
    {
      downFacet = facet;
      downFacetZ = CGAL::scalar_product(normal, downNormal);
    }
    {
      downFacet = facet;
      downNormal = normal;
    }
  }
  std::cout << "Down facet: " << downFacet->plane() << std::endl;
  std::cout << "Down normal: " << downNormal << std::endl;

  // The steepest edge cut tree contains the steepest edges of all vertices, except the vertex with maximal z-coordinate.

  // Find the vertex with maximal z-coordinate
  Polyhedron::Vertex_const_handle maxZVertex = unfoldedPolyhedron.vertices_begin();
  FT maxZ = 0;
  for (Polyhedron::Vertex_const_handle vertex = unfoldedPolyhedron.vertices_begin(); vertex != unfoldedPolyhedron.vertices_end(); ++vertex)
  {
    Point_3 vertexPoint = vertex->point();
    FT z = CGAL::scalar_product(Vector_3(vertexPoint.x(), vertexPoint.y(), vertexPoint.z()), downNormal);
    if (z > maxZ)
    {
      maxZVertex = vertex;
      maxZ = z;
    }
  }

  // Find the steepest edge for each vertex
  std::vector<Polyhedron::Halfedge_const_handle> steepestEdges;
  for (Polyhedron::Vertex_const_handle vertex = unfoldedPolyhedron.vertices_begin(); vertex != unfoldedPolyhedron.vertices_end(); ++vertex)
  {
    if (vertex == maxZVertex)
    {
      continue;
    }

    Polyhedron::Halfedge_const_handle steepestEdge = vertex->halfedge();
    FT steepestEdgeZ = 0;

    // TODO: does that actually cover all outgoing edges?
    for (Polyhedron::Halfedge_around_vertex_const_circulator edge = vertex->vertex_begin(); edge != vertex->vertex_begin(); edge = ++edge)
    {
      Vector_3 edgeVector = edge->vertex()->point() - edge->opposite()->vertex()->point();
      edgeVector /= CGAL::approximate_sqrt(edgeVector.squared_length());

      FT z = CGAL::scalar_product(edgeVector, downNormal);
      if (z > steepestEdgeZ)
      {
        steepestEdge = edge;
        steepestEdgeZ = z;
      }
    }
    steepestEdges.push_back(steepestEdge);
  }

  // TODO: create MST from steepest edges


  return P;
}
