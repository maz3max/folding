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

typedef std::vector<Polyhedron::Halfedge_const_handle> SpanningTree;

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

  Vector_3 facetPointVector = Vector_3(facetPoint.x(), facetPoint.y(), facetPoint.z());

  // Calculate translation vector
  Aff_transformation_3 translate(CGAL::TRANSLATION, projectedPoint - facetPoint);

  // Calculate rotation matrix
  // (first, translate the facetPoint to the origin, then rotate, then translate back to the original position)
  Vector_3 planeNormal = plane.orthogonal_vector();
  Aff_transformation_3 rotate =
      Aff_transformation_3(CGAL::TRANSLATION, facetPointVector) *
      createRotationTransformation(facetNormal, planeNormal) *
      Aff_transformation_3(CGAL::TRANSLATION, -facetPointVector);

  // Combine translation and rotation
  return translate * rotate;
}

Vector_3 getFaceNormal(const Polyhedron::Facet_const_handle facet)
{
  Polyhedron::Halfedge_const_handle h = facet->halfedge();
  Vector_3 normal = CGAL::cross_product(
      h->next()->vertex()->point() - h->vertex()->point(),
      h->next()->next()->vertex()->point() - h->next()->vertex()->point());
  normal /= CGAL::approximate_sqrt(normal.squared_length());

  return normal;
}

// generate transform to unfold the face on the opposing edge along the edges to match the face along the current edge
Aff_transformation_3 getFaceUnfoldTransformation(const Polyhedron::Halfedge_const_handle edge)
{
  Vector_3 targetNormal = getFaceNormal(edge->facet());
  Vector_3 currentNormal = getFaceNormal(edge->opposite()->facet());
  Point_3 edgePoint = edge->vertex()->point();
  Vector_3 edgePointVector = Vector_3(edgePoint.x(), edgePoint.y(), edgePoint.z());

  return Aff_transformation_3(CGAL::TRANSLATION, edgePointVector) *
         createRotationTransformation(currentNormal, targetNormal) *
         Aff_transformation_3(CGAL::TRANSLATION, -edgePointVector);
}

// The steepest edge cut tree contains the steepest edges of all vertices, except the vertex with maximal z-coordinate.
SpanningTree steepestEdgeCut(const Polyhedron &P, const Vector_3 downNormal)
{
  SpanningTree result;

  // Find the vertex with maximal z-coordinate
  Polyhedron::Vertex_const_handle maxZVertex = P.vertices_begin();
  FT maxZ = 0;
  for (Polyhedron::Vertex_const_handle vertex = P.vertices_begin(); vertex != P.vertices_end(); ++vertex)
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
  for (Polyhedron::Vertex_const_handle vertex = P.vertices_begin(); vertex != P.vertices_end(); ++vertex)
  {
    if (vertex == maxZVertex)
    {
      continue;
    }

    Polyhedron::Halfedge_const_handle steepestEdge = vertex->halfedge();
    FT steepestEdgeZ = 0;

    Polyhedron::Halfedge_around_vertex_const_circulator edge = vertex->vertex_begin();
    do
    {
      Vector_3 edgeVector = edge->vertex()->point() - edge->opposite()->vertex()->point();
      edgeVector /= CGAL::approximate_sqrt(edgeVector.squared_length());

      FT z = CGAL::scalar_product(edgeVector, downNormal);
      if (z > steepestEdgeZ)
      {
        steepestEdge = edge;
        steepestEdgeZ = z;
      }

      edge = ++edge;
    } while (edge != vertex->vertex_begin());

    steepestEdges.push_back(steepestEdge);
  }

  // TODO: create MST from steepest edges

  return result;
}

Polyhedron unfoldPolyhedron(const Polyhedron &P,
                            const SpanningTree &spanningTree)
{
  // Find the facet that is facing down most
  Polyhedron::Facet_const_handle downFacet = P.facets_begin();
  FT downFacetZ = 0;
  for (Polyhedron::Facet_const_handle facet = P.facets_begin(); facet != P.facets_end(); ++facet)
  {
    // Calculate and store the normal of the facet
    Vector_3 normal = getFaceNormal(facet);

    FT z = CGAL::scalar_product(normal, Vector_3(0, 0, -1));
    std::cout << "Normal: " << normal << "z: " << z << std::endl;
    if (z > downFacetZ)
    {
      downFacet = facet;
      downFacetZ = z;
    }
  }
  std::cout << "Down facet: " << downFacet->plane() << std::endl;

  // Populate ToDo list with facets that need to be transformed
  std::map<Point_3, Point_3> pointMap;
  std::map<Polyhedron::Facet_const_handle, Aff_transformation_3> transformMap;
  std::deque<Polyhedron::Facet_const_handle> todo = {downFacet};
  transformMap[downFacet] = getFallingDownTransformation(
      downFacet->halfedge()->vertex()->point(),
      downFacet->plane().orthogonal_vector(),
      Plane_3(Point_3(0, 0, 0), Vector_3(0, 0, -1)));

  while (todo.size() > 0)
  {
    // Fetch the current facet
    Polyhedron::Facet_const_handle facet = todo.front();
    todo.pop_front();

    // Transform the points of the current facet
    Polyhedron::Halfedge_const_handle h = facet->halfedge();
    do
    {
      // If the point has not been transformed yet, transform it
      if (pointMap.find(h->vertex()->point()) == pointMap.end())
      {
        // TODO: points need to be duplicated on cut edges
        pointMap[h->vertex()->point()] = transformMap[facet].transform(h->vertex()->point());
      }
      h = h->next();
    } while (h != facet->halfedge());

    // Add the children of the current facet to the ToDo list
    for (Polyhedron::Halfedge_const_handle tree_edge : spanningTree)
    {
      auto other_facet = tree_edge->opposite()->facet();
      if (tree_edge->facet() == facet &&
          transformMap.find(other_facet) == transformMap.end())
      {
        todo.push_back(other_facet);
        // The transformation of the child is the unfolding along the shared edge, followed by the parent's transformation
        transformMap[other_facet] = transformMap[facet] * getFaceUnfoldTransformation(tree_edge);
      }
    }
  }

  // TODO: construct the unfolded polyhedron
  Polyhedron unfolded;
  return unfolded;
}
