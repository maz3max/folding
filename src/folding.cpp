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

template <class T>
struct SimpleTree
{
  SimpleTree *parent;
  std::vector<SimpleTree> children;
  T data;
};

template <class T>
struct SimpleTree2
{
  std::vector<std::pair<size_t, T>> children;
};

typedef SimpleTree<Polyhedron::Halfedge_const_handle> SpanningTree;

Aff_transformation_3 getRotationTransformation(const Vector_3 &normal_a, const Vector_3 &normal_b)
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
      getRotationTransformation(facetNormal, planeNormal) *
      Aff_transformation_3(CGAL::TRANSLATION, -facetPointVector);

  // Combine translation and rotation
  return translate * rotate;
}

/** Compute Normal from three points of a face in ccw order */
Vector_3 getThreePointNormal(const Point_3 &a, const Point_3 &b, const Point_3 &c)
{
  Vector_3 normal = CGAL::cross_product(b - a, c - b);
  normal /= CGAL::approximate_sqrt(normal.squared_length());
  return normal;
}

/** Compute Normal of a polyhedron facet */
Vector_3 getFaceNormal(const Polyhedron::Facet_const_handle facet)
{
  Polyhedron::Halfedge_const_handle h = facet->halfedge();
  return getThreePointNormal(
    h->vertex()->point(),
    h->next()->vertex()->point(),
    h->next()->next()->vertex()->point()
  );
}

Vector_3 getFaceNormal(const std::vector<Point_3> &face)
{
  return getThreePointNormal(face[0], face[1], face[2]);
}

Vector_3 getFaceNormal(const std::vector<size_t> &face, const std::vector<Point_3> &vertices)
{
  return getThreePointNormal(vertices[face[0]], vertices[face[1]], vertices[face[2]]);
}

// generate transfrom from one plane to another by rotating around a line segment
Aff_transformation_3 getRotationAroundLineSegment(const Vector_3& currentNormal, const Vector_3& targetNormal, const Point_3& edgePoint)
{
  Vector_3 edgePointVector = Vector_3(edgePoint.x(), edgePoint.y(), edgePoint.z());
  return Aff_transformation_3(CGAL::TRANSLATION, edgePointVector) *
         getRotationTransformation(currentNormal, targetNormal) *
         Aff_transformation_3(CGAL::TRANSLATION, -edgePointVector);
}

// generate transform to unfold the face on the opposing edge along the edges to match the face along the current edge
Aff_transformation_3 getFaceUnfoldTransformation(const Polyhedron::Halfedge_const_handle edge)
{
  Vector_3 targetNormal = getFaceNormal(edge->facet());
  Vector_3 currentNormal = getFaceNormal(edge->opposite()->facet());
  Point_3 edgePoint = edge->vertex()->point();

  return getRotationAroundLineSegment(currentNormal, targetNormal, edgePoint);
}

Polyhedron::Vertex_const_handle getFurthestVertex(const Polyhedron &P, const Vector_3 downNormal)
{
  Polyhedron::Vertex_const_handle maxZVertex = P.vertices_begin();
  FT maxZ = -INFINITY;
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
  return maxZVertex;
}

std::vector<std::pair<Point_3, Point_3>> getSteepestEdges(const Polyhedron &P, const Vector_3 downNormal, Polyhedron::Vertex_const_handle maxZVertex)
{
  std::vector<std::pair<Point_3, Point_3>> steepestEdges;
  for (Polyhedron::Vertex_const_handle vertex = P.vertices_begin(); vertex != P.vertices_end(); ++vertex)
  {
    if (vertex == maxZVertex)
    {
      continue;
    }
    Polyhedron::Halfedge_const_handle steepestEdge = vertex->halfedge();
    FT steepestEdgeZ = -INFINITY;

    // iterate over edges pointing to vertex
    Polyhedron::Halfedge_around_vertex_const_circulator edge = vertex->vertex_begin();
    do
    {
      // normalized vector pointing away from the vertex
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

    steepestEdges.emplace_back(steepestEdge->vertex()->point(), steepestEdge->opposite()->vertex()->point());
  }
  return steepestEdges;
}

Polyhedron::Facet_const_handle findDownFacet(const Polyhedron &P, const Vector_3 downNormal)
{
  Polyhedron::Facet_const_handle downFacet = P.facets_begin();
  FT downFacetZ = -INFINITY;
  for (Polyhedron::Facet_const_handle facet = P.facets_begin(); facet != P.facets_end(); ++facet)
  {
    // Calculate and store the normal of the facet
    Vector_3 normal = getFaceNormal(facet);

    FT z = CGAL::scalar_product(normal, downNormal);
    if (z > downFacetZ)
    {
      downFacet = facet;
      downFacetZ = z;
    }
  }
  return downFacet;
}

std::vector<Point_3> getFaceVertices(const Polyhedron::Halfedge_const_handle commonEdge)
{
  Polyhedron::Halfedge_const_handle edge = commonEdge;
  std::vector<Point_3> faceVertices;
  do
  {
    faceVertices.push_back(edge->vertex()->point());
    edge = edge->next();
  } while (edge != commonEdge);
  return faceVertices;
}

std::vector<size_t> getFaceVertexIndices(const Polyhedron::Halfedge_const_handle commonEdge,
const std::vector<Point_3> &vertices)
{
  Polyhedron::Halfedge_const_handle edge = commonEdge;
  std::vector<size_t> faceVertices;
  do
  {
    auto index = std::find(vertices.begin(), vertices.end(), edge->vertex()->point()) - vertices.begin();
    faceVertices.push_back(index);
    edge = edge->next();
    assert(index < vertices.size());
  } while (edge != commonEdge);
  return faceVertices;
}

std::pair<SimpleTree2<std::vector<size_t>>, std::vector<Point_3>> constructSpanningTree(
  const Polyhedron &P,
  Polyhedron::Facet_const_handle startFace,
  std::vector<std::pair<Point_3, Point_3>> excludedEdges
)
{
  // create a list of all vertices
  std::vector<Point_3> vertices;
  for (auto &vertex : P.vertex_handles())
  {
    vertices.push_back(vertex->point());
  }

  SimpleTree2<std::vector<size_t>> tree;
  // todo list contains the index of the current node and the facet
  std::deque<std::pair<size_t, Polyhedron::Facet_const_handle>> todo{{0, startFace}};
  // visited facets
  std::set<Polyhedron::Facet_const_handle> visited{startFace};
  // add the first face to the tree (parent index, facet) (root is pointing to itself)
  tree.children.push_back({0, getFaceVertexIndices(startFace->halfedge(), vertices)});


  while (todo.size() > 0)
  {
    auto [current_index, currentFace] = todo.front();
    todo.pop_front();

    // Add the neighboring faces to the todo list
    Polyhedron::Halfedge_const_handle edgeIterator = currentFace->halfedge();
    do
    {
      // Skip edge if both points are in the excluded list
      Point_3 point1 = edgeIterator->vertex()->point();
      Point_3 point2 = edgeIterator->opposite()->vertex()->point();
      if (
        std::find(excludedEdges.begin(), excludedEdges.end(), std::make_pair(point1, point2)) != excludedEdges.end()
        || std::find(excludedEdges.begin(), excludedEdges.end(), std::make_pair(point2, point1)) != excludedEdges.end()
      )
      {
        edgeIterator = edgeIterator->next();
        continue;
      }

      // Add the neighbor to the todo list
      Polyhedron::Facet_const_handle neighbor = edgeIterator->opposite()->facet();
      if (visited.find(neighbor) == visited.end())
      {
        tree.children.push_back({current_index, getFaceVertexIndices(edgeIterator->opposite()->prev(), vertices)});
        visited.insert(neighbor);
        todo.push_back({tree.children.size() - 1, neighbor});
      }
      edgeIterator = edgeIterator->next();
    } while (edgeIterator != currentFace->halfedge());
  }
  return {tree, vertices};
}

// The steepest edge cut tree contains the steepest edges of all vertices, except the vertex with maximal z-coordinate.
std::pair<SimpleTree2<std::vector<size_t>>, std::vector<Point_3>> steepestEdgeCut(const Polyhedron &P, const Vector_3 normal)
{
  // Find the vertex with maximal z-coordinate
  auto maxZVertex = getFurthestVertex(P, -normal);

  // Find the steepest edge for each vertex
  auto steepestEdges = getSteepestEdges(P, normal, maxZVertex);

  // Find the facet that is facing down most
  auto downFacet = findDownFacet(P, Vector_3(0, 0, -1));

  // Construct the spanning tree
  return constructSpanningTree(P, downFacet, steepestEdges);
}

std::vector<Point_3> unfoldTree(const SimpleTree2<std::vector<size_t>> &tree, const std::vector<Point_3> &vertices)
{
  Aff_transformation_3 transformations[tree.children.size()];
  std::vector<Point_3> transformedVertices(vertices.size());

  // set up initial transformation
  transformations[0] = getFallingDownTransformation(
    vertices[tree.children[0].second[0]],
    getFaceNormal(tree.children[0].second, vertices),
    Plane_3(Point_3(0, 0, 0), Vector_3(0, 0, -1))
  );
  // set up initial face
  for (size_t point_idx : tree.children[0].second)
  {
    transformedVertices[point_idx] = vertices[point_idx].transform(transformations[0]);
  }

  for (size_t i = 1; i< tree.children.size(); i++)
  {
    // The parent of the current node is the parent of the parent of the current node
    size_t parent_idx = tree.children[i].first;
    auto &face = tree.children[i].second;
    Aff_transformation_3 &parent_transform = transformations[parent_idx];
    transformations[i] = parent_transform * getRotationAroundLineSegment(
      getFaceNormal(face, vertices),
      getFaceNormal(tree.children[parent_idx].second, vertices),
      vertices[face[0]] // edge point, first two points are parent edge
    );
    // transform non-parent points
    for (size_t j = 2; j < face.size(); j++)
    {
      size_t point_idx = face[j];
      transformedVertices[point_idx] = vertices[point_idx].transform(transformations[i]);
    }
  }
  return transformedVertices;
}



std::string treeToSVG(const SimpleTree2<std::vector<size_t>> &tree, const std::vector<Point_3> &vertices)
{
  // extract crease edges
  std::vector<std::pair<size_t, size_t>> creaseEdges;
  for (size_t i = 1; i < tree.children.size(); i++)
  {
    size_t parent_idx = tree.children[i].first;
    auto &face = tree.children[i].second;
    if (face[0] < face[1]) {
      creaseEdges.push_back({face[0], face[1]});
    } else {
      creaseEdges.push_back({face[1], face[0]});
    }
  }
  // sort crease edges
  std::sort(creaseEdges.begin(), creaseEdges.end());

  std::vector<std::pair<size_t, size_t>> outlineEdges;
  for (auto &[parent, face] : tree.children)
  {
    for (size_t i = 0; i < face.size(); i++)
    {
      size_t point1 = face[i];
      size_t point2 = face[(i + 1) % face.size()];
      auto pair = (point1 < point2) ? std::make_pair(point1, point2) : std::make_pair(point2, point1);
      if (!std::binary_search(creaseEdges.begin(), creaseEdges.end(), pair))
      {
        outlineEdges.push_back(std::make_pair(point1, point2));
      }
    }
  }

  for (auto &[point1, point2] : outlineEdges)
  {
    std::cout << "Outline edge: (" << point1 << ") -> (" << point2  << ")" << std::endl;
  }

  std::vector<size_t> outline;
  // start the outline with a leaf node
  outline.push_back(outlineEdges[0].first);
  outline.push_back(outlineEdges[0].second);
  outlineEdges.erase(outlineEdges.begin());
  while (outlineEdges.size() > 0)
  {
    for (size_t i = 0; i < outlineEdges.size(); i++)
    {
      if (outline.back() == outlineEdges[i].first)
      {
        std::cout << "Adding edge " << outlineEdges[i].first << " -> " << outlineEdges[i].second << std::endl;
        outline.push_back(outlineEdges[i].second);
        outlineEdges.erase(outlineEdges.begin() + i);
        break;
      }
    }
  }


  std::string result = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"200\">\n";
  // add outline
  result += "<polygon points=\"";
  // for each point in the outline, add it to the svg as "x,y "
  for (auto &point : outline)
  {
    result += std::to_string(vertices[point].x()) + "," + std::to_string(vertices[point].y()) + " ";
  }
  result += "\" fill=\"none\" stroke=\"red\" stroke-width=\"1\"/>\n";
  for (auto &[point1, point2] : creaseEdges)
  {
    result += "<line x1=\"" + std::to_string(vertices[point1].x()) + "\" y1=\"" + std::to_string(vertices[point1].y()) + "\" x2=\"" + std::to_string(vertices[point2].x()) + "\" y2=\"" + std::to_string(vertices[point2].y()) + "\" stroke=\"black\" stroke-width=\"1\"/>\n";
  }
  // add crease edges
  result += "</svg>";
  return result;
}

std::string treeToSVG2(const SimpleTree2<std::vector<size_t>> &tree, const std::vector<Point_3> &vertices)
{

  std::vector<std::pair<size_t, size_t>> outlineEdges;
  for (auto &[parent, face] : tree.children)
  {
    for (size_t i = 0; i < face.size(); i++)
    {
      auto &point1 = face[i];
      auto &point2 = face[(i + 1) % face.size()];
      outlineEdges.push_back(std::make_pair(point1, point2));
    }
  }


  std::string result = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"200\">\n";
  // add outline
  for (auto &[point1, point2] : outlineEdges)
  {
    result += "<line x1=\"" + std::to_string(vertices[point1].x()) + "\" y1=\"" + std::to_string(vertices[point1].y()) + "\" x2=\"" + std::to_string(vertices[point2].x()) + "\" y2=\"" + std::to_string(vertices[point2].y()) + "\" stroke=\"red\" stroke-width=\"0.05\"/>\n";
  }
  // add crease edges
  result += "</svg>";
  return result;
}
