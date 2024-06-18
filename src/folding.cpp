#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Aff_transformation_3.h>
#include <algorithm>
#include <iostream>
#include <cstdlib>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::FT FT;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polygon_2<Kernel> Polygon;

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

SimpleTree2<std::vector<Point_3>> constructSpanningTree(
  const Polyhedron &P,
  Polyhedron::Facet_const_handle startFace,
  std::vector<std::pair<Point_3, Point_3>> excludedEdges
)
{
  SimpleTree2<std::vector<Point_3>> tree;
  // todo list contains the index of the current node and the facet
  std::deque<std::pair<size_t, Polyhedron::Facet_const_handle>> todo{{0, startFace}};
  // visited facets
  std::set<Polyhedron::Facet_const_handle> visited{startFace};
  // add the first face to the tree (parent index, facet) (root is pointing to itself)
  tree.children.push_back({0, getFaceVertices(startFace->halfedge())});

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
        tree.children.push_back({current_index, getFaceVertices(edgeIterator->opposite()->prev())});
        visited.insert(neighbor);
        todo.push_back({tree.children.size() - 1, neighbor});
      }
      edgeIterator = edgeIterator->next();
    } while (edgeIterator != currentFace->halfedge());
  }
  return tree;
}

bool doPolygonsIntersect(const Polygon &a, const Polygon &b)
{
  for (size_t i = 0; i < a.size(); i++)
  {
    for (size_t j = 0; j < b.size(); j++)
    {
      if (CGAL::do_intersect(a.edge(i), b.edge(j)))
      {
        return true;
      }
    }
  }
  return false;
}

bool treeContainsIntersection(const SimpleTree2<std::vector<Point_3>> &tree)
{
    std::vector<Polygon> faces;
    for (auto &[parent, face] : tree.children)
    {
      Polygon polygon;
      for (auto &point : face)
      {
        polygon.push_back(Point_2(point.x(), point.y()));
      }
      faces.push_back(polygon);
    }
    // Check for intersections
    for(size_t i = 0; i < faces.size(); i++)
    {
      for(size_t j = i + 1; j < faces.size(); j++)
      {
        if (tree.children[i].first == j || tree.children[j].first == i)
        {
          /* adjacent faces are allowed to intersect */
          continue;
        }
        if (doPolygonsIntersect(faces[i], faces[j]))
        {
          return true;
        }
      }
    }
    return false;
}

SimpleTree2<std::vector<Point_3>> unfoldTree(const SimpleTree2<std::vector<Point_3>> &tree)
{
  SimpleTree2<std::vector<Point_3>> unfoldedTree; unfoldedTree.children.reserve(tree.children.size());
  Aff_transformation_3 transformations[tree.children.size()];

  // set up initial transformation
  transformations[0] = getFallingDownTransformation(
    tree.children[0].second[0],
    getFaceNormal(tree.children[0].second),
    Plane_3(Point_3(0, 0, 0), Vector_3(0, 0, -1))
  );
  // set up initial face
  unfoldedTree.children.push_back({0, {}});
  unfoldedTree.children[0].second.reserve(tree.children[0].second.size());
  for (auto &point : tree.children[0].second)
  {
    unfoldedTree.children[0].second.push_back(point.transform(transformations[0]));
  }

  for (size_t i = 1; i< tree.children.size(); i++)
  {
    // The parent of the current node is the parent of the parent of the current node
    size_t parent_idx = tree.children[i].first;
    auto &face = tree.children[i].second;
    Aff_transformation_3 &parent_transform = transformations[parent_idx];
    transformations[i] = parent_transform * getRotationAroundLineSegment(
      getFaceNormal(face),
      getFaceNormal(tree.children[parent_idx].second),
      face[0] // edge point, first two points are parent edge
    );
    unfoldedTree.children.push_back({parent_idx, {}});
    auto &transformed_face = unfoldedTree.children[i].second;
    transformed_face.reserve(face.size());
    // transform parent edge (double computation, but saves searching for the edge)
    transformed_face.push_back(face[0].transform(parent_transform));
    transformed_face.push_back(face[1].transform(parent_transform));
    // transform remaining points
    for (size_t j = 2; j < face.size(); j++)
    {
      transformed_face.push_back(face[j].transform(transformations[i]));
    }
  }
  return unfoldedTree;
}



// The steepest edge cut tree contains the steepest edges of all vertices, except the vertex with maximal z-coordinate.
std::pair<SimpleTree2<std::vector<Point_3>>,SimpleTree2<std::vector<Point_3>>> steepestEdgeCut(const Polyhedron &P)
{
  SimpleTree2<std::vector<Point_3>> tree;
  SimpleTree2<std::vector<Point_3>> unfolded;
  while (true) {
    // generate new random normal vector
    Vector_3 normal = Vector_3(std::rand(), std::rand(), std::rand());
    normal /= CGAL::approximate_sqrt(normal.squared_length());

    // Find the vertex with maximal z-coordinate
    auto maxZVertex = getFurthestVertex(P, -normal);

    // Find the steepest edge for each vertex
    auto steepestEdges = getSteepestEdges(P, normal, maxZVertex);

    // Find the facet that is facing down most
    auto downFacet = findDownFacet(P, Vector_3(0, 0, -1));
    // Construct the spanning tree
    tree = constructSpanningTree(P, downFacet, steepestEdges);
    // Unfold the tree
    unfolded = unfoldTree(tree);

    if (0)//(treeContainsIntersection(unfolded))
    {
      std::cout << "Found intersection, retrying" << std::endl;
    } else {
      break;
    }
  }
  return {tree, unfolded};
}

SimpleTree2<std::vector<size_t>> convertTreeToIndexed(const Polyhedron &P, const SimpleTree2<std::vector<Point_3>> &tree)
{
  SimpleTree2<std::vector<size_t>> indexedTree;
  std::vector<Point_3> vertices;
  for (const Point_3 &vertex : P.points())
  {
    vertices.push_back(vertex);
  }
  for (auto &[parent, face] : tree.children)
  {
    indexedTree.children.push_back({parent, {}});
    for (const Point_3 &vertex : face)
    {
      size_t index = std::find(vertices.begin(), vertices.end(), vertex) - vertices.begin();
      assert(index < vertices.size());
      indexedTree.children[indexedTree.children.size() - 1].second.push_back(index);
    }
  }
  return indexedTree;
}

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<size_t>> findOutlineIndexed(
  const SimpleTree2<std::vector<size_t>> &tree
)
{
  // crease edges that are not cut and therefore not part of the outline (pairs of point indices)
  std::vector<std::pair<size_t, size_t>> creaseEdges;
  // outline edges (pairs of point indices)
  std::vector<std::pair<size_t, size_t>> outlineEdges;
  // pairs of indices: face index and index within face to define outline
  std::vector<std::pair<size_t, size_t>> outline;
  // index of the other edge in outline that was originally the same edge
  std::vector<size_t> outlineCorrespondingEdge;

  // find crease edges
  for(size_t i = 1; i < tree.children.size(); i++)
  {
    size_t parent_idx = tree.children[i].first;
    auto &face = tree.children[i].second;
    creaseEdges.push_back((face[0] < face[1]) ? std::make_pair(face[0], face[1]) : std::make_pair(face[1], face[0]));
  }
  // sort crease edges
  std::sort(creaseEdges.begin(), creaseEdges.end());

  // find outline
  for (auto &[parent, face] : tree.children)
  {
    for (size_t i = 0; i < face.size(); i++)
    {
      size_t current = face[i];
      size_t next = face[(i + 1) % face.size()];
      std::pair<size_t, size_t> edge = (current < next) ? std::make_pair(face[0], face[1]) : std::make_pair(face[1], face[0]);
      if (!std::binary_search(creaseEdges.begin(), creaseEdges.end(), edge))
      {
        outlineEdges.push_back({current, next});
      }
    }
  }
  //print outline edges
  for (auto [first, second] : outlineEdges)
  {
    std::cout << first << " " << second << std::endl;
  }

  // sort outline
  for (size_t i = 0; i < outlineEdges.size() - 2; i++)
  {
    if (outlineEdges[i].second == outlineEdges[i + 1].first)
    {
      continue;
    }
    for (size_t j = i + 2; j < outlineEdges.size(); j++)
    {
      if (outlineEdges[i].second == outlineEdges[j].first)
      {
        std::swap(outlineEdges[i + 1], outlineEdges[j]);
        break;
      }
    }
  }

  //print outline edges
  for (auto [first, second] : outlineEdges)
  {
    std::cout << first << " " << second << std::endl;
  }

  // find pairs of outline edges that were originally the same edge
  outlineCorrespondingEdge.resize(outlineEdges.size(), -1);
  for (size_t i = 0; i < outlineEdges.size(); i++)
  {
    for (size_t j = 0; j < outlineEdges.size(); j++)
    {
      if (outlineEdges[i].first == outlineEdges[j].second && outlineEdges[i].second == outlineEdges[j].first)
      {
        outlineCorrespondingEdge[i] = j;
        outlineCorrespondingEdge[j] = i;
      }
    }
  }
  // convert outline edges to outline
  for (auto [first, second] : outlineEdges)
  {
    for (size_t i = 0; i < tree.children.size(); i++)
    {
      for (size_t j = 0; j < tree.children[i].second.size(); j++)
      {
        if (tree.children[i].second[j] == first
            && tree.children[i].second[(j + 1) % tree.children[i].second.size()] == second)
        {
          outline.push_back({i, j});
        }
      }
    }
  }

  return {outline, outlineCorrespondingEdge};
}



std::string treeToSVG(const Polyhedron &P, const SimpleTree2<std::vector<Point_3>> &tree, const SimpleTree2<std::vector<Point_3>> &unfolded)
{
  // extract crease edges
  std::vector<std::pair<Point_3, Point_3>> creaseEdges;
  for (size_t i = 1; i < unfolded.children.size(); i++)
  {
    size_t parent_idx = unfolded.children[i].first;
    auto &face = unfolded.children[i].second;
    creaseEdges.push_back({face[0], face[1]});
  }

  auto [outline, outlineCorrespondingEdge] = findOutlineIndexed(convertTreeToIndexed(P, tree));

  std::string result = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"200\">\n";
  // add outline
  result += "<polygon points=\"";
  // for each point in the outline, add it to the svg as "x,y "
  for (auto [faceIdx, pointIdx] : outline)
  {
    Point_3 point = unfolded.children[faceIdx].second[pointIdx];
    result += std::to_string(point.x()) + "," + std::to_string(point.y()) + " ";
  }
  result += "\" fill=\"none\" stroke=\"red\" stroke-width=\"0.05\"/>\n";
  for (auto &[point1, point2] : creaseEdges)
  {
    result += "<line x1=\"" + std::to_string(point1.x()) + "\" y1=\"" + std::to_string(point1.y()) + "\" x2=\"" + std::to_string(point2.x()) + "\" y2=\"" + std::to_string(point2.y()) + "\" stroke=\"black\" stroke-width=\"0.05\"/>\n";
  }
  // add crease edges
  result += "</svg>";
  return result;
}

std::string treeToSVG2(const SimpleTree2<std::vector<Point_3>> &tree)
{
  std::string result = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"200\">\n";
  for (auto &[parent, face] : tree.children)
  {
    result += "<polygon points=\"";
    for (auto &point : face)
    {
      result += std::to_string(point.x()) + "," + std::to_string(point.y()) + " ";
    }
    result += "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.05\"/>\n";
  }
  result += "</svg>";
  return result;
}
