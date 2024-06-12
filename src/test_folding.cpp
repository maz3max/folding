#include "folding.cpp"
#include <gtest/gtest.h>
#include <algorithm>

static const FT epsilon = 1e-6;

TEST(getRotationTransformationTest, IdentityTest) {
    //setup
    Vector_3 normal_a(1, 0, 0);
    Vector_3 normal_b(1, 0, 0);

    //run DUT
    Aff_transformation_3 transformation = getRotationTransformation(normal_a, normal_b);

    //verify
    ASSERT_EQ(transformation, Aff_transformation_3(CGAL::IDENTITY));
}

TEST(getRotationTransformationTest, Simple) {
    //setup
    Vector_3 normal_a(0, 1, 0);
    Vector_3 normal_b(0, 0, 1);

    //run DUT
    Aff_transformation_3 transformation = getRotationTransformation(normal_a, normal_b);

    //verify
    Vector_3 normal_a_transformed = normal_a.transform(transformation);
    ASSERT_LT((normal_a_transformed-=normal_b).squared_length(), epsilon);
}

TEST(getRotationTransformationTest, InvertedTest) {
    //setup
    Vector_3 normal_a(0.137907, -0.721527, 0.678514);
    Vector_3 normal_b(-0.137907, 0.721527, -0.678514);

    //run DUT
    Aff_transformation_3 transformation = getRotationTransformation(normal_a, normal_b);

    //verify
    Vector_3 normal_a_transformed = normal_a.transform(transformation);
    ASSERT_LT((normal_a_transformed-=normal_b).squared_length(), epsilon);
}

TEST(FallingDownTransformationTest, TypicalTriangle) {
    //setup
    const Point_3 a(-0.525731, 0.850651, 0);
    const Point_3 b(0.525731, 0.850651, 0);
    const Point_3 c(0, 0.525731, -0.850651);
    const Vector_3 normal = getThreePointNormal(a, b, c);
    const FT la = (a-b).squared_length();
    const FT lb = (b-c).squared_length();
    const FT lc = (c-a).squared_length();

    ASSERT_LT(CGAL::abs(la - lb), epsilon);
    ASSERT_LT(CGAL::abs(lb - lc), epsilon);
    ASSERT_LT(CGAL::abs(lc - la), epsilon);

    //run DUT
    const Aff_transformation_3 transformation = getFallingDownTransformation(
        a, // any point on the triangle
        normal, // normal of the triangle
        Plane_3(Point_3(0, 0, 0), Vector_3(0, 0, -1)) // the "down" plane
    );

    const Point_3 a_transform = a.transform(transformation);
    const Point_3 b_transform = b.transform(transformation);
    const Point_3 c_transform = c.transform(transformation);

    //verify
    // the transformed points should be on the "down" plane
    ASSERT_LT(CGAL::abs(a_transform.z()), epsilon);
    ASSERT_LT(CGAL::abs(b_transform.z()), epsilon);
    ASSERT_LT(CGAL::abs(c_transform.z()), epsilon);

    // the distance between the transformed points should be the same
    const FT la_transform = (a_transform - b_transform).squared_length();
    const FT lb_transform = (b_transform - c_transform).squared_length();
    const FT lc_transform = (c_transform - a_transform).squared_length();
    ASSERT_LT(CGAL::abs(la_transform - lb_transform), epsilon);
    ASSERT_LT(CGAL::abs(lb_transform - lc_transform), epsilon);
    ASSERT_LT(CGAL::abs(la_transform - la), epsilon);

    // the new normal should be the same as the plane normal
    const Vector_3 normal_transform = getThreePointNormal(a_transform, b_transform, c_transform);
    ASSERT_LT((normal_transform - Vector_3(0, 0, -1)).squared_length(), epsilon);
}

TEST(RotationAroundLineSegmentTest, IsocahedronFaces)
{
    //setup
    const Point_3 a(-0.525731, 0.850651, 0);
    const Point_3 b(0.525731, 0.850651, 0);
    const Point_3 c(0, 0.525731, -0.850651);
    const Point_3 d(-0.850651, 0, -0.525731);
    // triangle ABC is a face of the isocahedron
    const Vector_3 normal0 = getThreePointNormal(c, a, b);
    // triangle CDA is a face of the isocahedron
    const Vector_3 normal1 = getThreePointNormal(a, c, d);
    const FT la = (a-b).squared_length();
    const FT lb = (b-c).squared_length();
    const FT lc = (c-a).squared_length();
    const FT ld = (d-a).squared_length();

    ASSERT_LT(CGAL::abs(la - lb), epsilon);
    ASSERT_LT(CGAL::abs(lb - lc), epsilon);
    ASSERT_LT(CGAL::abs(lc - la), epsilon);
    ASSERT_LT(CGAL::abs(la - ld), epsilon);

    //run DUT
    const Aff_transformation_3 transform = getRotationAroundLineSegment(normal0, normal1, a);

    //verify
    // transform the point that is not shared by the two faces
    const Point_3 a_transform = a.transform(transform);
    const Point_3 b_transform = b.transform(transform);
    const Point_3 c_transform = c.transform(transform);
    // the transformed point should be on the plane of the second face
    ASSERT_LT(CGAL::squared_distance(Plane_3(a, normal1), b_transform), epsilon);
    // the line segment points should stay the same after transformation
    ASSERT_LT((a_transform - a).squared_length(), epsilon);
    ASSERT_LT((c_transform - c).squared_length(), epsilon);
    // the distance between the transformed points should be the same
    ASSERT_LT(CGAL::abs((a - b_transform).squared_length() - la), epsilon);
    ASSERT_LT(CGAL::abs((c - b_transform).squared_length() - la), epsilon);
}

TEST(FindMaxZVertexTest, Tet)
{
    //setup
    Polyhedron P;
    std::ifstream(CGAL::data_file_path("meshes/tetrahedron.off")) >> P;

    //run DUT
    const Point_3 v0 = findMaxZVertex(P, Vector_3(0, 0, 1))->point();
    const Point_3 v1 = findMaxZVertex(P, Vector_3(1.2, 0, 1.0))->point();

    //verify
    ASSERT_EQ(v0, Point_3(0.0, 0.0, 1.0));
    ASSERT_EQ(v1, Point_3(1.0, 0.0, 0.0));
}

TEST(FindDownFacetTest, Tet)
{
    //setup
    Polyhedron P;
    std::ifstream(CGAL::data_file_path("meshes/tetrahedron.off")) >> P;

    //run DUT
    auto f0 = findDownFacet(P, Vector_3(-1, -1, -1));

    //verify
    auto vertices = getFaceVertices(f0->halfedge());
    ASSERT_EQ(vertices.size(), 3);
    ASSERT_TRUE(std::find(vertices.begin(), vertices.end(), Point_3(1, 0, 0)) != vertices.end());
    ASSERT_TRUE(std::find(vertices.begin(), vertices.end(), Point_3(0, 1, 0)) != vertices.end());
    ASSERT_TRUE(std::find(vertices.begin(), vertices.end(), Point_3(0, 0, 1)) != vertices.end());
}

TEST(GetSteepestEdgesTest, Cube)
{
    //setup
    Polyhedron P;
    std::ifstream(CGAL::data_file_path("meshes/cube_quad.off")) >> P;
    Vector_3 normal(0.64486, 0.324763, 0.691871);

    //run DUT
    auto maxZVertex = findMaxZVertex(P, normal);
    auto steepestEdges = getSteepestEdges(P, normal, maxZVertex);
    auto startFacet = findDownFacet(P, normal);

    for (auto edge: steepestEdges)
    {
        std::cout << "Edge: (" << edge.first << "), (" << edge.second << ")" << std::endl;
    }
    auto tree = constructSpanningTree(P, startFacet, steepestEdges);
    for (auto[parent, face_vertices] : tree.children) {
        std::cout << "[" << parent << "] -> ";
        for (auto vertex: face_vertices) {
            std::cout << "(" << vertex << "), ";
        }
        std::cout << std::endl;
    }
    //TODO
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
