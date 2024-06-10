#include "folding.cpp"
#include <gtest/gtest.h>

TEST(CreateRotationTransformationTest, IdentityTest) {
    //setup
    Vector_3 normal_a(1, 0, 0);
    Vector_3 normal_b(1, 0, 0);

    //run DUT
    Aff_transformation_3 transformation = createRotationTransformation(normal_a, normal_b);

    //verify
    ASSERT_EQ(transformation, Aff_transformation_3(CGAL::IDENTITY));
}

TEST(CreateRotationTransformationTest, Simple) {
    //setup
    Vector_3 normal_a(0, 1, 0);
    Vector_3 normal_b(0, 0, 1);

    //run DUT
    Aff_transformation_3 transformation = createRotationTransformation(normal_a, normal_b);

    //verify
    Vector_3 normal_a_transformed = normal_a.transform(transformation);
    ASSERT_LT((normal_a_transformed-=normal_b).squared_length(), 1e-6);
}

TEST(CreateRotationTransformationTest, InvertedTest) {
    //setup
    Vector_3 normal_a(0.137907, -0.721527, 0.678514);
    Vector_3 normal_b(-0.137907, 0.721527, -0.678514);

    //run DUT
    Aff_transformation_3 transformation = createRotationTransformation(normal_a, normal_b);

    //verify
    Vector_3 normal_a_transformed = normal_a.transform(transformation);
    ASSERT_LT((normal_a_transformed-=normal_b).squared_length(), 1e-6);
}

TEST(FallingDownTransformationTest, TypicalTriangle) {
    //setup
    const FT epsilon = 1e-6;
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
    Aff_transformation_3 transformation = getFallingDownTransformation(
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

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
