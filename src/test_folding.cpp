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

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
