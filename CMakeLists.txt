cmake_minimum_required (VERSION 3.0)
project (myproject)

find_package(CGAL 5.4 REQUIRED) # OPTIONAL_COMPONENTS Qt5)

file(GLOB SOURCES "src/main.cpp")
add_executable (folding ${SOURCES})
target_link_libraries (folding CGAL::CGAL)

# if(CGAL_Qt5_FOUND)
#   target_link_libraries(folding CGAL::CGAL_Basic_viewer)
# endif()

# add CGAL_DATA_DIR inside project directory
target_compile_definitions(folding PRIVATE CGAL_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

# add gtest tests
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable (test_folding src/test_folding.cpp)
target_link_libraries (test_folding ${GTEST_LIBRARIES} CGAL::CGAL)
add_test (NAME test_folding COMMAND test_folding)

target_compile_definitions(test_folding PRIVATE CGAL_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
