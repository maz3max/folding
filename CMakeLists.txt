cmake_minimum_required (VERSION 3.0)
project (myproject)

find_package(CGAL 5.4 REQUIRED OPTIONAL_COMPONENTS Qt5)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable (folding src/folding.cpp)
target_link_libraries (folding CGAL::CGAL Eigen3::Eigen)

if(CGAL_Qt5_FOUND)
  target_link_libraries(folding CGAL::CGAL_Basic_viewer)
endif()
