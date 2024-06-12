#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/VRML.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

#include "folding.cpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

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
        std::ifstream in1((argc > 1) ? argv[1] : CGAL::data_file_path("meshes/icosahedron.off"));
        in1 >> P;

        Vector_3 normal(0.64486, 0.324763, 0.691871);
        normal = normal / CGAL::approximate_sqrt(normal.squared_length());

        // print the normal vector
        std::cout << "Normal vector: " << normal << std::endl;

        return EXIT_SUCCESS;
}
