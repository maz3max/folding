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

        SpanningTree mst;

        Polyhedron U = unfoldPolyhedron(P, mst);

        {
                std::ofstream out("original.wrl");
                CGAL::VRML_2_ostream S(out);
                S << P;
        }

        {
                std::ofstream out("unfolded.wrl");
                CGAL::VRML_2_ostream S(out);
                S << U;
        }

        return EXIT_SUCCESS;
}
