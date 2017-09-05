
#include "includes.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K>    Vb;
//typedef CGAL::Triangulation_data_structure_3< my_vertex_base<K> >           Tds;
typedef CGAL::Triangulation_data_structure_3< Vb >                      Tds;
//typedef CGAL::Delaunay_triangulation_3<K,Tds>        Delaunay;
typedef CGAL::Delaunay_triangulation_3<K,Tds,CGAL::Fast_location>        Delaunay;

typedef Delaunay::Point                          Point;
typedef Delaunay::Cell_handle                    Cell_handle;
typedef Delaunay::Vertex_handle                  Vertex_handle;
typedef Delaunay::Facet                          Face;

typedef CGAL::Polyhedron_3<K>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;
