//#ifndef MY_VERTEX_BASE
//#define MY_VERTEX_BASE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>


template < class GT, class Vb = CGAL::Triangulation_vertex_base_3<GT> >
class my_vertex_base
  : public Vb
{
public:
  typedef typename Vb::Point          Point;
  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef my_vertex_base<GT, Vb2>                        Other;
  };
  my_vertex_base() {}
  my_vertex_base(const Point& p)
    : Vb(p) {}
  my_vertex_base(const Point& p, int index)
    : Vb(p, index) {}

  int index;

};

