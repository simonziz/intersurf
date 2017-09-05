#include "includes.h"
#include "typedefs.cpp"


Delaunay makeInterfaceTriangulation(PdbImage* pdb, std::map<unsigned int, bool> is_interface, Delaunay d_t);

K::Point_3 calculateBarycenter(Delaunay interface_tr);

double calculateMaxDist(Delaunay interface_tr, K::Point_3 barycenter);

void calculateInterface(PdbImage *pdb, Delaunay interface_tr, K::Point_3 barycenter, double squared_max_dist,
                        std::map<K::Point_3, unsigned int> *surf_points, std::vector<std::vector<K::Point_3> > *all_faces_indexes);
