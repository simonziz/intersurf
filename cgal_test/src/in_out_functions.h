#include "includes.h"
#include "typedefs.cpp"


void savePointsOFF(const char* filename, Delaunay m_dt, PdbImage *pdb);

void saveSurface0FF(const char* filename, std::map<K::Point_3, unsigned int > surf_points, std::vector<std::vector<K::Point_3> > all_faces_indexes);

Delaunay readPdb(PdbImage *pdb, std::map<unsigned int, bool> is_interface);

Polyhedron makePolyhedron(const char* filename, bool show_warnings);
