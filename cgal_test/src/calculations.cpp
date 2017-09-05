#include "calculations.h"


/// Function to create a triangulation containing only the tetrahedra at the interface

Delaunay makeInterfaceTriangulation(PdbImage* pdb, std::map<unsigned int, bool> is_interface, Delaunay d_t)
{
    bool is_relevant = false;
    std::vector<std::pair< Point, unsigned> > reduced_vector;

    //iterating on the cells
    for (Delaunay::Finite_cells_iterator cit = d_t.finite_cells_begin(); cit != d_t.finite_cells_end(); ++cit) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                //checking if a cell contains 2 different chains
                if (strcmp(pdb->atom[cit->vertex(j)->info()].chain, pdb->atom[cit->vertex(k)->info()].chain) != 0) {
                    //std::cout<<strcmp(pdb->atom[cit->vertex(j)->info()].chain, pdb->atom[cit->vertex(k)->info()].chain)<<std::endl;
                    is_relevant = true; // if yes -> to be kept
                    //std::cout<<pdb->atom[cit->vertex(j)->info()].chain<<pdb->atom[cit->vertex(k)->info()].chain<<std::endl;
                }
            }
        }
        if (is_relevant == true) {
            for (int j = 0; j < 4; j++) {
                is_interface[cit->vertex(j)->info()] = true;
            }
        }
    }
    // Filling the vector with the coordinates of the atoms
    for (int i = 0; i < pdb->n_atoms; i++) {
        if (is_interface[i] == true) {
            reduced_vector.push_back(std::make_pair(Point(pdb->atom[i].pt.x, pdb->atom[i].pt.y, pdb->atom[i].pt.z), i) );
        }
    }

    // Creating Delaunay triangulation with the coordinates vector
    Delaunay interface_tr(reduced_vector.begin(), reduced_vector.end());

    return interface_tr;
}


/// Function to calculate the barycenter of the interface triangulation

K::Point_3 calculateBarycenter(Delaunay interface_tr)
{
    // Barycenter
    double x_coord = 0.0;
    double y_coord = 0.0;
    double z_coord = 0.0;
    for (Delaunay::Finite_vertices_iterator vit = interface_tr.finite_vertices_begin(); vit != interface_tr.finite_vertices_end(); ++vit) {
        x_coord += vit->point().x() / interface_tr.number_of_vertices();
        y_coord += vit->point().y() / interface_tr.number_of_vertices();
        z_coord += vit->point().z() / interface_tr.number_of_vertices();
    }
    K::Point_3 barycenter(x_coord,y_coord,z_coord);

    return barycenter;
}


/// Function to calculate the maximum distance from the barycenter
double calculateMaxDist(Delaunay interface_tr, K::Point_3 barycenter)
{
    // Max distance
    K::Vector_3 max_dist_vec; // 3D vector : distance to be normalized
    double squared_max_dist = 0.0;
    for (Delaunay::Finite_vertices_iterator vit = interface_tr.finite_vertices_begin(); vit != interface_tr.finite_vertices_end(); ++vit) {
        max_dist_vec = vit->point() - barycenter;
        if (max_dist_vec.squared_length() > squared_max_dist) {
            squared_max_dist = max_dist_vec.squared_length();
        }
    }

    return squared_max_dist;
}


/// Function to calculate the interface (2 arrays) from a given triangulation

void calculateInterface(PdbImage *pdb, Delaunay interface_tr, K::Point_3 barycenter, double squared_max_dist,
                        std::map<K::Point_3, unsigned int> *surf_points, std::vector<std::vector<K::Point_3> > *all_faces_indexes)
{
    // Initialize circulators
    Delaunay::Cell_circulator circ; // Circulator for cells around a given edge
    Delaunay::Cell_circulator circ_copy; // Copy to check when all cells have been tested

    // Initialize variables
    std::vector<K::Point_3> face_indexes; // points of one face
    K::Point_3 p; // 3D point : in the triangulation
    K::Vector_3 dist; // 3D vector : distance to be normalized
    bool is_used = true; // To check if the point is not too far from the barycenter
    bool is_finite_tetrahedron = true; // To check if a tetrahedron contains the infinite point
    //double size_lim = 34 * 34; // Max distance from the barycenter    ***** To be changed *****
    double size_lim = squared_max_dist + 0.3 * squared_max_dist; // Max distance from the barycenter

    // Iterate on the edges of the whole triangulation
    for (Delaunay::Finite_edges_iterator eit = interface_tr.finite_edges_begin(); eit != interface_tr.finite_edges_end(); ++eit) {
        // Check if the edge is at the interface
        if((strcmp(pdb->atom[eit->first->vertex(eit->second)->info()].chain, pdb->atom[eit->first->vertex(eit->third)->info()].chain) != 0)){
            //std::cout<<pdb->atom[eit->first->vertex(eit->second)->info()].chain<<"   "<<pdb->atom[eit->first->vertex(eit->third)->info()].chain<<std::endl;
            circ = interface_tr.incident_cells(*eit);
            circ_copy = circ;
            for (int i = 0; i < 4; i++) {
                if(interface_tr.is_infinite(circ->vertex(i))) { //check if the vertex is the infinite one
                    is_finite_tetrahedron = false; //if yes   the tetrahedron will not be stored
                    //std::cout<<circ->vertex(i)->point()<<std::endl; //display infinite point
                }
            }
            do {
                p = interface_tr.dual(circ); // Find the dual of the cell (Voronoï Diagram)
                // Checking the distance between the barycenter and the point is in the limit
                if ((barycenter - p).squared_length() <= size_lim) {
                    surf_points->insert(std::make_pair(p,0)); // Store the point in the map with 0 value for the index
                }
                face_indexes.push_back(p); // Put the point in the vector of the current face
                //std::cout<<"face : "<<face_indexes.empty()<<std::endl;
                // Choix du sens de circulation selon l'orientaion de l'arête (A->B ou B->A)
                if (strcmp(pdb->atom[eit->first->vertex(eit->second)->info()].chain, "A") == 0) {
                    circ ++;
                }
                else {
                    circ --;
                }
            } while(circ != circ_copy); // Circle around the edge until the first cell is reached again
            for (int j = 0; j < face_indexes.size(); j++) { // Go through the points of the current face
                dist = barycenter - face_indexes[j]; // calculate the distance between the barycenter and the current point
                if ( dist.squared_length() > size_lim) { // If the distance of one point is over the limit
                    is_used = false; // we will not store the face
                }
            }
            // Check is the face is in the limit and does not contain the infinte point
            if (is_used == true && is_finite_tetrahedron == true) {
                all_faces_indexes->push_back(face_indexes); // Add that face to the intersurf vector
            }
            is_used = true;
            is_finite_tetrahedron = true; // Put the booleans back to true (at the begining the face is assumed to be valid)
            face_indexes.clear(); // Clear the current face vector for the next iteration
        }
    }

    //std::cout<<"point : "<<surf_points->empty()<<std::endl;
    //std::cout<<"faces : "<<all_faces_indexes->empty()<<std::endl;
}


/// Function to smooth the surface with the Catmull/Clark subdivision algorithm


