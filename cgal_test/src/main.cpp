#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/IO/Writer_OFF.h>

#include <CGAL/Unique_hash_map.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "pdbs.h"

#include <iostream>
#include <fstream>
#include <ostream>
#include <cstring>
#include <vector>
#include <map>



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


typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;


static void savePointsOFF(const char* filename, Delaunay m_dt, PdbImage *pdb);

int main(void) {


    /// Read a .pdb file into an array

    //PdbImage *pdb = hex_readPdb("../data/toto.pdb", "new_protein");  // Reading a .pdb file
    PdbImage *pdb = hex_readPdb("../data/2n77.pdb", "new_protein");  // Reading a .pdb file


    std::vector<std::pair< Point, unsigned> > P;  // Vector for the atoms of the protein

    std::map<unsigned int, bool> is_interface;

    for (int i = 0; i < pdb->n_atoms; i++) {
        P.push_back(std::make_pair(Point(pdb->atom[i].pt.x, pdb->atom[i].pt.y, pdb->atom[i].pt.z), i) ); // Filling the vector with the coordinates of the atoms and an index
        is_interface.insert(std::make_pair(i,false)); // Boolean : is the atom in interface cell
        //std::cout<<i<<std::endl;
    }

    //std::cout<<P.at(0)<<std::endl;

    Delaunay d_t(P.begin(), P.end());  // Creating Delaunay triangulation with the coordinates vector

    std::cout<<"validité : "<<d_t.is_valid()<<std::endl; //Checking the validity of the Delaunay triangulation

    Delaunay::Finite_cells_iterator cit;
    //Delaunay::Finite_vertices_iterator vit;

    bool is_relevant = false;

    std::vector<std::pair< Point, unsigned> > reduced_vector;

    for (cit = d_t.finite_cells_begin(); cit != d_t.finite_cells_end(); ++cit) { //iterating on the cells
        is_relevant = false;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                if (strcmp(pdb->atom[cit->vertex(j)->info()].chain, pdb->atom[cit->vertex(k)->info()].chain) != 0) { //checking if a cell contains 2 different chains
                    //std::cout<<strcmp(pdb->atom[cit->vertex(j)->info()].chain, pdb->atom[cit->vertex(k)->info()].chain)<<std::endl;
                    is_relevant = true; // if yes    to be kept
                    //std::cout<<pdb->atom[cit->vertex(j)->info()].chain<<pdb->atom[cit->vertex(k)->info()].chain<<std::endl;
                }
            }
        }
        if (is_relevant == true) {
            for (int j = 0; j < 4; j++) {
                is_interface[cit->vertex(j)->info()] = true;
            }
        }
        //std::cout<<"is_relevant : "<<is_identical<<std::endl;
    }

    for (int i = 0; i < pdb->n_atoms; i++) { // Filling the vector with the coordinates of the atoms
        if (is_interface[i] == true) {
            reduced_vector.push_back(std::make_pair(Point(pdb->atom[i].pt.x, pdb->atom[i].pt.y, pdb->atom[i].pt.z), i) );
        }
    }

    Delaunay interface_tr(reduced_vector.begin(), reduced_vector.end());  // Creating Delaunay triangulation with the coordinates vector

    // Check the coordinates and the number of atoms
    /*int cpt2 = 1;
    for (auto p : P) {
        std::cout<<cpt2<<" : "<<p<<std::endl;
        cpt2 ++;
    }*/

    /*std::ofstream oFileT("output.off",std::ios::out);
    // writing file output;
    oFileT << d_t;*/

    /// Barycenter Calculation

    double x_coord = 0.0;
    double y_coord = 0.0;
    double z_coord = 0.0;
    for (Delaunay::Finite_vertices_iterator vit = interface_tr.finite_vertices_begin(); vit != interface_tr.finite_vertices_end(); ++vit) {
        x_coord += vit->point().x() / interface_tr.number_of_vertices();
        y_coord += vit->point().y() / interface_tr.number_of_vertices();
        z_coord += vit->point().z() / interface_tr.number_of_vertices();
    }
    K::Point_3 barycenter(x_coord,y_coord,z_coord);
    std::cout<<"Barycenter : "<<x_coord<<"  "<<y_coord<<"  "<<z_coord<<std::endl;

    ///Calculate interface

    // Initialize iterators
    Delaunay::Finite_edges_iterator eit; // Edge iterator
    Delaunay::Cell_circulator circ; // Circulator for cells around a given edge
    Delaunay::Cell_circulator circ_copy; // Copy to check when all cells have been tested

    // Initialize variables
    std::map<K::Point_3, unsigned int > surf_points; // link between point and index for faces
    std::vector<std::vector<K::Point_3> > all_faces_indexes; // vector of faces (each line stores n points)
    std::vector<K::Point_3> face_indexes; //points of one face
    K::Point_3 p; // 3D point : in the triangulation
    K::Vector_3 dist; // 3D vector : distance to be normalized
    bool is_used = true; // To check if the point is not too far from the barycenter
    bool is_finite_tetrahedron = true; // To check if a tetrahedron contains the infinite point
    double size_lim = 34 * 34; // Max distance from the barycenter    ***** To be changed *****
    interface_tr = d_t; // If we want to work on the whole protein

    for (eit = interface_tr.finite_edges_begin(); eit != interface_tr.finite_edges_end(); ++eit) { // Iterate on the edges of the whole triangulation
        if((strcmp(pdb->atom[eit->first->vertex(eit->second)->info()].chain, pdb->atom[eit->first->vertex(eit->third)->info()].chain) != 0)){ // Check if the edge is at the interface
            //std::cout<<pdb->atom[eit->first->vertex(eit->second)->info()].chain<<"   "<<pdb->atom[eit->first->vertex(eit->third)->info()].chain<<std::endl;
            circ = interface_tr.incident_cells(*eit);
            circ_copy = circ;
            for (int i = 0; i < 4; i++) {
                //if (fabs(circ->vertex(i)->point().x()) < 0.001 || fabs(circ->vertex(i)->point().y()) < 0.001 || fabs(circ->vertex(i)->point().z()) < 0.001) {
                if(interface_tr.is_infinite(circ->vertex(i))) { //check if the vertex is the infinite one
                    is_finite_tetrahedron = false; //if yes   the tetrahedron will not be stored
                    //std::cout<<circ->vertex(i)->point()<<std::endl; //display infinite point
                }
            }
            do {
                p = interface_tr.dual(circ); // Find the dual of the cell (Voronoï Diagram)
                //face.push_back(p);
                //if (abs(p.x()) < size_lim && abs(p.y()) < size_lim && abs(p.z()) < size_lim && is_finite_tetrahedron == true) {
                if ((barycenter - p).squared_length() <= size_lim) { // Checking the distance between the barycenter and the point is in the limit
                    surf_points[p] = 0; // Store the point in the map with 0 value for the index
                }
                face_indexes.push_back(p); // Put the point in the vector of the current face
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
                    is_used = false; //we will not store the face
                    //std::cout<<all_faces_indexes[i][j].x()<<"  "<<all_faces_indexes[i][j].y()<<"  "<<all_faces_indexes[i][j].z()<<std::endl;
                }
            }
            if (is_used == true && is_finite_tetrahedron == true) { // Check is the face is in the limit and does not contain the infinte point
                all_faces_indexes.push_back(face_indexes); // Add that face to the intersurf vector
            }
            is_used = true;
            is_finite_tetrahedron = true; // Put the booleans back to true (at the begining the face is assumed to be valid)
            face_indexes.clear(); // Clear the current face vector for the next iteration
        }
    }

    ///Write the surface into a .off file -> interface.off
    // Open file to write .off
    std::ofstream fout;
    fout.open( "interface.off" );
    // .OFF Header : Number of points, Number of faces
    fout<<"OFF"<<std::endl<<surf_points.size()<<" "<<all_faces_indexes.size()<<" "<<0<<std::endl;

    //Write vertices into .off file
    int cpt_ind = 0; // Index following the wrinting order of the points
    std::map<K::Point_3, unsigned int>::iterator it_surf; //Initialize iterator to go through the surface map

    for(it_surf = surf_points.begin(); it_surf != surf_points.end(); ++it_surf) { // For all the points of the surface map
        fout<<it_surf->first.x()<<" "<<it_surf->first.y()<<" "<<it_surf->first.z()<<std::endl; // Write the current point int the .off file
        surf_points[it_surf->first] = cpt_ind; // Give the index the value of the counter
        cpt_ind++; // Increment the index
    }

    // Write indexes into .off file
    for (int i = 0; i < all_faces_indexes.size(); i++) { // Go through the faces vector
        fout<<all_faces_indexes[i].size()<<" "; // Write the number of points for the current face
        for (int j = 0; j < all_faces_indexes[i].size(); j++) { // Go through the current face
            fout<<surf_points[all_faces_indexes[i][j]]; // Write the index corresponding to the current point in the surface map
            if (j != all_faces_indexes[i].size() - 1) { // if the current point is NOT the last one
                fout<<" "; // Add a spacing
            }
        }
        fout<<std::endl; // Break
    }
    fout.close(); // Close the file


    /// Smooth the obtained surface with the CGAL Polyhedron face
    int degree = 4; // Degree of the subdivision
    Polyhedron poly_surf;
    std::ifstream fout_2( "interface.off" ); // Open file for reading
    CGAL::scan_OFF(fout_2,poly_surf,true); // Write the surface into the polyhedron structure + check errors
    if(!fout_2) { // Check if reading was successful
        std::cerr << "OFF reading failed.\n";
    }
    fout_2.close(); // Close file
    std::cout<<"polysurf empty : "<<poly_surf.is_empty()<<std::endl; // Check if Polyhedron is empty

    CGAL::Subdivision_method_3::CatmullClark_subdivision(poly_surf,degree); // Apply Catmull/Clark subdividion method to the polyhedron

    std::ofstream fout_3; // Open a new file for the smoothed interface
    fout_3.open( "smoothed_interface.off" );
    fout_3 << poly_surf;
    fout_3.close(); // Close file


    savePointsOFF("output2.off",interface_tr, pdb); // Write the protein into a .off file

    // Protein and interface data
    /*std::cout<<"nb vertices : "<<d_t.number_of_vertices()<<std::endl;
    std::cout<<"nb edges : "<<d_t.number_of_edges()<<std::endl;
    std::cout<<"nb faces : "<<d_t.number_of_facets()<<std::endl;
    std::cout<<"nb cells : "<<d_t.number_of_cells()<<std::endl;

    std::cout<<"nb vertices : "<<interface_tr.number_of_vertices()<<std::endl;
    std::cout<<"nb edges : "<<interface_tr.number_of_edges()<<std::endl;
    std::cout<<"nb faces : "<<interface_tr.number_of_facets()<<std::endl;
    std::cout<<"nb cells : "<<interface_tr.number_of_cells()<<std::endl;*/

    return 0;
}



static void savePointsOFF(const char* filename, Delaunay m_dt, PdbImage *pdb)
{
  // Open file to write .off
  std::ofstream fout;
  fout.open( filename );
  if( !fout ) {
      std::cout<<"Error"<<std::endl;
    return;
  }

  Delaunay::size_type n_vertices = m_dt.number_of_vertices();
  fout<<"[C]OFF"<<std::endl<<n_vertices<<" "<<m_dt.number_of_edges()<<" "<<m_dt.number_of_facets()<<std::endl<<std::endl;

  std::vector<Vertex_handle> TV(n_vertices + 1);
  Delaunay::size_type i = 0;

  // write points (get from point array)
  for(Delaunay::Finite_vertices_iterator vit=m_dt.finite_vertices_begin();
      vit!=m_dt.finite_vertices_end(); ++vit) {
    K::Point_3& p = vit->point();
    fout<<p.x()<<" "<<p.y()<<" "<<p.z()<<" ";
    if (strcmp( pdb->atom[vit->info()].chain, "A") == 0) {
        fout<<CGAL::RED;
    }
    else {
        fout<<CGAL::BLUE;
    }
    fout<<std::endl;
    TV[i++] = vit;
  }
  fout<<std::endl;

  CGAL::Unique_hash_map<Vertex_handle, std::size_t > V;

  V[m_dt.infinite_vertex()] = 0;
  for (i=1; i <= n_vertices; i++) {
    V[TV[i]] = i;
  }

  // write faces (get from point array)
  for(Delaunay::Finite_facets_iterator fit=m_dt.finite_facets_begin();
      fit!=m_dt.finite_facets_end(); ++fit) {
    fout<<3<<" ";
    for (int i = 0; i < 4; i++) {
        if (i != fit->second) {
            fout<<V[fit->first->vertex(i)]<<" ";
        }
    }
    fout<<std::endl;
  }
  fout.close();
}

