#include "includes.h"
#include "typedefs.cpp"

#include "in_out_functions.h"
#include "calculations.h"


int main(void) {

    /// Read a .pdb file into an array
    PdbImage *pdb = hex_readPdb("../data/2n77.pdb", "new_protein");  // Reading a .pdb file

    std::map<unsigned int, bool> is_interface;
    Delaunay d_t = readPdb(pdb, is_interface);

    std::cout<<"validité : "<<d_t.is_valid()<<std::endl; // Checking the validity of the Delaunay triangulation

    Delaunay interface_tr = makeInterfaceTriangulation(pdb, is_interface, d_t);

    // Check the coordinates and the number of atoms
    /*int cpt2 = 1;
    for (auto p : P) {
        std::cout<<cpt2<<" : "<<p<<std::endl;
        cpt2 ++;
    }*/

    /// Barycenter + Max Distance Calculation
    // Barycenter
    K::Point_3 barycenter = calculateBarycenter(interface_tr);
    std::cout<<"Barycenter : "<<barycenter.x()<<"  "<<barycenter.y()<<"  "<<barycenter.z()<<std::endl;

    // Max distance
    double squared_max_dist = calculateMaxDist(interface_tr, barycenter);
    std::cout<<"max dist : "<<squared_max_dist<<std::endl;


    ///Calculate interface
    // Initialize variables
    std::map<K::Point_3, unsigned int > surf_points; // link between point and index for faces
    std::vector<std::vector<K::Point_3> > all_faces_indexes; // vector of faces (each line stores n points)

    //interface_tr = d_t; // If we want to work on the whole protein

    calculateInterface(pdb, d_t, barycenter, squared_max_dist, &surf_points, &all_faces_indexes);

    //std::cout<<"points end : "<<surf_points.empty()<<std::endl;
    //std::cout<<"faces end : "<<all_faces_indexes.empty()<<std::endl;

    ///Write the surface into a .off file -> interface.off
    saveSurface0FF("interface.off", surf_points, all_faces_indexes);


    /// Smooth the obtained surface with the CGAL Polyhedron face
    Polyhedron poly_surf = makePolyhedron("interface.off", false); // Make polyhedron from .off file
    std::cout<<"polysurf empty : "<<poly_surf.is_empty()<<std::endl; // Check if Polyhedron is empty

    // Apply Catmull/Clark subdivision method to the polyhedron
    int degree = 2; // Degree of the subdivision
    CGAL::Subdivision_method_3::CatmullClark_subdivision(poly_surf,degree);

    std::ofstream fout_3; // Open a new file for the smoothed interface
    fout_3.open( "smoothed_interface.off" );
    fout_3 << poly_surf;
    fout_3.close(); // Close file


    /// Remove useless triangles
    K::Point_3 p_surf;
    K::Point_3 p_triang;
    double min_dist = 20000;
    double current_dist = 0;
    for (Polyhedron::Vertex_iterator vit = poly_surf.vertices_begin(); vit != poly_surf.vertices_end(); ++vit) {
        //std::cout<<fit->halfedge()->vertex()->point()<<std::endl;
        p_surf = vit->point();
        for (Delaunay::Finite_vertices_iterator vit_tr = d_t.finite_vertices_begin();
             vit_tr != d_t.finite_vertices_end(); ++vit_tr) {
            p_triang = vit_tr->point();
            current_dist = sqrt((p_triang - p_surf).squared_length());
            min_dist = std::min(current_dist, min_dist);
        }
        if(min_dist > 12) {
            std::cout<<"min distance : "<<min_dist<<std::endl;
        }
        //if( ( p - barycenter).squared_length() >= squared_max_dist * 0.6) {
        if( min_dist > 12 ) {
            if(vit->halfedge()->is_border() == false)
            {
                poly_surf.erase_facet(vit->halfedge());
            }
        }
        min_dist = 20000;
    }
    CGAL::Subdivision_method_3::CatmullClark_subdivision(poly_surf);

    std::ofstream fout_4; // Open a new file for the reduced smoothed interface
    fout_4.open( "reduced_interface.off" );
    fout_4 << poly_surf;
    fout_4.close(); // Close file


    ///Write the protein into a .off file -> complex.off
    savePointsOFF("complex.off",interface_tr, pdb); // Write the protein into a .off file

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

