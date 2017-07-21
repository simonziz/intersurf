#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_writer_OFF.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/Writer_OFF.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

//#include <CGAL/triangulation_assertions.h>

//#include <CGAL/Complex_2_in_triangulation_3.h>
//#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
//#include <CGAL/make_surface_mesh.h>

#include <QGLViewer/qglviewer.h>


#include "pdbs.h"
#include "my_vertex_base.h"
#include "DelaunayMeshTriangulationGraphicsItem.h"

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
//typedef CGAL::Qt::TriangulationGraphicsItem<Delaunay> viewer;

typedef Delaunay::Point                          Point;
typedef Delaunay::Cell_handle                    Cell_handle;
typedef Delaunay::Facet                          Face;


//typedef CGAL::Complex_2_in_triangulation_3<Delaunay> C2t3;


static void savePointsOFF(const char* filename, Delaunay m_dt);

int main(void) {

    CGAL::Geomview_stream gv(CGAL::Bbox_3(-100, -100, -100, 60, 60, 60));
    gv.set_line_width(2);
    //gv.set_trace(true);
    gv.set_bg_color(CGAL::Color(200, 200, 200));
    // gv.clear();


    //PdbImage *pdb = hex_readPdb("../data/toto.pdb", "new_protein");  // Reading a .pdb file
    PdbImage *pdb = hex_readPdb("../data/2n77_reduced.pdb", "new_protein");  // Reading a .pdb file


    std::vector<std::pair< Point, unsigned> > P;  // Vector for the atoms of the protein

    std::map<unsigned int, bool> is_interface;

    for (int i = 0; i < pdb->n_atoms; i++) { // Filling the vector with the coordinates of the atoms
        P.push_back(std::make_pair(Point(pdb->atom[i].pt.x, pdb->atom[i].pt.y, pdb->atom[i].pt.z), i) );
        is_interface.insert(std::make_pair(i,false));
        //std::cout<<i<<std::endl;
    }

    //std::cout<<P.at(0)<<std::endl;

    Delaunay d_t(P.begin(), P.end());  // Creating Delaunay triangulation with the coordinates vector

    std::cout<<"validité : "<<d_t.is_valid()<<std::endl;

    Delaunay::Finite_cells_iterator cit;
    //Delaunay::Finite_vertices_iterator vit;

    bool is_relevant = false;

    std::vector<std::pair< Point, unsigned> > reduced_vector;

    for (cit = d_t.finite_cells_begin(); cit != d_t.finite_cells_end(); ++cit) {
        is_relevant = false;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                if (strcmp(pdb->atom[cit->vertex(j)->info()].chain, pdb->atom[cit->vertex(k)->info()].chain) != 0) {
                    //std::cout<<strcmp(pdb->atom[cit->vertex(j)->info()].chain, pdb->atom[cit->vertex(k)->info()].chain)<<std::endl;
                    is_relevant = true;
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
    int cpt2 = 1;
    /*for (auto p : P) {
        std::cout<<cpt2<<" : "<<p<<std::endl;
        cpt2 ++;
    }*/

    std::ofstream oFileT("output.off",std::ios::out);
    // writing file output;
    oFileT << d_t;

    savePointsOFF("output2.off",interface_tr);
    //CGAL::output_surface_facets_to_off(oFileT,c2t3);
    std::cout<<"nb vertices : "<<d_t.number_of_vertices()<<std::endl;
    std::cout<<"nb edges : "<<d_t.number_of_edges()<<std::endl;
    std::cout<<"nb faces : "<<d_t.number_of_facets()<<std::endl;
    std::cout<<"nb cells : "<<d_t.number_of_cells()<<std::endl;

    std::cout<<"nb vertices : "<<interface_tr.number_of_vertices()<<std::endl;
    std::cout<<"nb edges : "<<interface_tr.number_of_edges()<<std::endl;
    std::cout<<"nb faces : "<<interface_tr.number_of_facets()<<std::endl;
    std::cout<<"nb cells : "<<interface_tr.number_of_cells()<<std::endl;

    std::cout << "Drawing 3D Delaunay triangulation in wired mode.\n";
    gv.set_wired(true);
    gv << d_t;
    sleep(100);
    gv.clear();

    return 0;
}



static void savePointsOFF(const char* filename, Delaunay m_dt)
{
  std::ofstream fout;
  fout.open( filename );
  if( !fout ) {
    //showError( QObject::tr("Error: cannot open file %1 for writting.").arg(filename) );
      std::cout<<"Error"<<std::endl;
    return;
  }

  std::ostream *pOut = &fout;

  /* Use CGAL::File_writer_OFF to write points */
  // initialize header_OFF
  CGAL::File_header_OFF header(false,  // true: binary output; false: ASCII
                               false,  // true: no comments in file
                               false,  // true: Geomview SKEL format
                               true);  // true: verbosity on; false: verbosity off
  // a simpler way to initialize header_OFF
//  CGAL::File_header_OFF header(true);  // true: verbosity on
//                                       // (ASCII output, comments, no SKEL)
  CGAL::File_writer_OFF writer( header );

  // write header
  /*writer.write_header(*pOut,  // output ostream
                      m_dt.number_of_vertices,  // number of points/vertices
                      0,  // number of halfedges
                      0,  // number of facets
                      false);  // true: has normals*/

  writer.write_header(*pOut,m_dt.number_of_vertices(),m_dt.number_of_edges(),m_dt.number_of_facets(),false);

  // write points (get from point array)
  for(Delaunay::Finite_vertices_iterator vit=m_dt.finite_vertices_begin();
      vit!=m_dt.finite_vertices_end(); ++vit) {
    K::Point_3& p = vit->point();
    writer.write_vertex( p.x(), p.y(), p.z() );
  }

  writer.write_facet_header();
  // write faces (get from point array)
  for(Delaunay::Finite_facets_iterator fit=m_dt.finite_facets_begin();
      fit!=m_dt.finite_facets_end(); ++fit) {
    writer.write_facet_begin(3);
    for (int i = 0; i < 4; i++) {
        if (i != fit->second) {
            writer.write_facet_vertex_index(fit->first->vertex(i)->info());
        }
    }
    writer.write_facet_end();
    //std::cout<<"face : "<<fit->first->vertex(1)->info()<<std::endl;
  }
  /*for(Delaunay::Finite_cells_iterator cit=m_dt.finite_cells_begin();
      cit!=m_dt.finite_cells_end(); ++cit) {
    for (int i = 0; i < 4; i++) {
        writer.write_facet_begin(3);
        writer.write_facet_vertex_index(cit->vertex(i % 4)->info());
        writer.write_facet_vertex_index(cit->vertex((i + 1) % 4)->info());
        writer.write_facet_vertex_index(cit->vertex((i + 2) % 4)->info());
        writer.write_facet_end();
    }
  }*/
  // write footer
  writer.write_footer();
  fout.close();
}

