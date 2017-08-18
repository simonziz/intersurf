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

#include <CGAL/Unique_hash_map.h>

#include <QGLViewer/qglviewer.h>


#include "pdbs.h"
#include "my_vertex_base.h"
#include "DelaunayMeshTriangulationGraphicsItem.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

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
typedef Delaunay::Vertex_handle                  Vertex_handle;
typedef Delaunay::Facet                          Face;


typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;


//typedef CGAL::Complex_2_in_triangulation_3<Delaunay> C2t3;


static void savePointsOFF(const char* filename, Delaunay m_dt, PdbImage *pdb);

int main(void) {

    /*CGAL::Geomview_stream gv(CGAL::Bbox_3(-100, -100, -100, 60, 60, 60));
    gv.set_line_width(2);
    //gv.set_trace(true);
    gv.set_bg_color(CGAL::Color(200, 200, 200));
    // gv.clear();*/


    //PdbImage *pdb = hex_readPdb("../data/toto.pdb", "new_protein");  // Reading a .pdb file
    PdbImage *pdb = hex_readPdb("../data/2n77.pdb", "new_protein");  // Reading a .pdb file


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
    /*int cpt2 = 1;
    for (auto p : P) {
        std::cout<<cpt2<<" : "<<p<<std::endl;
        cpt2 ++;
    }*/

    /*std::ofstream oFileT("output.off",std::ios::out);
    // writing file output;
    oFileT << d_t;*/

    Delaunay::Finite_edges_iterator eit;
    Delaunay::Cell_circulator circ;
    Delaunay::Cell_circulator circ_copy;

    //std::vector<std::pair< unsigned int, std::vector<K::Point_3> > > intersurf;
    //std::vector<K::Point_3> face;
    //unsigned int nb_points_face;
    unsigned int index_of_point = 0;

    std::map<K::Point_3, unsigned int > surf_points;
    std::vector<std::vector<K::Point_3> > all_faces_indexes;
    std::vector<K::Point_3> face_indexes;
    int nb_points = -1;
    K::Point_3 p;
    K::Point_3 check_dist;
    K::Vector_3 dist;
    bool is_used = true;
    bool is_finite_tetrahedron = true;
    unsigned int size_lim = 50;

    for (eit = interface_tr.finite_edges_begin(); eit != interface_tr.finite_edges_end(); ++eit) {

        if((strcmp(pdb->atom[eit->first->vertex(eit->second)->info()].chain, pdb->atom[eit->first->vertex(eit->third)->info()].chain) != 0)){
            //std::cout<<pdb->atom[eit->first->vertex(eit->second)->info()].chain<<"   "<<pdb->atom[eit->first->vertex(eit->third)->info()].chain<<std::endl;
            circ = interface_tr.incident_cells(*eit);
            circ_copy = circ;
            for (int i = 0; i < 4; i++) {
                if (fabs(circ->vertex(i)->point().x()) < 0.001 || fabs(circ->vertex(i)->point().y()) < 0.001 || fabs(circ->vertex(i)->point().z()) < 0.001) {
                    is_finite_tetrahedron = false;
                    std::cout<<circ->vertex(i)->point()<<std::endl;
                }

            }

            //nb_points_face = 0;
            do {
                p = interface_tr.dual(circ);
                //fout<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
                //face.push_back(p);
                if (abs(p.x()) < size_lim && abs(p.y()) < size_lim && abs(p.z()) < size_lim && is_finite_tetrahedron == true) {
                    surf_points[p] = 0;//surf_points[p];
                }
                /*if (surf_points.size() != nb_points) {
                    surf_points[p] = surf_points.size() - 1;
                    //std::cout<<surf_points.size()<<std::endl;
                }*/
                face_indexes.push_back(p);
                nb_points = surf_points.size();
                circ ++;
                //nb_points_face ++;
            } while(circ != circ_copy);
            //intersurf.push_back(std::make_pair(nb_points_face, face));
            //face.clear();
            check_dist = face_indexes[0];
            for (int j = 0; j < face_indexes.size(); j++) {
                dist = check_dist - face_indexes[j];
                if ( dist.squared_length() > 600) {//abs(face_indexes[j].x()) >= size_lim || abs(face_indexes[j].y()) >= size_lim || abs(face_indexes[j].z()) > size_lim) {
                    is_used = false;
                    //std::cout<<all_faces_indexes[i][j].x()<<"  "<<all_faces_indexes[i][j].y()<<"  "<<all_faces_indexes[i][j].z()<<std::endl;
                }
                check_dist = face_indexes[j];
            }
            if (is_used == true && is_finite_tetrahedron == true) {
                all_faces_indexes.push_back(face_indexes);
            }
            is_used = true;
            is_finite_tetrahedron = true;
            face_indexes.clear();
        }
    }


    /*unsigned int total_vertices = 0;
    for (int i = 0; i < intersurf.size(); i++) {
        total_vertices += intersurf[i].first;
    }*/

    std::ofstream fout;
    fout.open( "interface.off" );
    //fout<<"OFF"<<std::endl<<total_vertices<<" "<<intersurf.size()<<" "<<0<<std::endl<<std::endl;
    fout<<"OFF"<<std::endl<<surf_points.size()<<" "<<all_faces_indexes.size()<<" "<<0<<std::endl<<std::endl;
    //Write vertices into .off file
    /*for (int i = 0; i < intersurf.size(); i++) {
        for (int j = 0; j < intersurf[i].second.size(); j++) {
            fout<<intersurf[i].second[j]<<std::endl;
        }
    }*/
    int cpt_ind = 0;
    std::map<K::Point_3, unsigned int>::iterator it_surf;

    for(it_surf = surf_points.begin(); it_surf != surf_points.end(); ++it_surf) {
        //if (abs(it_surf->first.x()) < size_lim && abs(it_surf->first.y()) < size_lim && abs(it_surf->first.z()) < size_lim) {
            fout<<it_surf->first.x()<<" "<<it_surf->first.y()<<" "<<it_surf->first.z()<<" "<<std::endl;
            surf_points[it_surf->first] = cpt_ind;
            cpt_ind++;
            //std::cout<<it_surf->second<<std::endl;
        //}

    }
    fout<<std::endl;
    /*for (int i = 0; i < surf_points.size(); i++) {
        std::cout<<surf_points.at(i)<<std::endl;
    }*/

    //Write indexes into .off file

    /*int cpt_index = 0;
    for (int i = 0; i < intersurf.size(); i++) {
        fout<<intersurf[i].first<<" ";
        for (int j = 0; j < intersurf[i].second.size(); j++) {
            fout<<cpt_index + j<<" ";
        }
        cpt_index += intersurf[i].first;
        fout<<std::endl;
    }*/

    for (int i = 0; i < all_faces_indexes.size(); i++) {
        //std::cout<<all_faces_indexes[i][1].x()<<std::endl;
        /*for (int j = 0; j < all_faces_indexes[i].size(); j++) {
            if (abs(all_faces_indexes[i][j].x()) >= size_lim || abs(all_faces_indexes[i][j].y()) >= size_lim || abs(all_faces_indexes[i][j].z()) > size_lim) {
                is_used = false;
                //std::cout<<all_faces_indexes[i][j].x()<<"  "<<all_faces_indexes[i][j].y()<<"  "<<all_faces_indexes[i][j].z()<<std::endl;
            }
        }*/
        //if (is_used == true) {
            fout<<all_faces_indexes[i].size()<<" ";
            for (int j = 0; j < all_faces_indexes[i].size(); j++) {
                fout<<surf_points[all_faces_indexes[i][j]]<<" ";
            }
            fout<<std::endl;
        //}
        //is_used = true;

    }

    fout.close();

    savePointsOFF("output2.off",interface_tr, pdb);
    //CGAL::output_surface_facets_to_off(oFileT,c2t3);
    std::cout<<"nb vertices : "<<d_t.number_of_vertices()<<std::endl;
    std::cout<<"nb edges : "<<d_t.number_of_edges()<<std::endl;
    std::cout<<"nb faces : "<<d_t.number_of_facets()<<std::endl;
    std::cout<<"nb cells : "<<d_t.number_of_cells()<<std::endl;

    std::cout<<"nb vertices : "<<interface_tr.number_of_vertices()<<std::endl;
    std::cout<<"nb edges : "<<interface_tr.number_of_edges()<<std::endl;
    std::cout<<"nb faces : "<<interface_tr.number_of_facets()<<std::endl;
    std::cout<<"nb cells : "<<interface_tr.number_of_cells()<<std::endl;

    /*std::cout << "Drawing 3D Delaunay triangulation in wired mode.\n";
    gv.set_wired(true);
    gv << d_t;
    sleep(5);
    gv.clear();*/

    return 0;
}



static void savePointsOFF(const char* filename, Delaunay m_dt, PdbImage *pdb)
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
                      0,  // number of facets
                      0,  // number of edges
                      false);  // true: has normals*/

  Delaunay::size_type n_vertices = m_dt.number_of_vertices();

  //writer.write_header(*pOut,n_vertices,m_dt.number_of_facets(),m_dt.number_of_edges(),false);
  //writer.write_header(*pOut,n_vertices,Delaunay::size_type(128),m_dt.number_of_edges(),false);
  fout<<"[C]OFF"<<std::endl<<n_vertices<<" "<<m_dt.number_of_edges()<<" "<<m_dt.number_of_facets()<<std::endl<<std::endl;

  std::vector<Vertex_handle> TV(n_vertices + 1);
  Delaunay::size_type i = 0;

  // write points (get from point array)
  for(Delaunay::Finite_vertices_iterator vit=m_dt.finite_vertices_begin();
      vit!=m_dt.finite_vertices_end(); ++vit) {
    K::Point_3& p = vit->point();
    //writer.write_vertex( p.x(), p.y(), p.z() );
    fout<<p.x()<<" "<<p.y()<<" "<<p.z()<<" ";
    //fout<<"1.000"<<" "<<"0.000"<<" "<<"0.000"<<" "<<0.75;
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

  //writer.write_facet_header();
  // write faces (get from point array)
  for(Delaunay::Finite_facets_iterator fit=m_dt.finite_facets_begin();
      fit!=m_dt.finite_facets_end(); ++fit) {
    //writer.write_facet_begin(3);
    fout<<3<<" ";
    for (int i = 0; i < 4; i++) {
        if (i != fit->second) {
            //writer.write_facet_vertex_index(V[fit->first->vertex(i)]);
            fout<<V[fit->first->vertex(i)]<<" ";
        }
    }
    //fout<<CGAL::RED;
    //writer.write_facet_end();
    fout<<std::endl;
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
  //writer.write_footer();
  fout.close();
}

