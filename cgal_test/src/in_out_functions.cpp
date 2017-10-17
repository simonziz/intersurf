#include "in_out_functions.h"

/// Function to write a CGAL::Delaunay_3 structure into a .off file
/// with the image of the .pdb to identify the chains

void savePointsOFF(const char* filename, Delaunay m_dt, PdbImage *pdb)
{
  // Open file to write .off
  std::ofstream fout;
  fout.open( filename );
  if( !fout ) { // Print message if failure to open
    std::cout<<"Error"<<std::endl;
    //return; // To kill the program if needed
  }

  Delaunay::size_type n_vertices = m_dt.number_of_vertices(); // Returns the number of vertices of the triangulation

  // Write the .off header
  fout<<"[C]OFF"<<std::endl<<n_vertices<<" "<<m_dt.number_of_edges()<<" "<<m_dt.number_of_facets()<<std::endl<<std::endl;

  std::vector<Vertex_handle> TV(n_vertices + 1); // Vector to store the indexes of the points
  Delaunay::size_type i = 0;

  // write points (get from point array)
  // Iterate on the vertices of the triangulation
  for(Delaunay::Finite_vertices_iterator vit=m_dt.finite_vertices_begin(); vit!=m_dt.finite_vertices_end(); ++vit) {
      K::Point_3& p = vit->point(); // Current vertex
      fout<<p.x()<<" "<<p.y()<<" "<<p.z()<<" "; // Write its coordinates
      // Choose color regarding the chain
      if (strcmp( pdb->atom[vit->info()].chain, "A") == 0) {
        fout<<CGAL::BLUE;
      }
      else {
        fout<<CGAL::RED;
      }
      fout<<std::endl; // Break
      TV[i++] = vit; // Store the vertex handle (pointer to the current vertex) at its index
  }
  fout<<std::endl; // Break

  // Iterate on the vertices of the triangulation
  for(Delaunay::Finite_vertices_iterator vit=m_dt.finite_vertices_begin(); vit!=m_dt.finite_vertices_end(); ++vit) {
      vit->info();
  }

  CGAL::Unique_hash_map<Vertex_handle, std::size_t > V; // Map with unique indexes
  V[m_dt.infinite_vertex()] = 0; // Same size as the triangulation

  for (i=1; i <= n_vertices; i++) { // Go through the vertices
    V[TV[i]] = i; // And give each point its corresponding index
  }

  // write faces (get from point array)
  // Iterate on the faces of the triangulation
  for(Delaunay::Finite_facets_iterator fit=m_dt.finite_facets_begin(); fit!=m_dt.finite_facets_end(); ++fit) {
    fout<<3<<" "; // 3 because there are only triangles + Spacing
    for (int i = 0; i < 4; i++) { // Go through the 4 points of the cell
        if (i != fit->second) { // Check which vertex is opposed to the current face
            fout<<V[fit->first->vertex(i)]<<" "; // Write the remaining veritces into the .off
        }
    }
    fout<<std::endl; // Break
  }
  fout.close(); // Close file
}

/// Function to write the surface into a .off file

void saveSurface0FF(const char* filename, std::map<K::Point_3, unsigned int > surf_points,
                    std::vector<std::vector<K::Point_3> > all_faces_indexes)
{
    // Open file to write .off
    std::ofstream fout;
    fout.open( filename );
    // .OFF Header : Number of points, Number of faces
    fout<<"OFF"<<std::endl<<surf_points.size()<<" "<<all_faces_indexes.size()<<" "<<0<<std::endl;

    //Write vertices into .off file
    int cpt_ind = 0; // Index following the wrinting order of the points

    // For all the points of the surface map
    for(std::map<K::Point_3, unsigned int>::iterator it_surf = surf_points.begin(); it_surf != surf_points.end(); ++it_surf) {
        // Write the current point int the .off file
        fout<<it_surf->first.x()<<" "<<it_surf->first.y()<<" "<<it_surf->first.z()<<std::endl;
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
}

/// Function to write a .pdb file into a CGAL::Delaunay_3 structure

Delaunay readPdb(PdbImage* pdb, std::map<unsigned int, bool> is_interface)
{
    std::vector<std::pair< Point, unsigned> > atoms_coordinates;  // Vector for the atoms of the protein    

    for (int i = 0; i < pdb->n_atoms; i++) {
        // Filling the vector with the coordinates of the atoms and an index
        atoms_coordinates.push_back(std::make_pair(Point(pdb->atom[i].pt.x, pdb->atom[i].pt.y, pdb->atom[i].pt.z), i) );
        is_interface.insert(std::make_pair(i,false)); // Boolean : is the atom in interface cell
        //std::cout<<i<<std::endl;
    }

    //std::cout<<atoms_coordinates.at(0)<<std::endl;

    // Creating Delaunay triangulation with the coordinates vector
    Delaunay d_t(atoms_coordinates.begin(), atoms_coordinates.end());
    return d_t;
}


/// Function to write the interface from a .off to a CGAL::Polyhedron

Polyhedron makePolyhedron(const char* filename, bool show_warnings)
{
    Polyhedron poly_surf;
    std::ifstream fout_2( filename ); // Open file for reading
    CGAL::scan_OFF(fout_2,poly_surf,show_warnings); // Write the surface into the polyhedron structure + check errors
    if(!fout_2) { // Check if reading was successful
        std::cerr << "OFF reading failed.\n";
    }
    fout_2.close(); // Close file

    return poly_surf;
}
