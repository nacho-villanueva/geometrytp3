//=============================================================================
//
//  Implementation of class Parametrizer
//
//=============================================================================

#include <iostream>
#include <math.h>
#include <assert.h>

#include "Parametrizer.h"

Parametrizer::Parametrizer(MyMesh &iMesh):
m_mesh(iMesh)
{
  OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_u");
  OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_h");

  OpenMesh::VProp<MyMesh::Scalar>(m_mesh, "theta");

  std::vector<OpenMesh::VertexHandle> boundaries;
  auto total_length = this->list_boundary(boundaries);
  this->set_boundary(boundaries, total_length);
}

void Parametrizer::solve_linear_system(gmm::dense_matrix<float> & iA, std::vector<float>& iB, std::vector<float> & oX)
{
  // for better performances, we should switch to sparse matrices and use solvers for
  // sparse systems but this requires the use of external third-party libraries...
  oX.clear();
  oX.resize(iB.size());
  gmm::lu_solve(iA, oX, iB);
}

void Parametrizer::set_boundary(std::vector<OpenMesh::VertexHandle> boundaries, float total_length) {
  auto vertex_theta = OpenMesh::VProp<MyMesh::Scalar>(m_mesh, "theta");

  float total_theta = 0;
  float tl2 = 0;
  auto it = boundaries.begin();
  auto lastVertex = *it;
  vertex_theta[*it] = 0;
  it++;

  for(;it != boundaries.end(); it++) {
    auto length = (m_mesh.point(lastVertex) - m_mesh.point(*it)).length();
    auto theta = length * ( (2 * M_PI) / total_length);
    tl2 += length;
    std::cout << "(" << lastVertex.idx() << ", " << it->idx() << ") Length: " << length << " - Theta: " << theta << " - Total Theta: " << total_theta << std::endl;
    total_theta += theta; 
    vertex_theta[*it] = total_theta;
    lastVertex = *it;
  }

  auto last_length = (m_mesh.point(boundaries.front()) - m_mesh.point(boundaries.back())).length();
  auto theta = last_length * ( (2 * M_PI) / total_length);
  std::cout << "(" << boundaries.back().idx() << ", " << boundaries.front().idx() << ") Length: " << last_length << " Total Theta: " << total_theta << std::endl;
}

/*float Parametrizer::list_boundary(std::vector<OpenMesh::VertexHandle> & oBoundary) {
  int nb_v = m_mesh.n_vertices();
  bool found = true;
  float total_length = 0;

  MyMesh::VertexIter v_iter;
  // Find first boundary vertex
  OpenMesh::VertexHandle currentVertex;
  for(v_iter = m_mesh.vertices().begin(); v_iter != m_mesh.vertices_end(); ++v_iter) {
    if(m_mesh.is_boundary(*v_iter)) {
      oBoundary.push_back(*v_iter);
      currentVertex = *v_iter;
      std::cout << "Found first boundary vertex: "<< v_iter -> idx() <<"\n";
      break;
    }
  }


  bool finished = false;
  MyMesh::VertexVertexIter vv_iter;
  while(!finished) {
    finished = true;
    for(vv_iter = m_mesh.vv_iter(currentVertex); vv_iter.is_valid(); ++vv_iter) {
      if(m_mesh.is_boundary(*vv_iter)) {
        std::cout << "Adding vertex to boundary. Id:" << (*vv_iter).idx() << "\n";
        if(std::find(oBoundary.begin(), oBoundary.end(), *vv_iter) == oBoundary.end()) {
          currentVertex = *vv_iter;
          total_length += (m_mesh.point(currentVertex) - m_mesh.point(oBoundary.back())).length();
          oBoundary.push_back(currentVertex);
          finished = false;
          break;
        }
      }
    }
  }

  std::cout << "Total boundary length: " << total_length << std::endl; 
  return total_length;
}*/

float Parametrizer::list_boundary(std::vector<OpenMesh::VertexHandle> & oBoundary) {
  int nb_v = m_mesh.n_vertices();
  bool found = true;
  float total_length = 0;
  
  MyMesh::VertexIter v_it;
  MyMesh::VertexVertexIter vv_it;

  for(v_it = m_mesh.vertices().begin(); v_it != m_mesh.vertices_end(); ++v_it) {
    if(m_mesh.is_boundary(*v_it)){
      oBoundary.push_back(*v_it);
      break;
    }
  }

  while (found){
    found = false;
    for (vv_it = m_mesh.vv_iter(oBoundary.back()); vv_it.is_valid(); ++v_it){
      if ((m_mesh.is_boundary(*vv_it)) && (std::find(oBoundary.begin(), oBoundary.end(), *vv_it) == oBoundary.end())){
        float length = (m_mesh.point(*vv_it) - m_mesh.point(oBoundary.back())).length();
        total_length += length;
        oBoundary.push_back(*vv_it);
        found = true;
        break;
      }
    }
  }

  return total_length;
}





void Parametrizer::calc_uniform_parameterization()
{
  // ------------- IMPLEMENT HERE ---------
  // TASK 3.1 Uniform map computation:
  // Search and order boundary vertices
  // Compute boundary parameters
  // Solve the linear system for internal vertex parameters using solve_linear_system()
  // Store the parameters in the vparam_u vertex property
  // ------------- IMPLEMENT HERE ---------

  auto vparam_u = OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_u");
  auto index = OpenMesh::VProp<MyMesh::Scalar>(m_mesh, "index");
  auto boundary_theta = OpenMesh::VProp<MyMesh::Scalar>(m_mesh, "theta");
  
  gmm::dense_matrix<float> W(m_mesh.n_vertices(), m_mesh.n_vertices());
  gmm::clear(W);

  std::vector<float> bx(m_mesh.n_vertices(), 0);
  std::vector<float> by(m_mesh.n_vertices(), 0);

  for (const auto& vi : m_mesh.vertices()) {
    if (m_mesh.is_boundary(vi)) {
      assert(0 <= (size_t)vi.idx() && (size_t)vi.idx() < m_mesh.n_vertices());
      W(vi.idx(), vi.idx()) = 1.0f;

      bx[vi.idx()] = cos(boundary_theta[vi]);
      by[vi.idx()] = sin(boundary_theta[vi]);
    }
    else {
      std::cout << "Valence of " << vi.idx() << ": " << -m_mesh.valence(vi) << std::endl;
      int valence = 0;
      for (const auto& eh : vi.edges()) valence++;
      W(vi.idx(), vi.idx()) = -valence;
      // -m_mesh.valence(vi)
    }
  }

  for (const auto& ei : m_mesh.edges()) {
    auto v0 = ei.v0();
    auto v1 = ei.v1();
    W(v0.idx(), v1.idx()) = 1.0f;
    W(v1.idx(), v0.idx()) = 1.0f;
  }

  //for(int j = 0; j < m_mesh.n_vertices(); j++) {
  //  for(int i = 0; i < m_mesh.n_vertices(); i++){ 
  //    std::cout << W(i, j) << " ";
  //  }
  //  std::cout << "\n";
  //}

  for(int j = 0; j < m_mesh.n_vertices(); j++) {
    //std::cout << "(" << bx[j] << ", " << by[j] << ")\n"; 
  }

  std::vector<float> x(m_mesh.n_vertices());
  std::vector<float> y(m_mesh.n_vertices());

  this->solve_linear_system(W, bx, x);
  this->solve_linear_system(W, by, y);

  for (const auto& vi : m_mesh.vertices()) {
    std::cout << "vi: " << vi.idx() << " (" << x[vi.idx()] << ", " << y[vi.idx()] << ")\n"; 
    vparam_u[vi][0] = x[vi.idx()];
    vparam_u[vi][1] = y[vi.idx()];
  }
}

void Parametrizer::calc_harmonic_parameterization()
{
  // ------------- IMPLEMENT HERE ---------
  // TASK 3.2 harmonic map computation:
  // Search and order boundary vertices
  // Compute boundary parameters
  // Solve the linear system for internal vertex parameters using solve_linear_system()
  // Store the parameters in the vparam_h vertex property
  // ------------- IMPLEMENT HERE ---------
  
}

void Parametrizer::computeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iRepeats, ParameterizationMode imode)
{
  // ------------- IMPLEMENT HERE ---------
  // TASK 3.3 Compute texture coordinates for textured mesh
  // rendering using a texture image of dimension iTextureWidth*iTextureHeights and iRepeats repeats
  // and store them in the standard texCoord property (hint: use the set_texcoord2D function of OpenMesh).
  // If imode is equals to Uniform, compute the texture coordinates using the
  // parameters stored in vparam_u.
  // If imode is equals to Harmonic, compute the texture coordinates using the
  // parameters stored in vparam_h.
  // ------------- IMPLEMENT HERE ---------
  
}

void Parametrizer::calc_distortion(ParameterizationMode imode)
{
  float angle_distortion=0., area_distortion=0.;
  
  // ------------- IMPLEMENT HERE ---------
  // TASK 3.4 Compute distortion of triangle areas and angles
  // and print it in the output window.
  // If imode is equals to Uniform, uniform map distortion has to be
  // computed.
  // If imode is equals to Harmonic, harmonic map distortion has to be
  // computed.
  // ------------- IMPLEMENT HERE ---------
  
  std::cout << "Parameterization distortion: " << std::endl;
  std::cout << (imode==Uniform? " * Uniform map: ": " * Harmonic map: ") << std::endl;
  std:: cout << "      Angle distortion: " << angle_distortion << " Area distortion: " << area_distortion << std::endl;
}

void Parametrizer::flattenMesh(ParameterizationMode imode)
{
  OpenMesh::VProp<OpenMesh::Vec2f> vparam_u = OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_u");
  OpenMesh::VProp<OpenMesh::Vec2f> vparam_h = OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_h");
  
  MyMesh::VertexIter v_it, v_end(m_mesh.vertices_end());
  for (v_it=m_mesh.vertices_begin(); v_it != v_end; ++v_it)
  {
    const OpenMesh::Vec2f UV = (imode==Uniform? vparam_u[v_it]: vparam_h[v_it]);
    m_mesh.set_point(v_it, OpenMesh::Vec3f(UV[0], UV[1], 0));
  }
  m_mesh.update_normals();
}

//=============================================================================
