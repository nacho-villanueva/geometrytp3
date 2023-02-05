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

  OpenMesh::VProp<OpenMesh::Scalar>(m_mesh, "theta");
  this->set_boundary();
}

void Parametrizer::solve_linear_system(gmm::dense_matrix<float> & iA, std::vector<float>& iB, std::vector<float> & oX)
{
  // for better performances, we should switch to sparse matrices and use solvers for
  // sparse systems but this requires the use of external third-party libraries...
  oX.clear();
  oX.resize(iB.size());
  gmm::lu_solve(iA, oX, iB);
}

// TODO: Add to header
void Parametrizer::set_boundary() {
  OpenMesh::Scalar total_length = 0;
  for (const auto& e : m_mesh.edges()) { 
    if (e.is_boundary()) {
      total_length += m_mesh.calc_edge_sqr_length(e);
    }
  }

  OpenMesh::SmartVertexHandle start;
  for (const auto& he : m_mesh.halfedges()) {
    if (m_mesh.is_boundary(he)) {
      start = he.from();
      break;
    }
  }

  auto vertex_theta = OpenMesh::VProp<OpenMesh::Scalar>(m_mesh, "theta");
  vertex_theta[start.from()] = 0;

  float total_theta = 0;

  auto boundaryhe = start;
  while (boundaryhe.to() != start.from()) { // TODO: Check if this is iterating through the boundary in a continous way
    auto length = m_mesh.calc_edge_sqr_length(boundaryhe);
    auto theta = length * (2 * M_PI / total_length);

    total_theta += theta;
    vertex_theta[boundaryhe.to()] = total_theta;

    boundaryhe = m_mesh.next_halfedge_handle(boundaryhe);
  } 
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
  auto boundary_theta = OpenMesh::VProp<OpenMesh::Scalar>(m_mesh, "theta");
  
  gmm::dense_matrix<float> W(m_mesh.n_vertices(), m_mesh.n_vertices());
  gmm::clear(W);

  std::vector<float> bx(m_mesh.n_vertices(), 0);
  std::vector<float> by(m_mesh.n_vertices(), 0);

  for (const auto& vi : m_mesh.vertices()) {
    if (m_mesh.is_boundary(vi)) {
      assert(0 <= vi.idx && vi.idx < m_mesh.n_vertices());
      W(vi.idx(), vi.idx()) = 1.0f;

      bx[vi.idx] = cos(boundary_theta[vi]);
      by[vi.idx] = sin(boundary_theta[vi]);
    }
    else {
      W(cur_idx, cur_idx) = -m_mesh.valence(vi);
    }
  }

  for (const auto& ei : m_mesh.edges()) {
    auto v0 = ei.v0();
    auto v1 = ei.v1();
    W(index[v0], index[v1]) = 1.0f;
    W(index[v1], index[v0]) = 1.0f;
  }

  std::vector<float> x(m_mesh.n_vertices());
  std::vector<float> y(m_mesh.n_vertices());

  this->solve_linear_system(W, bx, x);
  this->solve_linear_system(W, by, y);

  for (const auto& vi : m_mesh.vertices()) {
    vparam_u[vi][0] = x[vi.idx];
    vparam_u[vi][1] = y[vi.idx];
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
