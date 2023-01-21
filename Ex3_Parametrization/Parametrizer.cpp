//=============================================================================
//
//  Implementation of class Parametrizer
//
//=============================================================================

#include <iostream>

#include "Parametrizer.h"

Parametrizer::Parametrizer(MyMesh &iMesh):
m_mesh(iMesh)
{
  OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_u");
  OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_h");
}

void Parametrizer::solve_linear_system(gmm::dense_matrix<float> & iA, std::vector<float>& iB, std::vector<float> & oX)
{
  // for better performances, we should switch to sparse matrices and use solvers for
  // sparse systems but this requires the use of external third-party libraries...
  oX.clear();
  oX.resize(iB.size());
  gmm::lu_solve(iA, oX, iB);
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
