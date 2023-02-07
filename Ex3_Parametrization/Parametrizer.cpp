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
    // std::cout << "(" << lastVertex.idx() << ", " << it->idx() << ") Length: " << length << " - Theta: " << theta << " - Total Theta: " << total_theta << std::endl;
    total_theta += theta; 
    vertex_theta[*it] = total_theta;
    lastVertex = *it;
  }

  auto last_length = (m_mesh.point(boundaries.front()) - m_mesh.point(boundaries.back())).length();
  auto theta = last_length * ( (2 * M_PI) / total_length);
  // std::cout << "(" << boundaries.back().idx() << ", " << boundaries.front().idx() << ") Length: " << last_length << " Total Theta: " << total_theta << std::endl;
}

float Parametrizer::list_boundary(std::vector<OpenMesh::VertexHandle> & oBoundary) {
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
      // std::cout << "Found first boundary vertex: "<< v_iter -> idx() <<"\n";
      break;
    }
  }


  bool finished = false;
  MyMesh::VertexVertexIter vv_iter;
  while(!finished) {
    finished = true;
    for(vv_iter = m_mesh.vv_iter(currentVertex); vv_iter.is_valid(); ++vv_iter) {
      if(m_mesh.is_boundary(*vv_iter)) {
        // std::cout << "Adding vertex to boundary. Id:" << (*vv_iter).idx() << "\n";
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

  // std::cout << "Total boundary length: " << total_length << std::endl; 
  return total_length;
}

// float Parametrizer::list_boundary(std::vector<OpenMesh::VertexHandle> & oBoundary) {
//   int nb_v = m_mesh.n_vertices();
//   bool found = true;
//   float total_length = 0;
  
//   MyMesh::VertexIter v_it;
//   MyMesh::VertexVertexIter vv_it;

//   for(v_it = m_mesh.vertices().begin(); v_it != m_mesh.vertices_end(); ++v_it) {
//     if(m_mesh.is_boundary(*v_it)){
//       oBoundary.push_back(*v_it);
//       break;
//     }
//   }

//   while (found){
//     found = false;
//     for (vv_it = m_mesh.vv_iter(oBoundary.back()); vv_it.is_valid(); ++v_it){
//       if ((m_mesh.is_boundary(*vv_it)) && (std::find(oBoundary.begin(), oBoundary.end(), *vv_it) == oBoundary.end())){
//         float length = (m_mesh.point(*vv_it) - m_mesh.point(oBoundary.back())).length();
//         total_length += length;
//         oBoundary.push_back(*vv_it);
//         found = true;
//         break;
//       }
//     }
//   }

//   return total_length;
// }





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
  
  for (const auto& ei : m_mesh.edges()) {
    auto v0 = ei.v0();
    auto v1 = ei.v1();
    W(v0.idx(), v1.idx()) = 1.0f;
    W(v1.idx(), v0.idx()) = 1.0f;
  }

  for (const auto& vi : m_mesh.vertices()) {
    if (m_mesh.is_boundary(vi)) {
      assert(0 <= (size_t)vi.idx() && (size_t)vi.idx() < m_mesh.n_vertices());
      W(vi.idx(), vi.idx()) = 1.0f;

      bx[vi.idx()] = cos(boundary_theta[vi]);
      by[vi.idx()] = sin(boundary_theta[vi]);
    }
    else {
      // std::cout << "Valence of " << vi.idx() << ": " << -m_mesh.valence(vi) << std::endl;
      // int valence = 0;
      for (const auto& vj : m_mesh.vertices()) {
        if (vi.idx() != vj.idx()) {
          W(vi.idx(), vi.idx()) -= W(vi.idx(), vj.idx());
        }
      }
      // for (const auto& eh : vi.edges()) valence++;
      // W(vi.idx(), vi.idx()) = -valence;
      // -m_mesh.valence(vi)
    }
  }


  //for(int j = 0; j < m_mesh.n_vertices(); j++) {
  //  for(int i = 0; i < m_mesh.n_vertices(); i++){ 
  //    std::cout << W(i, j) << " ";
  //  }
  //  std::cout << "\n";
  //}

  // for(int j = 0; j < m_mesh.n_vertices(); j++) {
    //std::cout << "(" << bx[j] << ", " << by[j] << ")\n"; 
  // }

  std::vector<float> x(m_mesh.n_vertices());
  std::vector<float> y(m_mesh.n_vertices());

  this->solve_linear_system(W, bx, x);
  this->solve_linear_system(W, by, y);

  for (const auto& vi : m_mesh.vertices()) {
    // std::cout << "vi: " << vi.idx() << " (" << x[vi.idx()] << ", " << y[vi.idx()] << ")\n"; 
    vparam_u[vi][0] = x[vi.idx()];
    vparam_u[vi][1] = y[vi.idx()];
  }
}

//TODO: delete
std::string point_to_str (OpenMesh::DefaultTraits::Point p) {
  return "(" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) + ")";
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

  auto vparam_h = OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_h");
  auto index = OpenMesh::VProp<MyMesh::Scalar>(m_mesh, "index");
  auto boundary_theta = OpenMesh::VProp<MyMesh::Scalar>(m_mesh, "theta");
  
  gmm::dense_matrix<float> W(m_mesh.n_vertices(), m_mesh.n_vertices());
  gmm::clear(W);

  std::vector<float> bx(m_mesh.n_vertices(), 0);
  std::vector<float> by(m_mesh.n_vertices(), 0);
  
  for (const auto& ei : m_mesh.edges()) {
    auto v0 = ei.v0();
    auto v1 = ei.v1();
    auto p0 = m_mesh.point(v0);
    auto p1 = m_mesh.point(v1);
    auto he0 = ei.h0();
    auto he1 = ei.h1();
    
    auto he_alpha = he0.next();
    auto he_beta = he1.next();
    auto v_alpha = he_alpha.to();
    auto v_beta = he_beta.to();
    auto p_alpha = m_mesh.point(v_alpha);
    auto p_beta = m_mesh.point(v_beta);

    // std::cout << "v0: " << point_to_str(p0) 
    //   << " v1: " << point_to_str(p1) 
    //   << " v_alpha: " << point_to_str(p_alpha) 
    //   << " v_beta: " << point_to_str(p_beta) << "\n";

    auto alpha_vec1 = p0 - p_alpha; 
    auto alpha_vec2 = p1 - p_alpha; 
    auto beta_vec1 = p0 - p_beta; 
    auto beta_vec2 = p1 - p_beta; 
    alpha_vec1 = alpha_vec1 / alpha_vec1.length();
    alpha_vec2 = alpha_vec2 / alpha_vec2.length();
    beta_vec1 = beta_vec1 / beta_vec1.length();
    beta_vec2 = beta_vec2 / beta_vec2.length();

    auto cos_alpha = alpha_vec1.dot(alpha_vec2);
    auto cos_beta = beta_vec1.dot(beta_vec2);
    auto sin_alpha = alpha_vec1.cross(alpha_vec2).length();
    auto sin_beta = beta_vec1.cross(beta_vec2).length();
    auto cot_alpha = cos_alpha / sin_alpha;
    auto cot_beta = cos_beta / sin_beta;
    auto wij = (cot_alpha + cot_beta) / 2;

    W(v0.idx(), v1.idx()) = wij;
    W(v1.idx(), v0.idx()) = wij;
  }

  for (const auto& vi : m_mesh.vertices()) {
    if (m_mesh.is_boundary(vi)) {
      assert(0 <= (size_t)vi.idx() && (size_t)vi.idx() < m_mesh.n_vertices());
      W(vi.idx(), vi.idx()) = 1.0f;

      bx[vi.idx()] = cos(boundary_theta[vi]);
      by[vi.idx()] = sin(boundary_theta[vi]);
    }
    else {
      // std::cout << "Valence of " << vi.idx() << ": " << -m_mesh.valence(vi) << std::endl;
      // int valence = 0;
      for (const auto& vj : m_mesh.vertices()) {
        if (vi.idx() != vj.idx()) {
          W(vi.idx(), vi.idx()) -= W(vi.idx(), vj.idx());
        }
      }
      // for (const auto& eh : vi.edges()) valence++;
      // W(vi.idx(), vi.idx()) = -valence;
      // -m_mesh.valence(vi)
    }
  }


  // for(int j = 0; j < m_mesh.n_vertices(); j++) {
  //  for(int i = 0; i < m_mesh.n_vertices(); i++){ 
  //    std::cout << W(i, j) << " ";
  //  }
  //  std::cout << "\n";
  // }

  // for(int j = 0; j < m_mesh.n_vertices(); j++) {
    //std::cout << "(" << bx[j] << ", " << by[j] << ")\n"; 
  // }

  std::vector<float> x(m_mesh.n_vertices());
  std::vector<float> y(m_mesh.n_vertices());

  this->solve_linear_system(W, bx, x);
  this->solve_linear_system(W, by, y);

  for (const auto& vi : m_mesh.vertices()) {
    // std::cout << "vi: " << vi.idx() << " (" << x[vi.idx()] << ", " << y[vi.idx()] << ")\n"; 
    vparam_h[vi][0] = x[vi.idx()];
    vparam_h[vi][1] = y[vi.idx()];
  }

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

float angle3D(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2) {
  return acos(v1.dot(v2) / (v1.length() * v2.length()));
}

float angle2D(OpenMesh::Vec2f v1, OpenMesh::Vec2f v2) {
  return acos(v1.dot(v2) / (v1.length() * v2.length()));
}

float triangleAreaFromSides(float s1, float s2, float s3) {
  auto s = (s1 + s2 + s3) / 2;
  return sqrt(s * (s - s1) * (s - s2) * (s - s3));
}

float triangleArea3D(OpenMesh::Vec3f p0, OpenMesh::Vec3f p1, OpenMesh::Vec3f p2) {
  auto s1 = (p1 - p0).length();
  auto s2 = (p2 - p1).length();
  auto s3 = (p0 - p2).length();
  return triangleAreaFromSides(s1, s2, s3);
} 

float triangleArea2D(OpenMesh::Vec2f p0, OpenMesh::Vec2f p1, OpenMesh::Vec2f p2) {
  auto s1 = (p1 - p0).length();
  auto s2 = (p2 - p1).length();
  auto s3 = (p0 - p2).length();
  return triangleAreaFromSides(s1, s2, s3);
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
  auto vparam_h = OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_h");
  auto vparam_u = OpenMesh::VProp<OpenMesh::Vec2f>(m_mesh, "vparam_u");
  auto v_param = vparam_h;
  if (imode == Uniform) {
    v_param = vparam_u;
  }

  float totalArea3D = 0.0f;
  float totalArea2D = 0.0f;
  for (const auto& fh : m_mesh.faces()) {
    std::vector<OpenMesh::SmartVertexHandle> vhs;
    for (const auto& vh : fh.vertices()) {
      vhs.push_back(vh);
    }
    std::vector<float> angles3D;
    assert(vhs.size() == 3);
    auto p3D0 = m_mesh.point(vhs[0]);
    auto p3D1 = m_mesh.point(vhs[1]);
    auto p3D2 = m_mesh.point(vhs[2]);
    auto p2D0 = v_param[vhs[0]];
    auto p2D1 = v_param[vhs[1]];
    auto p2D2 = v_param[vhs[2]];
    totalArea3D += triangleArea3D(p3D0, p3D1, p3D2);
    totalArea2D += triangleArea2D(p2D0, p2D1, p2D2);
  }

  for (const auto& fh : m_mesh.faces()) {
    std::vector<OpenMesh::SmartVertexHandle> vhs;
    for (const auto& vh : fh.vertices()) {
      vhs.push_back(vh);
    }
    std::vector<float> angles3D;
    assert(vhs.size() == 3);
    auto p3D0 = m_mesh.point(vhs[0]);
    auto p3D1 = m_mesh.point(vhs[1]);
    auto p3D2 = m_mesh.point(vhs[2]);
    auto vec3Dv0v1 = p3D1 - p3D0;
    auto vec3Dv1v2 = p3D2 - p3D1;
    auto vec3Dv2v0 = p3D0 - p3D2;
    angles3D.push_back(angle3D(vec3Dv0v1, -vec3Dv2v0));
    angles3D.push_back(angle3D(-vec3Dv0v1, vec3Dv1v2));
    angles3D.push_back(angle3D(-vec3Dv1v2, vec3Dv2v0));
    auto area3D = triangleArea3D(p3D0, p3D1, p3D2);

    std::vector<float> angles2D;
    auto p2D0 = v_param[vhs[0]];
    auto p2D1 = v_param[vhs[1]];
    auto p2D2 = v_param[vhs[2]];
    auto vec2Dv0v1 = p2D1 - p2D0;
    auto vec2Dv1v2 = p2D2 - p2D1;
    auto vec2Dv2v0 = p2D0 - p2D2;
    angles2D.push_back(angle2D(vec2Dv0v1, -vec2Dv2v0));
    angles2D.push_back(angle2D(-vec2Dv0v1, vec2Dv1v2));
    angles2D.push_back(angle2D(-vec2Dv1v2, vec2Dv2v0));
    auto area2D = triangleArea2D(p2D0, p2D1, p2D2);

    for (int i = 0; i < 3; i++) {
      auto angle_diff = angles2D[i] - angles3D[i];
      auto area_diff = (area3D / totalArea3D) - (area2D / totalArea2D);
      angle_distortion += angle_diff * angle_diff;
      area_distortion += area_diff * area_diff;
    }
  }

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
