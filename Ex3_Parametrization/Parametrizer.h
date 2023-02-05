#ifndef Parametrizer_h
#define Parametrizer_h

#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Mesh/PolyMeshT.hh>

#include <gmm/gmm.h>


typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

//== CLASS DEFINITION =========================================================

class Parametrizer
{
public:
  enum ParameterizationMode {NoParameterization, Uniform, Harmonic};
  
  Parametrizer(MyMesh &iMesh);
  ~Parametrizer() {};
  
public:
  // calculate uniform parametrization
  void calc_uniform_parameterization();
  
  // calculate harmonic parametrization
  void calc_harmonic_parameterization();
  void calc_harmonic_weights();
  
  // calculate parametrizaton distortion
  void calc_distortion(ParameterizationMode imode);
  
  void computeTextureCoordinates(int iTextureWidth, int iTextureHeight, int iReapeats, ParameterizationMode imode);
  
  // flatten the mesh using the computed parametrization
  void flattenMesh(ParameterizationMode imode);
  
private:
  // solve linear system A * x = b
  void solve_linear_system(gmm::dense_matrix<float> & iA, std::vector<float>& iB, std::vector<float> & oX);

  void set_boundary();
  
private:
  MyMesh &  m_mesh;
};


//=============================================================================
#endif // Parametrizer_h defined
//=============================================================================



