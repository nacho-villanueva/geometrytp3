#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include "Parametrizer.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

int main(int argc, char **argv)
{
  // check command line options
  if (argc == 1 || (argc == 2 && strcmp(argv[1], "-h")==0))
  {
    // display manual
    std::cout << "Usage:  " << argv[0] << " inputfile outputfile mode outputformat [width] [height] [repeats]\n";
    std::cout << "Valid options for mode are: \n";
    std::cout << " -u:" << " compute uniform parametrization\n";
    std::cout << " -h:" << " compute harmonic parametrization\n";
    std::cout << "Valid options for outputformat are: \n";
    std::cout << " -f:" << " flattened mesh \n";
    std::cout << " -t:" << " 3D mesh with texture coordinates\n";
    
    return 1;
  }
  else if (argc !=5 && (argc !=8 && strcmp(argv[4], "-t")==0))
  {
    std::cerr << "Usage:  " << argv[0] << " inputfile outputfile mode outputformat [width] [height] [repeats]\n";
    return 1;
  }
  
  MyMesh  mesh;
  
  // read mesh from stdin
  if ( ! OpenMesh::IO::read_mesh(mesh, argv[1]) )
  {
    std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
    return 1;
  }
  
  mesh.request_vertex_texcoords2D();
  
  // file options
  OpenMesh::IO::Options wopt;
  
  Parametrizer::ParameterizationMode parametrizationMode = Parametrizer::NoParameterization;
  
  Parametrizer parametrizer(mesh);
  
  char * option = argv[3];
  if (strcmp(option, "-u")==0)
  {
    parametrizationMode = Parametrizer::Uniform;
    parametrizer.calc_uniform_parameterization();
    parametrizer.calc_distortion(parametrizationMode);
  }
  else if (strcmp(option, "-h")==0)
  {
    parametrizationMode = Parametrizer::Harmonic;
    parametrizer.calc_harmonic_parameterization();
    parametrizer.calc_distortion(parametrizationMode);
  }
  else
  {
    std::cerr << "Invalid option" << std::endl;
    std::cout << "Valid options for mode are: \n";
    std::cout << " -u:" << " compute uniform parametrization\n";
    std::cout << " -h:" << " compute harmonic parametrization\n";
    return 1;
  }
                  
  char * output_mode = argv[4];
  if (strcmp(output_mode, "-f")==0)
  {
    parametrizer.flattenMesh(parametrizationMode);
  }
  else if (strcmp(output_mode, "-t")==0)
  {
    int texture_width = atoi(argv[5]);
    int texture_height = atoi(argv[6]);
    int repeats = atoi(argv[7]);
    parametrizer.computeTextureCoordinates(texture_width, texture_height, repeats, parametrizationMode);
    wopt = OpenMesh::IO::Options::VertexTexCoord;
  }
  else
  {
    std::cerr << "Invalid option" << std::endl;
    std::cout << "Valid options for outputformat are: \n";
    std::cout << " -f:" << " flattened mesh \n";
    std::cout << " -t:" << " 3D mesh with texture coordinates\n";
    return 1;
  }
  
  // write mesh to stdout
  if ( ! OpenMesh::IO::write_mesh(mesh, argv[2], wopt) )
  {
    std::cerr << "Error: cannot write mesh to " << argv[2] << std::endl;
    return 1;
  }
  
  return 0;
}


