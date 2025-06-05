#include "mesh.hpp"

int main(int argc, char ** argv) {
  // open file for reading
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return EXIT_FAILURE;
  }

  Mesh msh(argv[1]);  // create mesh instance from filename
  msh.createNodes();
  msh.printNodes();           // print the nodes
  msh.simpleTriangulation();  // meshes the domain using simple triangulation
  msh.printElements();        // print the elements

  return EXIT_SUCCESS;
}
