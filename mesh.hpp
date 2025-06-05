#include <cmath>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "node.hpp"
#include "triangle.hpp"

typedef std::vector<std::pair<double, double> > surface;

/**
 * @class Mesh
 * @brief Represents a rectangular mesh.
 *
 * The user constructs their mesh directly through this class and its public functions.
 * The inteded workflow in using this class to make a mesh is to first construct an instance
 * of Mesh by properly supplying a formatted input file, then creating the nodes, then
 * meshing the geometry through triangulation. Lastly, to view the created mesh, the user has
 * print functions.
 **/
class Mesh {
 private:
  double hx, hy;            // stores the max element size in the x and y directions
  surface corners;          // stores the corners that create the rectangular geometry
  std::vector<Node> nodes;  // stores the nodes of the mesh
  std::vector<Triangle> elements;  // stores the elements of the mesh

  void parseFile(const char * filename);
  void isRectangle();
  void performTriangulation(bool ifDelaunay);
  void replaceExistingElement(Triangle & triangle, std::vector<Triangle> & new_elements);
  void delaunayElementCheck(std::vector<Triangle> & new_elements,
                            Node & inserted_node,
                            Triangle & original_element);
  std::vector<int> findOppositeEdge(Triangle & triangle, Node & inserted_node);
  Node * getNodeById(int id);

 public:
  Mesh(const char * filename);
  void createNodes();
  void createNodes(double radius);
  void simpleTriangulation();
  void delaunayTriangulation();

  void printNodes();
  void printElements();
  void printMeshInfo();
};
