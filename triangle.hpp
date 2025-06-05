#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "node.hpp"

/**
 * @class Triangle
 * @brief Represents a triangular element in the mesh.
 * 
 * A triangle is defined by three node pointers and has an associated element id.
 * It provides methods to check point containment, perform splitting, and access
 * edge and vertex information.
 */
class Triangle {
 private:
  int id;  // the element id
  std::vector<Node *>
      vertices;  // stores pointers to the nodes that make the verteces of the triangle

  bool isValidElement();

 public:
  Triangle(std::vector<Node *> nodes) : id(-1), vertices(nodes) {}
  Triangle(int elem_id, std::vector<Node *> nodes) : id(elem_id), vertices(nodes) {}
  int get_id() { return id; };
  void set_id(int new_id) { id = new_id; };
  Node * get_v1() { return vertices[0]; }
  Node * get_v2() { return vertices[1]; }
  Node * get_v3() { return vertices[2]; }
  bool isInside(Node & interestNode);
  std::vector<Triangle> splitElement(Node & interestNode, int & elem_id_count);
  bool hasEdge(int id1, int id2) const;
  Node * getOppositeVertex(Node & a, Node & b) const;
};
