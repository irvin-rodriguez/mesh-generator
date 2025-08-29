#include "triangle.hpp"

/**
 * @brief Checks if a node is inside the triangle or lies on its boundary.
 * 
 * This method determines whether a given node lies strictly inside the triangle
 * using barycentric coordinates. However if the node is marked as being on the boundary,
 * it checks if it lies exactly on one of the triangle's edges (vertical or horizontal).
 *
 * @param interestNode The node to test for inclusion in the triangle.
 * @return true if the node is strictly inside the triangle or lies on its edge; false otherwise.
 * 
 * @note The specific check for nodes on the boundary is to avoid falsely excluding nodes that are exactly 
 * on the boundary due to rounding, just to be extra safe.
 **/
bool Triangle::isInside(Node & interestNode) {
  // Get the coordinates of the node we want to check
  Point p = interestNode.get_coord();

  // Special handling for nodes on the domain boundary
  if (interestNode.onBoundary()) {
    // Get the endpoints of one edge of the triangle
    for (int i = 0; i < 3; ++i) {
      Point vi = vertices[i]->get_coord();
      Point vj = vertices[(i + 1) % 3]->get_coord();  // wrap around

      // Check for vertical edge: x-coordinates must match
      if (vi.x == vj.x && p.x == vi.x) {
        // Check if the point's y coord lies between the endpoints' y coords
        if ((p.y >= std::min(vi.y, vj.y)) && (p.y <= std::max(vi.y, vj.y))) {
          return true;
        }
      }

      // Check for horizontal edge: y-coordinates must match
      if (vi.y == vj.y && p.y == vi.y) {
        // Check if the point's x coord lies between the endpoints' x coords
        if ((p.x >= std::min(vi.x, vj.x)) && (p.x <= std::max(vi.x, vj.x))) {
          return true;
        }
      }
    }
  }

  // For all nodes (not just boundaries), use barycentric coordinates
  // to determine if the point is inside the triangle

  // Get coordinates of the triangle vertices
  Point v1 = get_v1()->get_coord();
  Point v2 = get_v2()->get_coord();
  Point v3 = get_v3()->get_coord();

  // Compute the denominator of the barycentric coordinate formulas
  double denom = (v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y);

  // Compute barycentric coordinates (alpha, beta, gamma)
  double alpha = ((v2.y - v3.y) * (p.x - v3.x) + (v3.x - v2.x) * (p.y - v3.y)) / denom;
  double beta = ((v3.y - v1.y) * (p.x - v3.x) + (v1.x - v3.x) * (p.y - v3.y)) / denom;
  double gamma = 1.0 - alpha - beta;

  // If all barycentric coordinates are non-negative, the point lies inside
  // or on the edge of the triangle (inclusive check)
  return (alpha >= 0) && (beta >= 0) && (gamma >= 0);
}

/**
 * @brief Checks if the supplied vertices create a valid triangle.
 *
 * This method checks if the vertices which were used to construct
 * an instance of Triangle, actually forms a valid triangle. A valid
 * triangle must have an area larger than 0.
 *
 * @return true if the Triangle is valid, false if the Triangle is not.
 */
bool Triangle::isValidElement() {
  // Grab the coordinates of the of the vertices
  Point v1 = get_v1()->get_coord();
  Point v2 = get_v2()->get_coord();
  Point v3 = get_v3()->get_coord();

  // Compute the area of the triangle
  double area =
      0.5 * std::abs(v1.x * (v2.y - v3.y) + v2.x * (v3.y - v1.y) + v3.x * (v1.y - v2.y));

  // Check if the area is not equal to 0
  if (area != 0) {
    return true;
  }
  return false;
}

/**
 * @brief Given a Node inside the Triangle, split the Triangle into smaller ones.
 *
 * Given a node that is inside the Triangle, this method will split the Triangle
 * into smaller Triangles that together form the original Triangle. This method
 * does not handle replacement of the original Triangle, except it only returns 
 * the Triangles which it creates. As the smaller triangles are created, they are
 * being checked if they create valid triangles. They are also assigned their 
 * Triangle id, ensuring that one Triangle created has the same id as the original
 * Triangle that was split.
 * 
 * @param interestNode The Node that within the Triangle that will be used to split.
 * @param elem_id_count The global element id count within the Mesh.
 * @return A vector of Triangles that were created by splitting the original Triangle.
 * 
 **/
std::vector<Triangle> Triangle::splitElement(Node & interestNode, int & elem_id_count) {
  // Ensure Node given is inside the triangle
  if (!isInside(interestNode)) {
    std::cerr << "Error: the node given is not within the element we are trying to split"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Create a vector tor store new elements created
  std::vector<Triangle> new_elements;

  // Grab the vertices of the element we are splitting
  Node * v1 = vertices[0];
  Node * v2 = vertices[1];
  Node * v3 = vertices[2];

  // Create first element
  std::vector<Node *> new_verticies;
  new_verticies.push_back(v1);
  new_verticies.push_back(v2);
  new_verticies.push_back(&interestNode);
  Triangle t1(new_verticies);
  // Check if its a valid element
  if (t1.isValidElement()) {
    // Insert it into the vector of elements
    new_elements.push_back(t1);
  }
  // Clear the vector of vertices to prepare for the next one
  new_verticies.clear();

  // Create the second element
  new_verticies.push_back(v2);
  new_verticies.push_back(v3);
  new_verticies.push_back(&interestNode);
  Triangle t2(new_verticies);
  // Check that its a valid element
  if (t2.isValidElement()) {
    // Insert it into the vector of elements
    new_elements.push_back(t2);
  }
  // Clear the vector of vertices to prepare for the next one
  new_verticies.clear();

  // Create the third element
  new_verticies.push_back(v1);
  new_verticies.push_back(v3);
  new_verticies.push_back(&interestNode);
  Triangle t3(new_verticies);
  // Check that its a valid element
  if (t3.isValidElement()) {
    // Check that its a valid element
    new_elements.push_back(t3);
  }

  // Assign the new elements their ids
  for (size_t i = 0; i < new_elements.size(); ++i) {
    // Ensure the first one has the same element id as the original Triangle
    if (i == 0) {
      new_elements[i].set_id(id);
    }
    // Other Triangles, can have the next new element id in the Mesh
    else {
      new_elements[i].set_id(elem_id_count);
      ++elem_id_count;
    }
  }

  return new_elements;
}

/**
 * @brief Given the the id of two nodes that form an edge, check if the Triangle contains that same edge.
 * 
 * @param id1 The id of a Node.
 * @param id2 The id of another Node.
 * @return true if the Triangle contains the same edge, false if it doesn't
 * 
 **/
bool Triangle::hasEdge(int id1, int id2) const {
  int a = vertices[0]->get_id();
  int b = vertices[1]->get_id();
  int c = vertices[2]->get_id();

  return (a == id1 && b == id2) || (b == id1 && a == id2) || (b == id1 && c == id2) ||
         (c == id1 && b == id2) || (c == id1 && a == id2) || (a == id1 && c == id2);
}

/**
 * @brief Given two nodes in the Triangle, return a pointer to the third Node.
 * 
 * @param a A referece to a Node in the Triangle
 * @param b A reference to the other Node in the Triangle
 * @return A pointer to the third Node in the Triangle.
 * 
 **/
Node * Triangle::getOppositeVertex(Node & a, Node & b) const {
  // Loop through the Nodes in the Triangle
  for (size_t i = 0; i < 3; ++i) {
    Node * v = vertices[i];
    // If Node doesn't share the id of the Nodes we passed in, this is the opposite vertex
    if (v->get_id() != a.get_id() && v->get_id() != b.get_id()) {
      return v;
    }
  }

  std::cerr << "Error: a node given is not within the element" << std::endl;
  exit(EXIT_FAILURE);
}
