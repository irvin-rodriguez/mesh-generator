#ifndef NODE_HPP
#define NODE_HPP

/**
 * @struct Point
 * @brief Contains the x and y coordinates of specific location.
 **/
struct Point {
  double x, y;
};

/**
 * @class Node
 * @brief Represents a mesh node with an id, coordinates, and boundary status.
 * 
 * Each node stores its unique id, position as a Point, and whether it's located
 * on the domain boundary.
 */
class Node {
 private:
  int id;            // node id
  Point coord;       // coordinates of the node
  bool on_boundary;  // whether or not the node lies on the domain boundary

 public:
  Node(int node_id, double x_coord, double y_coord, bool boundary) :
      id(node_id), on_boundary(boundary) {
    coord.x = x_coord;
    coord.y = y_coord;
  }
  Point get_coord() { return coord; };
  int get_id() { return id; };
  bool onBoundary() { return on_boundary; };
};

#endif
