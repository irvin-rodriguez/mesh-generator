#include "mesh.hpp"

/**
 * @brief Within a circular region, defined by a radius, provide a random dx and dy to shift a point by.
 *
 * @param radius The radius of the circular region you wish to possibly displace the point by.
 * 
 * @return A pair of doubles indicating the dx and dy values.
 * 
 **/
std::pair<double, double> randomShiftInRadius(double radius) {
  // Generate random seed
  static bool seeded = false;
  if (!seeded) {
    std::srand(std::time(0));
    seeded = true;
  }

  // Compute random theta and radius within the region
  double theta = ((double)std::rand() / RAND_MAX) * 2.0 * M_PI;
  double r = radius * std::sqrt((double)std::rand() / RAND_MAX);

  // Compute the dx and dy to shift by
  double dx = r * std::cos(theta);
  double dy = r * std::sin(theta);

  return std::make_pair(dx, dy);
}

/**
 * @brief Count the number of lines in a file.
 *
 * @param file The file which will be counted.
 * 
 * @return An integer of the number of lines in the file.
 * 
 **/
int countLines(std::ifstream & file) {
  // For every line in the file, count it
  int count = 0;
  std::string line;
  while (std::getline(file, line)) {
    ++count;
  }

  // Reset the global variable that tracks where we are in the file
  file.clear();
  file.seekg(0, std::ios::beg);

  return count;
}

/**
 * @brief In a specific distance, compute the element size and number of elements that fit within that distance.
 *
 * @param distance The distance in a direction of the region we are working with.
 * @param max_h The maximum the element size can be. 
 * 
 * @return A pair containing the number of the elements that fit int the distance, and the size of the elements.
 * 
 **/
std::pair<int, double> computeElemSize(double distance, double max_h) {
  int n_elem = std::ceil(distance / max_h);
  double new_h = distance / n_elem;
  return std::make_pair(n_elem, new_h);
}

/**
 * @brief Parse an input file to obtain mesh parameters.
 *
 * @param filename The file we wish to parse through.
 * 
 * 
 **/
void Mesh::parseFile(const char * filename) {
  std::ifstream inputfile(filename);

  // Check that the file exists
  if (!inputfile) {
    std::cerr << "Error: Can't open file." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Ensure that the file has the correct number of lines desired for inputs
  if (countLines(inputfile) != 6) {
    std::cerr << "Error: Input file does not have proper amount of lines." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  size_t lineCount = 1;
  // Loop through every line in the file
  while (std::getline(inputfile, line)) {
    std::stringstream iss(line);

    // The first line must contain just a single values
    if (lineCount == 1) {
      if (!(iss >> hx) || !(iss >> std::ws).eof()) {
        std::cerr << "Error: hx in line 1 is not formatted properly." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    // The second line must contain a single values
    else if (lineCount == 2) {
      if (!(iss >> hy) || !(iss >> std::ws).eof()) {
        std::cerr << "Error: hy in line 2 is not formatted properly." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    // The third line must contain just two values
    else {
      double x, y;
      if (!(iss >> x >> y) || !(iss >> std::ws).eof()) {
        std::cerr << "Error: x and y coordinates for a vertex is not formatted properly."
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      corners.push_back(std::make_pair(x, y));
    }
    ++lineCount;
  }
}

/**
 * @brief Ensure that the corners provided create a valid Rectangle.
 *
 **/
void Mesh::isRectangle() {
  // Check that there are four corners
  if (corners.size() != 4) {
    std::cerr << "Error: there is not enough points to make a rectangle." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Extract the x and y points and put them into sets
  std::set<double> x_coords;
  std::set<double> y_coords;
  for (size_t i = 0; i < corners.size(); i++) {
    x_coords.insert(corners[i].first);
    y_coords.insert(corners[i].second);
  }

  // Check that the sets contain two values indicating they were unique
  if (x_coords.size() != 2 || y_coords.size() != 2) {
    std::cerr
        << "Error: a valid rectangle needs to have 2 unique x and 2 unique y values."
        << std::endl;
    exit(EXIT_FAILURE);
  }

  // Knowing they are unique, we must ensure they form the expected rectangle
  double x_min = *x_coords.begin();
  double x_max = *x_coords.rbegin();
  double y_min = *y_coords.begin();
  double y_max = *y_coords.rbegin();

  std::set<std::pair<double, double> > expected_corners;
  expected_corners.insert(std::make_pair(x_min, y_min));
  expected_corners.insert(std::make_pair(x_min, y_max));
  expected_corners.insert(std::make_pair(x_max, y_min));
  expected_corners.insert(std::make_pair(x_max, y_max));

  // This is our actual rectangle
  std::set<std::pair<double, double> > actual_corners(corners.begin(), corners.end());

  // Compare them to ensure we have a rectangle
  if (actual_corners != expected_corners) {
    std::cerr << "Error: the corners provided do not form a rectangle." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Order the corners so its in the right order
  corners.clear();
  corners.push_back(std::make_pair(x_min, y_min));
  corners.push_back(std::make_pair(x_min, y_max));
  corners.push_back(std::make_pair(x_max, y_min));
  corners.push_back(std::make_pair(x_max, y_max));
}

/**
 * @brief Creates a uniform set of Nodes in the domain.
 * 
 * The way the nodes are created follow this order. The first four
 * nodes are the corners of the domain, then from the lowest leftmost
 * corner we begin adding nodes bottom up and left to right.
 *
 **/
void Mesh::createNodes() {
  // If the nodes are not empty, clear them to prepare for the new nodes
  if (!nodes.empty()) {
    nodes.clear();
  }

  // Compute the distance
  double x_distance = std::abs(corners[3].first - corners[0].first);
  double y_distance = std::abs(corners[3].second - corners[0].second);

  // Compute the element size and number of elements
  std::pair<int, double> x_result = computeElemSize(x_distance, hx);
  std::pair<int, double> y_result = computeElemSize(y_distance, hy);

  int nx_elem = x_result.first;
  int ny_elem = y_result.first;
  hx = x_result.second;
  hy = y_result.second;

  // Add corners into nodes
  int id = 1;
  for (int i = 0; i < 4; i++) {
    nodes.push_back(Node(id, corners[i].first, corners[i].second, true));
    ++id;
  }

  // Define the origin
  Node origin = nodes[0];
  // Performing algorithm to place the nodes in the domain
  for (int i = 0; i <= nx_elem; i++) {
    double x = origin.get_coord().x + i * hx;
    for (int j = 0; j <= ny_elem; j++) {
      double y = origin.get_coord().y + j * hy;

      // Check if the node is at the corner if so don't place it
      bool is_corner_node = (i == 0 && j == 0) || (i == 0 && j == ny_elem) ||
                            (i == nx_elem && j == 0) || (i == nx_elem && j == ny_elem);
      if (is_corner_node) {
        continue;
      }
      // Check if the node is on the boundary, if so identify that
      bool on_boundary = (i == 0 || i == nx_elem || j == 0 || j == ny_elem);
      nodes.push_back(Node(id, x, y, on_boundary));

      ++id;
    }
  }
}

/**
 * @overload
 * @brief Creates randomly pertubed nodes with a circular region.
 * 
 * @param radius The radius of the circular region you wish you to possibly
 * pertube a Node by. This radius must be smaller than the element size
 * in any given direction.
 * 
 */
void Mesh::createNodes(double radius) {
  // If the nodes are not empty, clear them to prepare for the new nodes
  if (!nodes.empty()) {
    nodes.clear();
  }

  // Compute the distance
  double x_distance = std::abs(corners[3].first - corners[0].first);
  double y_distance = std::abs(corners[3].second - corners[0].second);

  // Compute the element size and number of elements
  std::pair<int, double> x_result = computeElemSize(x_distance, hx);
  std::pair<int, double> y_result = computeElemSize(y_distance, hy);

  int nx_elem = x_result.first;
  int ny_elem = y_result.first;
  hx = x_result.second;
  hy = y_result.second;

  // Check that the radius is not larger than the element size in any direction
  if (radius > hx || radius > hy) {
    std::cerr
        << "Error: radius can not be larger than effective element size in any direction"
        << std::endl;
    exit(EXIT_FAILURE);
  }

  // Add corners into nodes
  int id = 1;
  for (int i = 0; i < 4; i++) {
    nodes.push_back(Node(id, corners[i].first, corners[i].second, true));
    ++id;
  }

  // Set the origin
  Node origin = nodes[0];
  // Performing algorithm to place the nodes in the domain
  for (int i = 0; i <= nx_elem; i++) {
    double x = origin.get_coord().x + i * hx;
    for (int j = 0; j <= ny_elem; j++) {
      double y = origin.get_coord().y + j * hy;

      // Check if the node is at the corner if so don't place it
      bool is_corner_node = (i == 0 && j == 0) || (i == 0 && j == ny_elem) ||
                            (i == nx_elem && j == 0) || (i == nx_elem && j == ny_elem);
      if (is_corner_node) {
        continue;
      }

      // Check if the node is on the boundary, if so identify that and do not perturb it
      bool on_boundary = (i == 0 || i == nx_elem || j == 0 || j == ny_elem);
      if (on_boundary) {
        nodes.push_back(Node(id, x, y, on_boundary));
      }
      // If the node is not on the boundary, identify that and perturb it
      else {
        std::pair<double, double> shift = randomShiftInRadius(radius);
        nodes.push_back(Node(id, x + shift.first, y + shift.second, on_boundary));
      }

      ++id;
    }
  }
}

/**
 * @brief Constructor for an instance of Mesh from an input file.
 * 
 * Parses an input file with mesh information, and validates that
 * we have a valid rectangular domain.
 *
 **/
Mesh::Mesh(const char * filename) {
  parseFile(filename);  // parse file
  isRectangle();        // check domain is a valid rectangle
}

/**
 * @brief Prints the element id and coordinates of the nodes of the mesh.
 *
 **/
void Mesh::printNodes() {
  std::cout << "$nodes" << std::endl;
  for (size_t i = 0; i < nodes.size(); i++) {
    Node & node = nodes[i];
    std::cout << node.get_id() << " " << node.get_coord().x << " " << node.get_coord().y
              << std::endl;
  }
}

/**
 * @brief Replaces an existing element in the mesh with another that shares the same element id.
 * 
 * This method serves to replace an existing element with an element labeled with the same element id.
 * This method is intended to be used after the existing element has been split into smaller elements.
 *
 **/
void Mesh::replaceExistingElement(Triangle & element,
                                  std::vector<Triangle> & new_elements) {
  int original_id = element.get_id();
  bool replaced = false;
  // For every new element check if it's id is the same as the element we wish to replace
  for (size_t i = 0; i < new_elements.size(); ++i) {
    // If so, replace that element and remove it from the vector of new elements
    if (new_elements[i].get_id() == original_id) {
      element = new_elements[i];                     // Replace in-place
      new_elements.erase(new_elements.begin() + i);  // Remove it from new_elements
      replaced = true;  // Indicate that the replacement has been done
      break;
    }
  }
  // If no replacment was done, notify that this was the case
  if (!replaced) {
    std::cerr
        << "Error: splitting element did not return a matching element id with original"
        << std::endl;
    exit(EXIT_FAILURE);
  }
}

/**
 * @brief For a Triangle and a Node that is a vertex for this triangle, return the edge across from that Node.
 * 
 * This method essentially checks the combinations of edges in the Triangle, and returns the
 * edge which sits across from Node of interest.
 * 
 * @param triangle This is the element we are working with.
 * @param inserted_node This is the node we are interested in finding its opposite edge.
 * 
 * @return A vector of Node id's representing an edge.
 *
 **/
std::vector<int> Mesh::findOppositeEdge(Triangle & triangle, Node & inserted_node) {
  int inserted_id = inserted_node.get_id();
  std::vector<int> edge;

  // Find the id of the vertices of the Triangle
  int id1 = triangle.get_v1()->get_id();
  int id2 = triangle.get_v2()->get_id();
  int id3 = triangle.get_v3()->get_id();

  // Find which vertex matches with the Node we are interested in
  // once found add this as an edge
  if (id1 == inserted_id) {
    edge.push_back(id2);
    edge.push_back(id3);
  }
  else if (id2 == inserted_id) {
    edge.push_back(id1);
    edge.push_back(id3);
  }
  else if (id3 == inserted_id) {
    edge.push_back(id1);
    edge.push_back(id2);
  }
  else {
    std::cerr << "Error: inserted node not found in triangle" << std::endl;
    exit(EXIT_FAILURE);
  }

  return edge;
}

/**
 * @brief Returns a node in the mesh based on the id given.
 * 
 * @param id The id of a node we wish to grab.
 * 
 * @return A pointer to a specific Node.
 *
 **/
Node * Mesh::getNodeById(int id) {
  // Ensure nodes have been computed
  if (nodes.empty()) {
    std::cerr << "Error: nodes have not yet been created yet" << std::endl;
    exit(EXIT_FAILURE);
  }
  // Check every node in the mesh for the one that matches the id given
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i].get_id() == id) {
      return &nodes[i];
    }
  }

  std::cerr << "Error: Node with ID " << id << " not found." << std::endl;
  exit(EXIT_FAILURE);
}

/**
 * @brief Computes the angle between 3 nodes.
 * 
 * Computes the angle formed between nodes a, b and c.
 *   a
 *  /
 * b --- c
 * 
 * @param a The first node.
 * @param b The second node.
 * @param c The third node.
 * 
 * @return A angle in radians.
 *
 **/
double angleBetween(Node * a, Node * b, Node * c) {
  // Compute the vector formed from the Nodes
  double ux = a->get_coord().x - b->get_coord().x;
  double uy = a->get_coord().y - b->get_coord().y;
  double vx = c->get_coord().x - b->get_coord().x;
  double vy = c->get_coord().y - b->get_coord().y;

  // Computes angle at b between points a-b-c
  double dot = ux * vx + uy * vy;
  double len_u = std::sqrt(ux * ux + uy * uy);
  double len_v = std::sqrt(vx * vx + vy * vy);
  double cos_theta = dot / (len_u * len_v);
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));  // place between [-1, 1]

  return std::acos(cos_theta);  // in radians
}

/**
 * @brief Identifies whether an edge should be flipped or not.
 * 
 * Identifies if an edge should be flipped based on whether the swap results in a 
 * larger minimum interior angle of the two triangles that comprise it. If it does, swap it.
 * 
 * @param a The first node.
 * @param b The second node.
 * @param c The third node.
 * @param d The fourth node.
 * 
 * @return true if the edge should be swapped, false if it shouldn't
 * 
 * @note Triangles: abc and bad (shared edge ab, opposite points c and d)
 *
 **/
bool shouldFlipEdge(Node * a, Node * b, Node * c, Node * d) {
  // Angles in original triangles
  double angle_c_ab = angleBetween(c, a, b);
  double angle_d_ba = angleBetween(d, b, a);
  double min_original = std::min(angle_c_ab, angle_d_ba);

  // Angles in proposed triangles: c-d-a and d-c-b (after flip to diagonal cd)
  double angle_a_cd = angleBetween(a, c, d);
  double angle_b_dc = angleBetween(b, d, c);
  double min_proposed = std::min(angle_a_cd, angle_b_dc);

  return min_proposed > min_original;
}

/**
 * @brief Performs an edge flip between two adjacent triangles.
 * 
 * This function replaces two adjacent triangles sharing edge ab with two new triangles
 * that instead share the opposite diagonal cd. It updates the input triangle objects
 * t1 and t2 using their original id's and new vertex order.
 * 
 * The resulting triangles are:
 * - Triangle 1: vertices c, d, a
 * - Triangle 2: vertices d, c, b
 * 
 * 
 * @param t1 Pointer to the first triangle.
 * @param t2 Pointer to the second triangle.
 * @param a First shared node of the original edge.
 * @param b Second shared node of the original edge.
 * @param c Opposite node from triangle t1.
 * @param d Opposite node from triangle t2.
 * 
 */
void flipEdge(Triangle * t1, Triangle * t2, Node * a, Node * b, Node * c, Node * d) {
  int id1 = t1->get_id();
  int id2 = t2->get_id();

  // Create new triangles with diagonal cd
  // First triangle: c-d-a
  std::vector<Node *> new_vertices;
  new_vertices.push_back(c);
  new_vertices.push_back(d);
  new_vertices.push_back(a);
  Triangle new_t1(id1, new_vertices);
  new_vertices.clear();

  // Second triangle: d-c-b
  new_vertices.push_back(d);
  new_vertices.push_back(c);
  new_vertices.push_back(b);
  Triangle new_t2(id2, new_vertices);
  new_vertices.clear();

  // Replace the old triangles
  *t1 = new_t1;
  *t2 = new_t2;
}

/**
 * @brief Checks if any of the new elements need to have their edge flipped.
 * 
 * For each new triangle formed after inserting a node, this checks if the opposite edge 
 * (from the inserted node) is shared with an existing triangle. If so, it forms a 
 * quadrilateral and checks if flipping the edge would improve the triangulation.
 * 
 * If the edge should be flipped, it performs the flip.
 * 
 * @param new_elements Vector of triangles created by connecting to the inserted node.
 * @param inserted_node The node that was just added to the mesh.
 * @param original_element The triangle that the inserted node originally fell inside.
 * 
 */
void Mesh::delaunayElementCheck(std::vector<Triangle> & new_elements,
                                Node & inserted_node,
                                Triangle & original_element) {
  // For each proposed new element
  for (size_t i = 0; i < new_elements.size(); ++i) {
    // Find the opposite edge from the node we are inserting
    std::vector<int> edge = findOppositeEdge(new_elements[i], inserted_node);

    // Now find the elements in the current mesh that contain the that edge
    for (std::vector<Triangle>::iterator it_neighbor = elements.begin();
         it_neighbor != elements.end();
         ++it_neighbor) {
      // We are not interested in the element that we are currently working with
      if (it_neighbor->get_id() == original_element.get_id()) {
        continue;
      }
      // If we find another element that has the shared edge, this is the neighboring element
      if (it_neighbor->hasEdge(edge[0], edge[1])) {
        // Grab the quadralaterial formed
        Node * a = getNodeById(edge[0]);
        Node * b = getNodeById(edge[1]);
        Node * c = new_elements[i].getOppositeVertex(*a, *b);
        Node * d = it_neighbor->getOppositeVertex(*a, *b);

        // Check if we need to swap the element
        if (shouldFlipEdge(a, b, c, d)) {
          flipEdge(&new_elements[i], &(*it_neighbor), a, b, c, d);
        }
      }
    }
  }
}

/**
 * @brief Builds the mesh using either simple or Delaunay triangulation.
 * 
 * Starts by manually creating initial triangles from the first 5 nodes.
 * Then it loops through the rest of the nodes and inserts them one by one, 
 * splitting existing elements that contain the new node.
 * 
 * If Delaunay triangulation is requested, it checks whether any edges should 
 * be flipped to satisfy the Delaunay condition and performs those flips before 
 * updating the mesh.
 * 
 * @param ifDelauney Whether to apply Delaunay edge flipping during triangulation.
 */
void Mesh::performTriangulation(bool ifDelauney) {
  if (nodes.empty()) {
    std::cerr << "Error: nodes have not yet been created it" << std::endl;
    exit(EXIT_FAILURE);
  }
  // When we just have the corners at nodes we can only create two triangular elements
  if (nodes.size() == 4) {
    std::vector<std::vector<int> > element_node_ids;

    int tri1[] = {0, 1, 3};  // nodes 1, 2, 4
    int tri2[] = {0, 3, 2};  // nodes 1, 4, 3

    element_node_ids.push_back(std::vector<int>(tri1, tri1 + 3));
    element_node_ids.push_back(std::vector<int>(tri2, tri2 + 3));

    for (size_t i = 0; i < element_node_ids.size(); ++i) {
      std::vector<Node *> corners;
      for (size_t j = 0; j < element_node_ids[i].size(); ++j) {
        corners.push_back(&nodes[element_node_ids[i][j]]);
      }
      elements.push_back(Triangle(i + 1, corners));
    }
    return;
  }
  if (nodes.size() < 4) {
    std::cerr << "Error: not enough nodes created to perform simple triangulation"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!elements.empty()) {
    elements.clear();
  }

  // Manually assemply triangulation for the for first 5 nodes
  std::vector<std::vector<int> > element_node_ids;
  int tri1[] = {0, 2, 4};
  int tri2[] = {1, 3, 4};
  int tri3[] = {2, 3, 4};
  element_node_ids.push_back(std::vector<int>(tri1, tri1 + 3));
  element_node_ids.push_back(std::vector<int>(tri2, tri2 + 3));
  element_node_ids.push_back(std::vector<int>(tri3, tri3 + 3));
  for (size_t i = 0; i < element_node_ids.size(); ++i) {
    std::vector<Node *> corners;
    for (size_t j = 0; j < element_node_ids[i].size(); ++j) {
      corners.push_back(&nodes[element_node_ids[i][j]]);
    }
    elements.push_back(Triangle(i + 1, corners));
  }

  // Begin assempbling remaining nodes
  int elem_id_count = 4;
  for (std::vector<Node>::iterator it_node = nodes.begin() + 5; it_node != nodes.end();
       ++it_node) {
    // For every element in the mesh, check if the node we are inserting falls inside of it, if so store it
    std::vector<Triangle *> interestElements;
    for (std::vector<Triangle>::iterator it_element = elements.begin();
         it_element != elements.end();
         ++it_element) {
      if (it_element->isInside(*it_node)) {
        interestElements.push_back(&*it_element);
      }
    }

    // For every element that contains the node we are inserting
    for (size_t i = 0; i < interestElements.size(); ++i) {
      // Split the element with the node we are inserting, and obtain resulting elementss
      std::vector<Triangle> new_elements =
          interestElements[i]->splitElement(*it_node, elem_id_count);

      // If delauney trinagulation is requested, check if edges must be swapped
      if (ifDelauney) {
        delaunayElementCheck(new_elements, *it_node, *interestElements[i]);
      }

      // Replace existing element with the one with its same id
      replaceExistingElement(*interestElements[i], new_elements);
      // Add remaining new elements to mesh
      for (size_t i = 0; i < new_elements.size(); ++i) {
        elements.push_back(new_elements[i]);
      }
    }
  }
}

/**
 * @brief Performs a basic triangulation of the mesh.
 * 
 * 
 */
void Mesh::simpleTriangulation() {
  performTriangulation(false);
}

/**
 * @brief Performs a Delaunay triangulation of the mesh.
 * 
 * 
 */
void Mesh::delaunayTriangulation() {
  performTriangulation(true);
}

/**
 * @brief Prints all elements in the mesh.
 * 
 * Outputs the element ids and the ids of their three Node vertices.
 * 
 */
void Mesh::printElements() {
  if (elements.empty()) {
    std::cerr << "Error: geometry has not yet been meshed with a meshing algorithm"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << "$elements" << std::endl;
  for (size_t i = 0; i < elements.size(); i++) {
    Triangle & element = elements[i];
    std::cout << element.get_id() << " " << element.get_v1()->get_id() << " "
              << element.get_v2()->get_id() << " " << element.get_v3()->get_id()
              << std::endl;
  }
}
