# Mesh Generator

For this mini-project option, you will write a program to generate a mesh, such as could be used to solve a finite element problem, given a description of a two-dimensional surface.

There are four steps for this assignment (each worth 20 points of the functionality score):

Step 1: Read in mesh area from a file and generate nodes.
Step 2: Perform a simple triangulation.
Step 3: Perform a Delaunay triangulation.
Step 4: Improve the mesh points by adding randomness.

Permitted domain knowledge resources (in addition to the other permitted resources like AoP and cplusplus.com/reference):

 - A paper on Delaunay triangulation algorithms by Lee & Schachter: <https://link.springer.com/article/10.1007%2FBF00977785>

## Step 1

This step involves designing your graph structure and reading in a description of surface to generate a mesh for.

This step should create the program `mesh-step1`, which takes one command line argument specifying the name of a file to read a surface description from. For this project, we will assume the surface is a rectangle for simplicity. 

This input file has the following format:

```
1.0
1.35
0.0 0.0
3.0 0.0
0.0 2.5
3.0 2.5
```

where the first line gives the maximum x-element size, the second gives the maximum y-element size, and the remaining lines each have a vertex number, x-coordinate, and y-coordinate, which together describe the rectangular surface. 

Eventually, your program will select mesh points within this area, and triangulate them to choose triangular elements that represent a mesh of the surface. It will output a file that lists all of the chosen nodes and elements. You will begin by choosing the nodes in a simple way. In subsequent steps, you will add triangulation and sophistication to the node selection.

For a first pass at choosing the mesh points, determine the element size that gives points at regular intervals on the rectangle that do not exceed the maximum element size. For the above example, the x-size would be 1.0 and the y-size would be 1.25. Your program would select and print the resulting nodes (after a "$nodes" header):

  $nodes
  1  0.0 0.0
  2  0.0 2.5
  3  3.0 0.0
  4  3.0 2.5
  5  0.0 1.25
  6  1.0 0.0
  7  1.0 1.25
  8  1.0 2.5
  9  2.0 0.0
  10 2.0 1.25
  11 2.0 2.5
  12 3.0 1.25

where each node has a node ID, followed by its x- and y-coordinate.

You are free to implement the graph structure however you like, but here are some suggestions.

  - A Node class
  - An Mesh class

A rectangular mesh has four corners. Once you generate the nodes, it has nodes that are not connected yet. In subsequent steps, you will add edges and triangular elements.

We recommend writing your code in a re-usable way for later parts, i.e. you should write as little code as possible in the source file with this step's main.

Once you have thoroughly tested this step, add/commit/push before continuing to the next step.

## Step 2

In this step, write the program `mesh-step2`, which performs a simple triangulation and prints a description of the mesh.

Next, you need to add each of the nodes to the mesh and do the triangulation. For this step, the algorithm is as follows: 
  
 1. From all of the mesh points, remove the four vertices of the rectangle, and add them to the mesh, connecting them with edges that form a rectangle.
 2. Order the remaining mesh points, first by x-coordinate, then by y-coordinate if two nodes have the same x-coordinate.
 3. Remove the next mesh point, and add it to the triangulation, connecting it to each of the rectangle's vertices. Add each triangle that this creates to the triangulation.
 4. For each of the remaining mesh points,
     - Add the point to the triangulation.
     - Connect the point to the vertices of its enclosing triangle. (If the point is on an external edge, connect it to the nodes on either side of that edge.)
     - Add each triangle that this creates to the triangulation.

Note that this will create a triangulation, but not one optimized to have evenly sized elements or triangles without very small angles. You will improve the triangulation in a future step. 
  
For this step, your program should output a file with this node and element information in this format--for the nodes, give the node ID, x-coordinate, and y-coordinate; for the elements, give the element ID, and node IDs of the vertices of each element: 

  $nodes
  1  0.0  0.0
  2  0.0  2.5
  3  3.0  0.0
  4  3.0  2.5
  5  0.0  1.25
  6  1.0  0.0
  7  1.0  1.25
  8  1.0  2.5
  9  2.0  0.0
  10 2.0  1.25
  11 2.0  2.5
  12 3.0  1.25
  $elements
  1  1  5  6
  2  3  5  8
  3  5  8  11
  4  5  6  9
  5  5  11 4
  6  5  9  2
  7  5  7  4
  8  5  7  2
  9  7  4  10
  10 7  10 2
  11 10 4  12
  12 10 12 2

The number of decimal places displayed is up to you, as long as ID values are decimal integers, and the numbers are delimited by spaces.
  
To visualize your output and help you check your work, you should also plot the rectangle along with each mesh point you have chosen. We provide a simple plotting script.

Again, you are free to implement this in any way you prefer. We suggest:

  - A Triangle class with a vector of nodes or a vector of edges, which can determine if a given point is in its interior and calculate its minimum angle (for the next part)
  - A TriangulationStrategy abstract class with concrete subclasses for a simple or Delaunay triangulation.

You should write as little code as possible in the source file with this step's main, writing most of it in your classes or other source files for re-usability.

Once you have thoroughly tested this step, make sure Step 1 still works well and add/commit/push before continuing to the next step.

## Step 3

For this step, you will improve the triangulation algorithm by doing a Delaunay triangulation, which ensures no point in the interior of a triangular element is within the circumcircle of a different element. This also has the effect of maximizing the minimum interior angles of the elements. Write the program `mesh-step3`, which implements this algorithm.

The new algorithm is similar to the triangulation you implemented in Step 2, except that after each mesh point is added, you will check to see if swapping edges improves the result. Add to the algorithm (at the end of items 3 and 4): 
  
       - For each new triangle you have just formed, consider the
         quadrilateral formed by it and its neighboring triangle (if
         one exists). For each quadrilateral, there are two possible
         diagonals. If swapping the existing diagonal with the other
         possibility results in a larger minimum interior angle of
         the two triangles that comprise it, swap the diagonal edge.

A reference for this (and another) algorithm are given here: <https://link.springer.com/article/10.1007%2FBF00977785>

## Step 4

There are many ways to select mesh points for a triangulation, but evenly spaced points may not be optimal in all applications. For this step, you will add randomness to the way mesh points are generated. Write this step in `mesh-step4`.

You may do this in a way of your choosing, while preserving the number of mesh points indicated in Step 1. Two possibilities (not exhaustive) are:

 1. Randomly generate twice as many mesh points as you need in the area of the surface. Then, remove half of them that are closest to other points. 

 2. Generate points as in Step 1, but then for each evenly spaced mesh point, move it by some radius in a random direction.

See if you can implement this algorithm with relatively small changes to your existing code.

Once you have thoroughly tested this step, make sure Steps 1, 2, and 3 still work well and add/commit/push before running grade.

Review the overall README and make sure your `TESTING.txt` is polished and you have considered all of the elements of code quality.
