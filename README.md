Rasterization Basics
===============

### Using Eigen

ITo simplify the code, I used [Eigen](http://eigen.tuxfamily.org/).

Ex.1: Load and Render a 3D model
------------------

The provided code rasterizes a correct silhouette of the object rendered. The cameras take into account the size of the framebuffer, properly adapting the aspect ratio to not distort the image whenever the framebuffer is resized. 


Ex.2: Wireframe Render
----------------

Render of the bunny using a wireframe;

Ex.3: Flat Shading
-------------

Each triangle is rendered using a unique normal.

Ex.4: Per-Vertex Shading
-------------

The normals are specified on the vertices of the mesh, the color is computed for each vertex, and then interpolated in the interior of the triangle.


Ex.5: Object Transformation
----------------------------------

Adds an option to rotate the object around the y-axis centered on its barycenter.