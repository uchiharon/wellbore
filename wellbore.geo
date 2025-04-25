// Gmsh script to create a 2D mesh around a wellbore with Quad9 elements
// Parameters
well_radius = 0.5;        // Radius of the wellbore (m)
domain_size = 10.0;       // Size of the square domain (m)
lc_well = 0.2;           // Mesh size near wellbore
lc_domain = 2.0;          // Mesh size at domain boundary

// Points for wellbore
Point(1) = {0, 0, 0, lc_well};              // Center of wellbore
Point(2) = {well_radius, 0, 0, lc_well};    // Wellbore boundary (right)
Point(3) = {0, well_radius, 0, lc_well};    // Wellbore boundary (top)
Point(4) = {-well_radius, 0, 0, lc_well};   // Wellbore boundary (left)
Point(5) = {0, -well_radius, 0, lc_well};   // Wellbore boundary (bottom)

// Points for outer domain
Point(6) = {-domain_size/2, -domain_size/2, 0, lc_domain}; // Bottom-left
Point(7) = {domain_size/2, -domain_size/2, 0, lc_domain};  // Bottom-right
Point(8) = {domain_size/2, domain_size/2, 0, lc_domain};   // Top-right
Point(9) = {-domain_size/2, domain_size/2, 0, lc_domain};  // Top-left

// Lines for wellbore (circle approximated by 4 arcs)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Lines for outer domain
Line(5) = {6, 7}; // Bottom edge
Line(6) = {7, 8}; // Right edge
Line(7) = {8, 9}; // Top edge
Line(8) = {9, 6}; // Left edge

// Line loops
Line Loop(1) = {1, 2, 3, 4};           // Wellbore boundary
Line Loop(2) = {5, 6, 7, 8};           // Outer domain boundary

// Surface
Plane Surface(1) = {2, 1};             // Domain with wellbore as hole

// Physical groups for boundaries
Physical Curve("wellbore", 1) = {1, 2, 3, 4}; // Wellbore boundary
Physical Curve("top_bottom", 2) = {5, 7}; // Boundary
Physical Surface("domain", 1) = {1};          // Domain surface

// Mesh settings
Mesh.Algorithm = 8;        // Frontal-Delaunay for quads
Mesh.ElementOrder = 2;     // Second-order elements (for Quad9)
Mesh.RecombineAll = 1;     // Recombine triangles into quadrilaterals
Mesh.SubdivisionAlgorithm = 1; // Ensure proper recombination
Mesh 2;                    // Generate 2D mesh
