///////////////////////////////////////////////////////////////////////////////
// Construct a cube with edge length equals to Pi and centered at (0,0,0).
// Each face of the cube is a physical face so that different boundary
// conditions can be imposed.
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Paramters
///////////////////////////////////////////////////////////////////////////////

// Characteristic length
char_length = 1.0;

// Half edge length
d = Pi;

// Number of points on the edge (used for meshing)
n_points = 10;


///////////////////////////////////////////////////////////////////////////////
// Geometry
///////////////////////////////////////////////////////////////////////////////

// Points that define the face that will be extruded 
p1 = newp; Point(p1) = {-d,-d,d,char_length};
p2 = newp; Point(p2) = {+d,-d,d,char_length};
p3 = newp; Point(p3) = {+d,+d,d,char_length};
p4 = newp; Point(p4) = {-d,+d,d,char_length};

// Lines
l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

// Line loops
ll1 = newl; Line Loop(ll1) = {l1,l2,l3,l4};

// Surfaces
s1 = news; Plane Surface(s1) = {ll1};

//Tranfinite surfaces
Transfinite Surface{s1} = {p1,p2,p3,p4};

// Transfinite lines
Transfinite Line{l1,l2,l3,l4} = n_points;

// Recombine surfaces
Recombine Surface{s1};

// Make a list of all surfaces which need to be extruded
allParts[] = {s1};

// Extrude geometry and quadrilateral mesh in Z direction //
zdir[] = Extrude{0, 0, -2*d} { Surface{allParts[]}; Layers{n_points-1}; Recombine;};


///////////////////////////////////////////////////////////////////////////////
// Boundary conditions and interior domain
///////////////////////////////////////////////////////////////////////////////

// x-y plane positive z
Physical Surface(40) = {6};

// x-y plane negative z
Physical Surface(41) = {28};

// y-z plane positive x
Physical Surface(42) = {19};

// y-z plane negative x
Physical Surface(43) = {27};

// x-z positive y
Physical Surface(44) = {23};

// x-z negative y
Physical Surface(45) = {15};

// Interior domain
Physical Volume(46) = {1};
