// Define parameters
 charLength = 1.0;

// allow changing transfinite constraints for Blossom
Mesh.FlexibleTransfinite = 1;
Mesh.RecombineAll        = 1;
Mesh.Algorithm           = 8;
// Subdivides triangles in 3 quads
Mesh.SubdivisionAlgorithm= 1;

Lay  = 2 ;
Nval = 1 ;
fac =  1.000     ; //  Scaling factor
dy  =  fac*0.1 ;
dz  =  dy      ;
p   = 9 ;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//                               INLET
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Begin Points Inlet 
p001   = newp; Point(p001) = { 0.0 , 0.0 ,  0.0,charLength};
p002   = newp; Point(p002) = { 1.0 , 0.0 ,  0.0,charLength};
p003   = newp; Point(p003) = { 1.0 , dy  ,  0.0,charLength};
p004   = newp; Point(p004) = { 0.0 , dy  ,  0.0,charLength};

// Vertical lines
l001  = newl; Line(l001) = {p001,p002};
l002  = newl; Line(l002) = {p002,p003};
l003  = newl; Line(l003) = {p003,p004};
l004  = newl; Line(l004) = {p004,p001};

// Line Loops
ll001 = newl; Line Loop(ll001) = {l001, l002, l003, l004};

// Inlet  Bottom, middle and Top walls
Transfinite Line{l001,l003} = 1 + 2^p Using Progression 1.0;
Transfinite Line{l002,l004} = 0       Using Progression 1.0;

s001 = news; Plane Surface(s001) = {ll001};

Transfinite Surface{s001} = {p001,p002,p003,p004};

Recombine Surface{s001};

allParts[] = {s001} ;
zdir[] = Extrude{0, 0, dz} { Surface{allParts[]}; Layers{Lay}; Recombine;};

//Transfinite Surface "*" ;
Recombine   Surface "*" ;
//Transfinite Volume  "*" ;
Recombine   Volume  "*" ;

Physical Surface(29) = {27, 19};
Physical Surface(41) = {28};
Physical Surface(42) = {6};
Physical Surface(43) = {23};
Physical Surface(44) = {15};
Physical Volume(30) = {1};
