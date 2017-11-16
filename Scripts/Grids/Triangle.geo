// Define parameters
 charLength = 1.0;

// allow changing transfinite constraints for Blossom
Mesh.FlexibleTransfinite = 1;
Mesh.RecombineAll        = 1;
Mesh.Algorithm           = 8;
// Subdivides triangles in 3 quads
Mesh.SubdivisionAlgorithm= 1;

Lay  = 1 ;
Nval = 1 ;
fac =  1.000     ; //  Scaling factor
far  = 3.0*fac ;
near = 1.0*fac ;
aspc = 0.1*fac ;
Sr3 = Sqrt[3] ;

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//                               INLET
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Begin Points Inlet 
p001   = newp; Point(p001) = {-1.0 ,-1.0 ,  0.0,charLength};
p002   = newp; Point(p002) = {+0.0 ,-1.0 ,  0.0,charLength};
p003   = newp; Point(p003) = {+1.0 ,-1.0 ,  0.0,charLength};
p004   = newp; Point(p004) = {+0.5 ,Sqrt[3]/2-1,  0.0,charLength};
p005   = newp; Point(p005) = {+0.0 ,Sqrt[3]/1-1,  0.0,charLength};
p006   = newp; Point(p006) = {-0.5 ,Sqrt[3]/2-1,  0.0,charLength};

p007   = newp; Point(p007) = {+0.0 ,-1.0+1/Sqrt[3] ,  0.0,charLength};

// End Points Inlet 

// Vertical lines
l001  = newl; Line(l001) = {p001,p002};
l002  = newl; Line(l002) = {p002,p003};
l003  = newl; Line(l003) = {p003,p004};
l004  = newl; Line(l004) = {p004,p005};
l005  = newl; Line(l005) = {p005,p006};
l006  = newl; Line(l006) = {p006,p001};
l007  = newl; Line(l007) = {p002,p007};
l008  = newl; Line(l008) = {p007,p004};
l009  = newl; Line(l009) = {p007,p006};

// Line Loops
ll001 = newl; Line Loop(ll001) = {l001, l007, l009, l006};
ll002 = newl; Line Loop(ll002) = {l002, l003,-l008,-l007};
ll003 = newl; Line Loop(ll003) = {l004, l005,-l009, l008};

Transfinite Line{l001,l002,l003,l004,l005,l006,l007,l008,l009} = 1 Using Progression 1.0;
//

s001 = news; Plane Surface(s001) = {ll001};
s002 = news; Plane Surface(s002) = {ll002};
s003 = news; Plane Surface(s003) = {ll003};
Transfinite Surface{s001} = {p001,p002,p007,p006};
Transfinite Surface{s002} = {p002,p003,p004,p007};
Transfinite Surface{s003} = {p004,p005,p006,p007};

Recombine Surface{s001,s002,s003};

allParts[] = {s001:s003} ;
zdir[] = Extrude{0, 0, 1} { Surface{allParts[]}; Layers{Lay}; Recombine;};

//Transfinite Surface "*" ;
Recombine   Surface "*" ;
//Transfinite Volume  "*" ;
Recombine   Volume  "*" ;
Physical Surface(21) = {46, 50, 24, 68, 36, 72};
Physical Surface(22) = {15, 14, 13};
Physical Surface(23) = {81, 59, 37};
Physical Volume(31) = {3, 2, 1};
