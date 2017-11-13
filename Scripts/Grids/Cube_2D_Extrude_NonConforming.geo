// Define parameters
charLength = 1.0;

// allow changing transfinite constraints for Blossom
Mesh.FlexibleTransfinite  = 1;
Mesh.SubdivisionAlgorithm = 1;

Geometry.Tolerance = 1e-4 ;
Mesh.Algorithm3D=6;

// Edge length square
len= 0.50 ;
ep = 0.25 ;
pH = 0.99 ;
pM = 1.00 ;
pV = 0.99 ;

// nn: 1,2,3,4,6,8,12,16,24,32
nn  = 2 ;
nH  = nn+1 ;  //  Horizontal connector spacing
nM  = nn+3 ;  //  Middle     connector spacing
nV  = nn+1 ;  //  Vertical   connector spacing

p01  = newp; Point(p01) = {-1.0   , -1.0   ,  0.0,charLength};
p02  = newp; Point(p02) = { 0.0   , -1.0-ep,  0.0,charLength};
p03  = newp; Point(p03) = {+1.0   , -1.0   ,  0.0,charLength};
p04  = newp; Point(p04) = {-1.0-ep,  0.0   ,  0.0,charLength};
p05  = newp; Point(p05) = { 0.0   ,  0.0   ,  0.0,charLength};
p06  = newp; Point(p06) = {+1.0+ep,  0.0   ,  0.0,charLength};
p07  = newp; Point(p07) = {-1.0   , +1.0   ,  0.0,charLength};
p08  = newp; Point(p08) = { 0.0   , +1.0+ep,  0.0,charLength};
p09  = newp; Point(p09) = {+1.0   , +1.0   ,  0.0,charLength};

// Horizontal lines:   moving Left -> Right
lH01 = newl; Line(lH01) = {p01,p02};
lH02 = newl; Line(lH02) = {p02,p03};
lH03 = newl; Line(lH03) = {p04,p05};
lH04 = newl; Line(lH04) = {p05,p06};
lH05 = newl; Line(lH05) = {p07,p08};
lH06 = newl; Line(lH06) = {p08,p09};

// Vertical lines:   moving Bottom -> Top
lV01 = newl; Line(lV01) = {p01,p04};
lV02 = newl; Line(lV02) = {p04,p07};
lV03 = newl; Line(lV03) = {p02,p05};
lV04 = newl; Line(lV04) = {p05,p08};
lV05 = newl; Line(lV05) = {p03,p06};
lV06 = newl; Line(lV06) = {p06,p09};

Transfinite Line{lH01,-lH02,lH05,-lH06}   =  nH Using Progression pH;
Transfinite Line{lH03,-lH04,lV03,-lV04}   =  nM Using Progression pM;
Transfinite Line{lV01,-lV02,lV05,-lV06}   =  nV Using Progression pV;

ll01 = newl; Line Loop(ll01) = {lH01,lV03,-lH03,-lV01};
ll02 = newl; Line Loop(ll02) = {lH02,lV05,-lH04,-lV03};
ll03 = newl; Line Loop(ll03) = {lH03,lV04,-lH05,-lV02};
ll04 = newl; Line Loop(ll04) = {lH04,lV06,-lH06,-lV04};

s01 = news; Plane Surface(s01) = {ll01}; Recombine Surface{s01} ;
s02 = news; Plane Surface(s02) = {ll02}; Recombine Surface{s02} ;
s03 = news; Plane Surface(s03) = {ll03}; Recombine Surface{s03} ;
s04 = news; Plane Surface(s04) = {ll04}; Recombine Surface{s04} ;

//Transfinite Surface{s01} = {p01,p02,p05,p04};
//Transfinite Surface{s02} = {p02,p03,p06,p05};
//Transfinite Surface{s03} = {p04,p05,p08,p07};
//Transfinite Surface{s04} = {p05,p06,p09,p08};

//Surface{s01} = {p01,p02,p05,p04};
//Surface{s02} = {p02,p03,p06,p05};
//Surface{s03} = {p04,p05,p08,p07};
//Surface{s04} = {p05,p06,p09,p08};

// Recombine surfaces
//Transfinite Surface "*" ;
Recombine   Surface "*" ;

Coherence Mesh ;

// Make a list of all surfaces which need to be extruded
allParts[] = {s01,s02,s03,s04};

// Extrude geometry and quadrilateral mesh in Z direction //
zdir[] = Extrude{0, 0, -len} { Surface{allParts[]}; Layers{1}; Recombine;};

//Transfinite Volume  "*" ;
Recombine   Volume  "*" ;
Physical Surface(31) = {29, 41, 51, 55, 103, 99, 85, 81};
Physical Surface(32) = {64, 108, 42, 86};
Physical Surface(33) = {18, 20, 17, 19};
Physical Volume(41) = {4, 2, 3, 1};
