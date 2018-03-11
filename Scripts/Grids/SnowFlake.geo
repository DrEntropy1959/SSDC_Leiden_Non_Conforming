// Define parameters
charLength = 1.0;
len        = 0.5;

// allow changing transfinite constraints for Blossom
Mesh.FlexibleTransfinite  = 1;
Mesh.SubdivisionAlgorithm = 1;

Geometry.Tolerance = 1e-4 ;
Mesh.Algorithm3D=6;

p01 = newp; Point(01) = {0, 0, 0, charLength};
p02 = newp; Point(02) = {0, -0.5, 0, charLength};
p03 = newp; Point(03) = {0, 0.5, 0, charLength};
p04 = newp; Point(04) = {-0.5, -0.3, 0, charLength};
p05 = newp; Point(05) = {0.5, 0.3, 0, charLength};
p06 = newp; Point(06) = {-0.5, 0.3, 0, charLength};
p07 = newp; Point(07) = {0.5, -0.3, 0, charLength};
p08 = newp; Point(08) = {-1.0, 0, 0, charLength};
p09 = newp; Point(09) = {+1.0, 0, 0, charLength};
p10 = newp; Point(10) = {-0.5, -0.9, 0, charLength};
p11 = newp; Point(11) = {-0.5, 0.9, 0, charLength};
p12 = newp; Point(12) = {0.5, -0.9, 0, charLength};
p13 = newp; Point(13) = {0.5, 0.9, 0, charLength};

l01 = newl; Line(l01) = {p01,p05};
l02 = newl; Line(l02) = {p03,p13};
l03 = newl; Line(l03) = {p05,p09};
l04 = newl; Line(l04) = {p07,p01};
l05 = newl; Line(l05) = {p01,p04};
l06 = newl; Line(l06) = {p02,p12};
l07 = newl; Line(l07) = {p01,p06};
l08 = newl; Line(l08) = {p06,p08};
l09 = newl; Line(l09) = {p08,p04};
l10 = newl; Line(l10) = {p03,p01};
l11 = newl; Line(l11) = {p11,p03};
l12 = newl; Line(l12) = {p06,p11};
l13 = newl; Line(l13) = {p13,p05};
l14 = newl; Line(l14) = {p09,p07};
l15 = newl; Line(l15) = {p10,p02};
l16 = newl; Line(l16) = {p04,p10};
l17 = newl; Line(l17) = {p07,p12};
l18 = newl; Line(l18) = {p02,p01};

ll01 = newl ; Line Loop(ll01) = {l01,-l13,-l02, l10}; s01 = news; Plane Surface(s01) = {ll01};
ll02 = newl ; Line Loop(ll02) = {l14, l04, l01, l03}; s02 = news; Plane Surface(s02) = {ll02};
ll03 = newl ; Line Loop(ll03) = {l06,-l17, l04,-l18}; s03 = news; Plane Surface(s03) = {ll03};
ll04 = newl ; Line Loop(ll04) = {l15, l18, l05, l16}; s04 = news; Plane Surface(s04) = {ll04};
ll05 = newl ; Line Loop(ll05) = {l09,-l05, l07, l08}; s05 = news; Plane Surface(s05) = {ll05};
ll06 = newl ; Line Loop(ll06) = {l07, l12, l11, l10}; s06 = news; Plane Surface(s06) = {ll06};

Transfinite Line{l01,l02,l03,l04,l05,l06,l07,l08,l09,l10,l11,l12,l13,l14,l15,l16,l17,l18}   =  3 Using Progression 1.0;

Transfinite Surface{s01} = {p01,p05,p13,p03};
Transfinite Surface{s02} = {p05,p09,p07,p01};
Transfinite Surface{s03} = {p01,p02,p07,p12};
Transfinite Surface{s04} = {p10,p02,p01,p04};
Transfinite Surface{s05} = {p08,p04,p01,p06};
Transfinite Surface{s06} = {p06,p01,p03,p11};
Recombine   Surface "*" ;

// Make a list of all surfaces which need to be extruded
allParts[] = {s01,s02,s03,s04,s05,s06};

// Extrude geometry and quadrilateral mesh in Z direction //
zdir[] = Extrude{0, 0, -len} { Surface{allParts[]}; Layers{1}; Recombine;};

Physical Surface(30) = {157, 47, 153, 43, 139, 73, 61, 127, 83, 87, 105, 117};
Physical Surface(31) = {52, 162, 74, 140, 96, 118};
Physical Surface(32) = {30, 20, 22, 28, 26, 24};
Physical Volume(40) = {1, 6, 2, 5, 3, 4};
