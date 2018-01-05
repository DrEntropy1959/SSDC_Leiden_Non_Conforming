
// Define parameters
charLength = 1.0;
fac = 2/3 ;
//fac = 1/5 ;

// allow changing transfinite constraints for Blossom
Mesh.FlexibleTransfinite = 1;

// Grid boundary layer size = d2-d1
d1 = 0.5;
d2 = 1.0;

// Wake
pCenter = 1.0;

dwS1 =  3.2;
dwS2 = 10.0;
dwS3 = -0.5;
dwS4 =  1.0;

//  Mach 1.5  => shockx = 0.0  
shockx = 0.0 ;
dw1 = 2.15;
dw2 = 8.0;
dw3 = 14.0;
dw4 = 2.15;

dw5 = 16.0;
dw6 = 09.0;

dw7 =  1.0;


// Far-field
xfff = 20.0;
yfff = 25.0;
xffb = 50.0;
yffb = 25.0;

// Use to define 4 corners of the domain
xMax = 20.0;
xMin = -8.0;
yMax = +8.0;

// Points square cylinder
p1 = newp; Point(p1) = {-d1,+d1,0.0,charLength};
p2 = newp; Point(p2) = {-d1,-d1,0.0,charLength};
p3 = newp; Point(p3) = {+d1,-d1,0.0,charLength};
p4 = newp; Point(p4) = {+d1,+d1,0.0,charLength};

// Points BL
p5 = newp; Point(p5) = {-d2,+d2,0.0,charLength};
p6 = newp; Point(p6) = {-d2,-d2,0.0,charLength};
p7 = newp; Point(p7) = {+d2,-d2,0.0,charLength};
p8 = newp; Point(p8) = {+d2,+d2,0.0,charLength};

// Lines square cylinder
l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p1,p4};

// Lines BL
l5  = newl; Line(l5 ) = {p1,p5};
l6  = newl; Line(l6 ) = {p2,p6};
l7  = newl; Line(l7 ) = {p3,p7};
l8  = newl; Line(l8 ) = {p4,p8};
l9  = newl; Line(l9 ) = {p5,p6};
l10 = newl; Line(l10) = {p6,p7};
l11 = newl; Line(l11) = {p7,p8};
l12 = newl; Line(l12) = {p5,p8};
 

// Points wake
p9  = newp; Point(p9)  = {+dw7,+dw4,0.0,charLength};
p10 = newp; Point(p10) = {+dw7,-dw4,0.0,charLength};
p11 = newp; Point(p11) = {+dw3,-dw2,0.0,charLength};
p12 = newp; Point(p12) = {+dw3,+dw2,0.0,charLength};

p111 = newp; Point(p111) = {+dw5,-dw6,0.0,charLength};
p112 = newp; Point(p112) = {+dw5,+dw6,0.0,charLength};

p113 = newp; Point(p113) = {+dw3,-d2,0.0,charLength};
p114 = newp; Point(p114) = {+dw3,+d2,0.0,charLength};

// Lines wake
l14 = newl; Line(l14) = {p10,p11};
l16 = newl; Line(l16) = {p9,p12};

l114 = newl; Line(l114) = {p11,p111};
l115 = newl; Line(l115) = {p12,p112};
l116 = newl; Line(l116) = {p111,p112};

l117 = newl; Line(l117) = {p11,p113};
l118 = newl; Line(l118) = {p113,p114};
l119 = newl; Line(l119) = {p12,p114};
l120 = newl; Line(l120) = {p7,p113};
l121 = newl; Line(l121) = {p8,p114};

// Connect wake to BL
l17 = newl; Line(l17) = {p8,p9};
l18 = newl; Line(l18) = {p7,p10};

// Far-field points
p13 = newp; Point(p13) = {-xfff,yfff,0.0,charLength};
p14 = newp; Point(p14) = {-xfff,-yfff,0.0,charLength};
p15 = newp; Point(p15) = {xffb,-yffb,0.0,charLength};
p16 = newp; Point(p16) = {xffb,yffb,0.0,charLength};

// Far-field lines
l19 = newl; Line(l19) = {p13,p14};
l20 = newl; Line(l20) = {p14,p15};
l21 = newl; Line(l21) = {p15,p16};
l22 = newl; Line(l22) = {p13,p16};

// Points second butterfly layer
p17 = newp; Point(p17) = {-dw1+shockx,+dw1,0.0,charLength};
p18 = newp; Point(p18) = {-dw1+shockx,-dw1,0.0,charLength};


// Lines second butterfly layer
l23 = newl; Line(l23) = {p9,p17};

l24 = newl; Line(l24) = {p17,p18};

l25 = newl; Line(l25) = {p10,p18};
l26 = newl; Line(l26) = {p5,p17};
l27 = newl; Line(l27) = {p6,p18};

pCenter = newp; Point(pCenter) = {6.5,0.0,0.0,charLength};
p19 = newp; Point(p19) = {-dwS1+shockx,-dwS1,0.0,charLength};
p20 = newp; Point(p20) = {-dwS1+shockx,+dwS1,0.0,charLength};

l28 = newl; Circle(l28) = {p19,pCenter,p20};
l29 = newl; Line(l29) = {p18,p19};
l30 = newl; Line(l30) = {p17,p20};

p21 = newp; Point(p21) = {+dwS3+shockx,-dwS2,0.0,charLength};
p22 = newp; Point(p22) = {+dwS3+shockx,+dwS2,0.0,charLength};
p23 = newp; Point(p23) = {+dwS4+shockx,-dwS2,0.0,charLength};
p24 = newp; Point(p24) = {+dwS4+shockx,+dwS2,0.0,charLength};

l31 = newl; Line(l31) = {p20,p22};
l32 = newl; Line(l32) = {p22,p24};
l33 = newl; Line(l33) = {p17,p24};

l34 = newl; Line(l34) = {p19,p21};
l35 = newl; Line(l35) = {p21,p23};
l36 = newl; Line(l36) = {p18,p23};

// Line loops
ll1 = newl; Line Loop(ll1) = {-l1,l5,l9,-l6};
ll2 = newl; Line Loop(ll2) = {-l2,l6,l10,-l7};
ll3 = newl; Line Loop(ll3) = {-l3,l7,l11,-l8};
ll4 = newl; Line Loop(ll4) = { l4,l8,-l12,-l5};
ll7 = newl; Line Loop(ll7) = {l19,l20,l21,-l22,-l28,l34,l35,-l36,-l25,l14,l114,l116,-l115,-l16,l23,l33,-l32,-l31};
ll8 = newl; Line Loop(ll8) = {l23,-l26,l12,l17};
ll9 = newl; Line Loop(ll9) = {l26,l24,-l27,-l9};
ll10 = newl; Line Loop(ll10) = {l27,-l25,-l18,-l10};
ll116 = newl; Line Loop(ll116) = {l119,-l118,-l117,l114,l116,-l115};

ll117 = newl; Line Loop(ll117) = {-l28,-l29,-l24,l30};
ll118 = newl; Line Loop(ll118) = {-l30,l33,-l32,-l31};
ll119 = newl; Line Loop(ll119) = { l29,l34,l35,-l36};

ll120 = newl; Line Loop(ll120) = { l14,l117,-l120,l18};
ll121 = newl; Line Loop(ll121) = { l120,l118,-l121,-l11};
ll122 = newl; Line Loop(ll122) = { l121,-l119,-l16,-l17};

// Surfaces
s1 = news; Plane Surface(s1) = {ll1};
s2 = news; Plane Surface(s2) = {ll2};
s3 = news; Plane Surface(s3) = {ll3};
s4 = news; Plane Surface(s4) = {ll4};
s7 = news; Plane Surface(s7) = {ll7};
s8 = news; Plane Surface(s8) = {ll8};
s9 = news; Plane Surface(s9) = {ll9};
s10 = news; Plane Surface(s10) = {ll10};
s11 = news; Plane Surface(s11) = {ll116};

s12 = news; Plane Surface(s12) = {ll117};
s13 = news; Plane Surface(s13) = {ll118};
s14 = news; Plane Surface(s14) = {ll119};

s17 = news; Plane Surface(s17) = {ll120};
s18 = news; Plane Surface(s18) = {ll121};
s19 = news; Plane Surface(s19) = {ll122};

Transfinite Surface{s18} = {p7,p113,p114,p8};

Transfinite Surface{s1} = {p1,p5,p6,p2};
Transfinite Surface{s2} = {p2,p6,p7,p3};
//Transfinite Surface{s3} = {p3,p7,p8,p4};
//Transfinite Surface{s9} = {p17,p18,p6,p5};
Transfinite Surface{s4} = {p4,p8,p5,p1};
Transfinite Surface{s12} = {p19,p18,p17,p20};

Transfinite Surface{s13} = {p20,p17,p24,p22};
Transfinite Surface{s14} = {p19,p21,p23,p18};

//// DISTRIBUTION OF THE POINTS: The game has to be played here ;-)
///////////////////////////////////////////////////////////////////
// North, South, East square walls
Transfinite Line{l3}          = 35*fac Using Bump 0.6;
Transfinite Line{l2,-l4}      = 25*fac Using Bump 0.5;
Transfinite Line{l10,l12}     = 25*fac Using Progression 0.97;
Transfinite Line{l11,l118}    = 30*fac Using Bump 0.90;

// West square wall
Transfinite Line{l1,l9}           = 20*fac Using Bump 0.5;

// Four lines which connect the square and the end of the BL
Transfinite Line{-l5,-l6,-l7,-l8} = 18*fac Using Progression 0.90;

// First divergent part of the wake: sides
Transfinite Line{-l17,-l18} = 13*fac Using Progression 0.92;

// End first and second divergent part of the wake
Transfinite Line{l116} = 40*fac Using Bump 3.0;

// End third Straight  part of the wake
Transfinite Line{l114,l115} = 10*fac Using Progression 1.0;

// Second divergent part of wake: sides
//Transfinite Line{-l14,-l16} = 70 Using Progression 0.985;
Transfinite Line{-l14,-l16} = 50*fac Using Progression 0.95;


// N and S 3rd horizontal connectors
Transfinite Line{-l23,-l25} = 24*fac Using Progression 0.97;
// NW and SW diagonal connectors
Transfinite Line{-l26,-l27} = 17*fac Using Progression 0.95;

//Transfinite Line{ l117, l119} =  40 Using Progression 0.97;
Transfinite Line{ l117, l119} =  20*fac Using Progression 0.97;
//Transfinite Line{-l120,-l121} = 170 Using Progression 1.00;
Transfinite Line{-l120,-l121} = 080*fac Using Progression 0.98;

// Inlet
Transfinite Line{l19}        = 10*fac Using Progression 1.0;

// Outlet
Transfinite Line{l21}        = 15*fac Using Progression 1.0;

// Far-field sides
Transfinite Line{-l22,-l20}  = 35*fac Using Progression 1.0;

// Shock
Transfinite Line{-l28, l24}  = 35*fac Using Bump        0.99;
Transfinite Line{ l29, l30}  = 10*fac Using Progression 1.0;
Transfinite Line{-l32,+l35}  = 10*fac Using Progression 1.0;
Transfinite Line{-l31,-l34}  = 20*fac Using Progression 0.92;
Transfinite Line{-l33,-l36}  = 20*fac Using Progression 0.89;


//// Subdivides triangles in 3 quads
//Mesh.SubdivisionAlgorithm=1;

// Recombine surfaces
Recombine Surface{s1};
Recombine Surface{s2};
Recombine Surface{s3};
Recombine Surface{s4};
Recombine Surface{s7};
Recombine Surface{s8};
Recombine Surface{s9};
Recombine Surface{s10};
Recombine Surface{s11};

Recombine Surface{s12};

Recombine Surface{s13};
Recombine Surface{s14};

Recombine Surface{s17};
Recombine Surface{s18};
Recombine Surface{s19};

// Make a list of all surfaces which need to be extruded
allParts[] = {s1,s2,s3,s4,s7,s8,s9,s10,s11,s12,s13,s14,s17,s18,s19};

// Extrude geometry and quadrilateral mesh in Z direction //
//First Grid//zdir[] = Extrude{0, 0, 0.5} { Surface{allParts[]}; Layers{2}; Recombine;};
zdir[] = Extrude{0, 0, 1.0} { Surface{allParts[]}; Layers{2}; Recombine;};

// 3) Set background field
//    --------------------

//Field[5] = Min;
//Field[5].FieldsList = {2,4};
//Background Field = 5;

// Far-field
Physical Surface(483) = {183, 195, 191, 187};
// Solid wall
Physical Surface(484) = {81, 147, 125, 103};
// Front symmetry plane
Physical Surface(485) = {252, 350, 460, 482, 438, 394, 416, 372, 296, 318, 274, 160, 138, 116, 94};
// Back symmetry plane
Physical Surface(486) = {62, 66, 71, 72, 70, 69, 67, 68, 64, 63, 65, 58, 61, 59, 60};
Physical Volume(487) = {5, 9, 14, 15, 13, 12, 11, 10, 7, 8, 2, 1, 4, 3, 6};
