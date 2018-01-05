// parameters
Mesh.FlexibleTransfinite  = 1;
Mesh.SubdivisionAlgorithm = 1;

scale = 4 ;
charLength = 1.0;

xDcR = +0.5;
xDmR = +0.25;
xDfR = +0.17;

xDcL = +0.5;
xDmL = +0.3;
xDfL = +0.07;
xDSh = +0.035;

yDc  = +0.5;
yDm  = +0.25;
yDf  = +0.15;

nXShock         = 07*scale ;

nYCoarse2Medium = 06*scale ;
nYMedium2Finest = 04*scale ;

nXCoarseR       = 07*scale ;
nXCoarseL       = 04*scale ;
nYCoarseR       = 08*scale ;
nYCoarseL       = 06*scale ;

nXMediumR       = 07*scale ;
nXMediumL       = 10*scale ;
nYMediumR       = 15*scale ;
nYMediumL       = 08*scale ;

nXFinestR       = 09*scale ;
nXFinestL       = 03*scale ;
nYFinest        = 14*scale ;

// define Coarse_ perimeter points
pCoarse_NE  = newp; Point(pCoarse_NE)  = {+xDcR,+yDc ,0.0,charLength};
pCoarse_SE  = newp; Point(pCoarse_SE)  = {+xDcR,-yDc ,0.0,charLength};
pCoarse_NW  = newp; Point(pCoarse_NW)  = {-xDcL,+yDc ,0.0,charLength};
pCoarse_SW  = newp; Point(pCoarse_SW)  = {-xDcL,-yDc ,0.0,charLength};

pCoarse_ShockNE  = newp; Point(pCoarse_ShockNE)  = {+xDSh,+yDc ,0.0,charLength};
pCoarse_ShockNW  = newp; Point(pCoarse_ShockNW)  = {-xDSh,+yDc ,0.0,charLength};
pCoarse_ShockSE  = newp; Point(pCoarse_ShockSE)  = {+xDSh,-yDc ,0.0,charLength};
pCoarse_ShockSW  = newp; Point(pCoarse_ShockSW)  = {-xDSh,-yDc ,0.0,charLength};

l01 = newl; Line(l01) = {pCoarse_SW     , pCoarse_ShockSW};
l02 = newl; Line(l02) = {pCoarse_ShockSW, pCoarse_ShockSE};
l03 = newl; Line(l03) = {pCoarse_ShockSE, pCoarse_SE     };
l04 = newl; Line(l04) = {pCoarse_SE     , pCoarse_NE     };
l05 = newl; Line(l05) = {pCoarse_NE     , pCoarse_ShockNE};
l06 = newl; Line(l06) = {pCoarse_ShockNE, pCoarse_ShockNW};
l07 = newl; Line(l07) = {pCoarse_ShockNW, pCoarse_NW     };
l08 = newl; Line(l08) = {pCoarse_NW     , pCoarse_SW     };

// define Medium perimeter points
pMedium_NE  = newp; Point(pMedium_NE)  = {+xDmR,+yDm ,0.0,charLength};
pMedium_SE  = newp; Point(pMedium_SE)  = {+xDmR,-yDm ,0.0,charLength};
pMedium_NW  = newp; Point(pMedium_NW)  = {-xDmL,+yDm ,0.0,charLength};
pMedium_SW  = newp; Point(pMedium_SW)  = {-xDmL,-yDm ,0.0,charLength};

pMedium_ShockNE  = newp; Point(pMedium_ShockNE)  = {+xDSh,+yDm ,0.0,charLength};
pMedium_ShockNW  = newp; Point(pMedium_ShockNW)  = {-xDSh,+yDm ,0.0,charLength};
pMedium_ShockSE  = newp; Point(pMedium_ShockSE)  = {+xDSh,-yDm ,0.0,charLength};
pMedium_ShockSW  = newp; Point(pMedium_ShockSW)  = {-xDSh,-yDm ,0.0,charLength};

l11 = newl; Line(l11) = {pMedium_SW     , pMedium_ShockSW};
l12 = newl; Line(l12) = {pMedium_ShockSW, pMedium_ShockSE};
l13 = newl; Line(l13) = {pMedium_ShockSE, pMedium_SE     };
l14 = newl; Line(l14) = {pMedium_SE     , pMedium_NE     };
l15 = newl; Line(l15) = {pMedium_NE     , pMedium_ShockNE};
l16 = newl; Line(l16) = {pMedium_ShockNE, pMedium_ShockNW};
l17 = newl; Line(l17) = {pMedium_ShockNW, pMedium_NW     };
l18 = newl; Line(l18) = {pMedium_NW     , pMedium_SW     };

// define Finest perimeter points
pFinest_NE  = newp; Point(pFinest_NE)  = {+xDfR,+yDf ,0.0,charLength};
pFinest_SE  = newp; Point(pFinest_SE)  = {+xDfR,-yDf ,0.0,charLength};
pFinest_NW  = newp; Point(pFinest_NW)  = {-xDfL,+yDf ,0.0,charLength};
pFinest_SW  = newp; Point(pFinest_SW)  = {-xDfL,-yDf ,0.0,charLength};

pFinest_ShockNE  = newp; Point(pFinest_ShockNE)  = {+xDSh,+yDf ,0.0,charLength};
pFinest_ShockNW  = newp; Point(pFinest_ShockNW)  = {-xDSh,+yDf ,0.0,charLength};
pFinest_ShockSE  = newp; Point(pFinest_ShockSE)  = {+xDSh,-yDf ,0.0,charLength};
pFinest_ShockSW  = newp; Point(pFinest_ShockSW)  = {-xDSh,-yDf ,0.0,charLength};

l21 = newl; Line(l21) = {pFinest_SW     , pFinest_ShockSW};
l22 = newl; Line(l22) = {pFinest_ShockSW, pFinest_ShockSE};
l23 = newl; Line(l23) = {pFinest_ShockSE, pFinest_SE     };
l24 = newl; Line(l24) = {pFinest_SE     , pFinest_NE     };
l25 = newl; Line(l25) = {pFinest_NE     , pFinest_ShockNE};
l26 = newl; Line(l26) = {pFinest_ShockNE, pFinest_ShockNW};
l27 = newl; Line(l27) = {pFinest_ShockNW, pFinest_NW     };
l28 = newl; Line(l28) = {pFinest_NW     , pFinest_SW     };

l30 = newl; Line(l30) = {pCoarse_ShockSW, pMedium_ShockSW};
l31 = newl; Line(l31) = {pCoarse_ShockSE, pMedium_ShockSE};
l32 = newl; Line(l32) = {pCoarse_ShockNW, pMedium_ShockNW};
l33 = newl; Line(l33) = {pCoarse_ShockNE, pMedium_ShockNE};

l34 = newl; Line(l34) = {pMedium_ShockSW, pFinest_ShockSW};
l35 = newl; Line(l35) = {pMedium_ShockSE, pFinest_ShockSE};
l36 = newl; Line(l36) = {pMedium_ShockNW, pFinest_ShockNW};
l37 = newl; Line(l37) = {pMedium_ShockNE, pFinest_ShockNE};

l38 = newl; Line(l38) = {pFinest_ShockSW, pFinest_ShockNW};
l39 = newl; Line(l39) = {pFinest_ShockSE, pFinest_ShockNE};

//  TFI Shock domains (bottom to top, counter-clockwise starting on lower segment)
ll1 = newl;  Line Loop(ll1) = {+l02,+l31,-l12,-l30};
ll2 = newl;  Line Loop(ll2) = {+l12,+l35,-l22,-l34};
ll3 = newl;  Line Loop(ll3) = {+l22,+l39,+l26,-l38};
ll4 = newl;  Line Loop(ll4) = {-l26,-l37,+l16,+l36};
ll5 = newl;  Line Loop(ll5) = {-l16,-l33,+l06,+l32};

s1 = news; Plane Surface(s1) = {ll1};
s2 = news; Plane Surface(s2) = {ll2};
s3 = news; Plane Surface(s3) = {ll3};
s4 = news; Plane Surface(s4) = {ll4};
s5 = news; Plane Surface(s5) = {ll5};

Transfinite Surface{s1} = {pCoarse_ShockSW,pCoarse_ShockSE,pMedium_ShockSE,pMedium_ShockSW} ;
Transfinite Surface{s2} = {pMedium_ShockSW,pMedium_ShockSE,pFinest_ShockSE,pFinest_ShockSW} ;
Transfinite Surface{s3} = {pFinest_ShockSW,pFinest_ShockSE,pFinest_ShockNE,pFinest_ShockNW} ;
Transfinite Surface{s4} = {pFinest_ShockNW,pFinest_ShockNE,pMedium_ShockNE,pMedium_ShockNW} ;
Transfinite Surface{s5} = {pMedium_ShockNW,pMedium_ShockNE,pCoarse_ShockNE,pCoarse_ShockNW} ;

Transfinite Line{+l02,+l12,+l22,-l26,-l16,-l06}  = nXShock            Using Bump        1.65;
Transfinite Line{+l30,+l31,+l32,+l33}            = nYCoarse2Medium    Using Progression 0.95;
Transfinite Line{+l34,+l35,+l36,+l37}            = nYMedium2Finest    Using Progression 0.96;
Transfinite Line{+l38,+l39}                      = nYFinest           Using Progression 1.00;

//  TFI Shock domains (Right and Left of shock, counter-clockwise starting on lower segment)
ll6 = newl;  Line Loop(ll6) = {+l23,+l24,+l25,-l39};
ll7 = newl;  Line Loop(ll7) = {+l21,+l38,+l27,+l28};

s6 = news; Plane Surface(s6) = {ll6};
s7 = news; Plane Surface(s7) = {ll7};

Transfinite Surface{s6} = {pFinest_ShockSE,pFinest_SE,pFinest_NE,pFinest_ShockNE} ;
Transfinite Surface{s7} = {pFinest_SW,pFinest_ShockSW,pFinest_ShockNW,pFinest_NW} ;

Transfinite Line{+l24,-l28}                     = nYFinest           Using Progression 1.00;
Transfinite Line{-l23,+l25}                     = nXFinestR          Using Progression 0.985;
Transfinite Line{+l21,-l27}                     = nXFinestL          Using Progression 0.985;

//  Unstructured Intermediate and coarse Shock domains (Right and Left of shock, counter-clockwise starting on lower segment)

//  Medium
ll8 = newl;  Line Loop(ll8) = {+l13,+l14,+l15,+l37,-l25,-l24,-l23,-l35};
ll9 = newl;  Line Loop(ll9) = {+l11,+l34,-l21,-l28,-l27,-l36,+l17,+l18};

s8  = news; Plane Surface(s8)  = {ll8};
s9  = news; Plane Surface(s9)  = {ll9};

Transfinite Line{+l14}                          = nYMediumR          Using Bump        1.10;
Transfinite Line{+l18}                          = nYMediumL          Using Bump        2.50;
Transfinite Line{-l13,+l15}                     = nXMediumR          Using Progression 0.95;
Transfinite Line{+l11,-l17}                     = nXMediumL          Using Progression 0.95;

//  Coarse
ll10 = newl;  Line Loop(ll10) = {+l03,+l04,+l05,+l33,-l15,-l14,-l13,-l31};
ll11 = newl;  Line Loop(ll11) = {+l01,+l30,-l11,-l18,-l17,-l32,+l07,+l08};

s10 = news; Plane Surface(s10) = {ll10};
s11 = news; Plane Surface(s11) = {ll11};

Transfinite Line{+l04}                          = nYCoarseR          Using Bump        1.10;
Transfinite Line{+l08}                          = nYCoarseL          Using Bump        1.10;
Transfinite Line{-l03,+l05}                     = nXCoarseR          Using Progression 0.94;
Transfinite Line{+l01,-l07}                     = nXCoarseL          Using Progression 0.94;


Recombine Surface(s1);
Recombine Surface(s2);
Recombine Surface(s3);
Recombine Surface(s4);
Recombine Surface(s5);
Recombine Surface(s6);
Recombine Surface(s7);
Recombine Surface(s8);
Recombine Surface(s9);
Recombine Surface(s10);
Recombine Surface(s11);


// Make a list of all surfaces which need to be extruded
allParts[] = {s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11};

// Extrude geometry and quadrilateral mesh in Z direction //
zdir[] = Extrude{0, 0, 0.05} { Surface{allParts[]}; Layers{1}; Recombine;};

// Z = 0 ;  Periodic Surface
Physical Surface(379) = {44, 43, 42, 47, 48, 52, 51, 41, 56, 55, 40};
// Z = 1 ;  Periodic Surface
Physical Surface(380) = {166, 144, 210, 122, 188, 252, 294, 100, 336, 378, 78};
// Y = -L ;  Periodic Surface
Physical Surface(381) = {373, 161, 315};
// Y = +L ;  Periodic Surface
Physical Surface(382) = {349, 65, 307};
// X = +0 ;  Dirichlet inflow Surface
Physical Surface(383) = {377};
// X = +1 ;  Dirichlet outflow Surface
Physical Surface(384) = {311};
// Volumes
Physical Volume(385) = {5, 11, 9, 4, 7, 3, 6, 8, 2, 10, 1};
