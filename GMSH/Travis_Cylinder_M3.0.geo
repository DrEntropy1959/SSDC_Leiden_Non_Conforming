// parameters
Mesh.FlexibleTransfinite  = 1;
Mesh.SubdivisionAlgorithm = 1;


sqrt2Inv   = 1.0/Sqrt(2)  ;
charLength = 1.0 ;

////////////////////////////////////   POINTS   /////////////////////////////////////////////////////

xCenter  =  0.00;
yCenter  =  0.00;
radius   =  0.50;
bndSize  =  0.15;
shkSize1 = +0.00;   //  X position of M shock connector (larger + moves to left)
shkSize3 =  1.25;   //  Y position of M shock connector
shkSize2 =  0.50;   //  X position of W shock connector
shkSize4 =  1.50;   //  Y position of W shock connector
shkSize5 =  0.35;
shkSize6 =  0.50;

radPbndS  = radius + bndSize;
radPbndS1 = radius + shkSize1;
radPbndS2 = radius + shkSize2;
radPbndS3 = radius + shkSize3;
radPbndS4 = radius + shkSize4;
radPbndS5 = radius + shkSize5;
radPbndS6 = radius + shkSize6;

xMin = -3;
xMax = +20;
yMin = -3;
yMax = +3;

nBL     =  12 ;
ncylW   =  17 ;
ncylE   =  12 ;
ncylNS  =  15 ;
nNearY  =  20 ;
nInX    =  20 ;
nInY    =  25 ;
nFarX   = 100 ;
nFarXb2 =  70 ;
nFarY   =  18 ;
nFarYb2 =  10 ;
nShock  =  08 ;

pSE = 1.1  ;
pSM = 1.0  ;
pSW = 0.05 ;
xM  = xCenter ;
yM  = yCenter ;
xSE = xCenter + pSE ;
xSM = xCenter + pSM ;
xSW = xCenter + pSW ;
ySE = yCenter ;
ySM = yCenter ;
ySW = yCenter ;
// define points
pCenter   = newp; Point(pCenter  ) = {xM  ,yM  ,0.0,charLength};  //p1
pCenterSE = newp; Point(pCenterSE) = {xSE ,ySE ,0.0,charLength};  //p2
pCenterSM = newp; Point(pCenterSM) = {xSM ,ySM ,0.0,charLength};  //p3
pCenterSW = newp; Point(pCenterSW) = {xSW ,ySW ,0.0,charLength};  //p4

pCirc1  = newp; Point(pCirc1)  = {xCenter+sqrt2Inv*radius,yCenter-sqrt2Inv*radius,0.0,charLength};      //p5
pCirc2  = newp; Point(pCirc2)  = {xCenter+sqrt2Inv*radius,yCenter+sqrt2Inv*radius,0.0,charLength};      //p6
pCirc3  = newp; Point(pCirc3)  = {xCenter-sqrt2Inv*radius,yCenter+sqrt2Inv*radius,0.0,charLength};      //p7
pCirc4  = newp; Point(pCirc4)  = {xCenter-sqrt2Inv*radius,yCenter-sqrt2Inv*radius,0.0,charLength};      //p8

pBndL1  = newp; Point(pBndL1)  = {xCenter+sqrt2Inv*radPbndS,yCenter-sqrt2Inv*radPbndS,0.0,charLength};  //p9
pBndL2  = newp; Point(pBndL2)  = {xCenter+sqrt2Inv*radPbndS,yCenter+sqrt2Inv*radPbndS,0.0,charLength};  //p10
pBndL3  = newp; Point(pBndL3)  = {xCenter-sqrt2Inv*radPbndS,yCenter+sqrt2Inv*radPbndS,0.0,charLength};  //p11
pBndL4  = newp; Point(pBndL4)  = {xCenter-sqrt2Inv*radPbndS,yCenter-sqrt2Inv*radPbndS,0.0,charLength};  //p12

x1  = xCenter - sqrt2Inv*radPbndS1 ;
y1m = yCenter - sqrt2Inv*radPbndS3 ;
y1p = yCenter + sqrt2Inv*radPbndS3 ;
x2  = xCenter - sqrt2Inv*radPbndS2 ;
y2m = yCenter - sqrt2Inv*radPbndS4 ;
y2p = yCenter + sqrt2Inv*radPbndS4 ;

//  other quadratic root,  may be useful?
//x5 = (xM + xSE + yM - ySE + Sqrt(2*x1*x1 - xM*xM - 4*x1*xSE + 2*xM*xSE + xSE*xSE + 2*y1p*y1p - 2*xM*yM + 2*xSE*yM - yM*yM + 2*(xM - xSE - 2*y1p + yM)*ySE + ySE*ySE))/2 ;
//y5 = (xM - xSE + yM + ySE - Sqrt(2*x1*x1 - xM*xM - 4*x1*xSE + 2*xM*xSE + xSE*xSE + 2*y1p*y1p - 2*xM*yM + 2*xSE*yM - yM*yM + 2*(xM - xSE - 2*y1p + yM)*ySE + ySE*ySE))/2 ;

x5p = (xM + xSE + yM - ySE - Sqrt(2*x1*x1 - xM*xM - 4*x1*xSE + 2*xM*xSE + xSE*xSE + 2*y1p*y1p - 2*xM*yM + 2*xSE*yM - yM*yM + 2*(xM - xSE - 2*y1p + yM)*ySE + ySE*ySE))/2 ;
y5p = (xM - xSE + yM + ySE + Sqrt(2*x1*x1 - xM*xM - 4*x1*xSE + 2*xM*xSE + xSE*xSE + 2*y1p*y1p - 2*xM*yM + 2*xSE*yM - yM*yM + 2*(xM - xSE - 2*y1p + yM)*ySE + ySE*ySE))/2 ;
x5m = +x5p ;
y5m = -y5p ;

x6p = (xM + xSM + yM - ySM - Sqrt(2*x2*x2 - xM*xM - 4*x2*xSM + 2*xM*xSM + xSM*xSM + 2*y2p*y2p - 2*xM*yM + 2*xSM*yM - yM*yM + 2*(xM - xSM - 2*y2p + yM)*ySM + ySM*ySM))/2 ;
y6p = (xM - xSM + yM + ySM + Sqrt(2*x2*x2 - xM*xM - 4*x2*xSM + 2*xM*xSM + xSM*xSM + 2*y2p*y2p - 2*xM*yM + 2*xSM*yM - yM*yM + 2*(xM - xSM - 2*y2p + yM)*ySM + ySM*ySM))/2 ;
x6m = +x6p ;
y6m = -y6p ;

//  Middle shock connector
pBndL5  = newp; Point(pBndL5)   = {x1,y1p,0.0,charLength} ;                                          //p13    
pBndL6  = newp; Point(pBndL6)   = {x5p,y5p,0.0,charLength};                                          //p14
pBndL7  = newp; Point(pBndL7)   = {x5m,y5m,0.0,charLength};                                          //p15
pBndL8  = newp; Point(pBndL8)   = {x1,y1m,0.0,charLength} ;                                          //p16

//  West shock connector
pBndL9  = newp; Point(pBndL9)   = {x2,y2p,0.0,charLength} ;                                          //p17
pBndL10 = newp; Point(pBndL10)  = {x6p,y6p,0.0,charLength};                                          //p18
pBndL11 = newp; Point(pBndL11)  = {x6m,y6m,0.0,charLength};                                          //p19
pBndL12 = newp; Point(pBndL12)  = {x2,y2m,0.0,charLength} ;                                          //p20

// Four corners of the larger domain
pFarF1  = newp; Point(pFarF1)  = {xMax              ,yMin                   ,0.0,charLength};        //p21
pFarF2  = newp; Point(pFarF2)  = {xMax              ,yMax                   ,0.0,charLength};        //p22
pFarF3  = newp; Point(pFarF3)  = {xMin              ,yMax                   ,0.0,charLength};        //p23
pFarF4  = newp; Point(pFarF4)  = {xMin              ,yMin                   ,0.0,charLength};        //p24


// Upper and lower shock bounce coordinates 
BxW  = 1.00;
BxE  = 1.70;
BWid = 0.2 ;
BHgh = 0.35 ;
pB1  = newp; Point(pB1)  = {xCenter+BxW,+yMax,0.0,charLength};                                       //p25
pB2  = newp; Point(pB2)  = {xCenter+BxW,-yMax,0.0,charLength};                                       //p26
pB3  = newp; Point(pB3)  = {xCenter+BxE,+yMax,0.0,charLength};                                       //p27
pB4  = newp; Point(pB4)  = {xCenter+BxE,-yMax,0.0,charLength};                                       //p28

pB5  = newp; Point(pB5)  = {xCenter+BxW+BWid,+yMax-BHgh,0.0,charLength};                             //p29
pB6  = newp; Point(pB6)  = {xCenter+BxW+BWid,-yMax+BHgh,0.0,charLength};                             //p30
pB7  = newp; Point(pB7)  = {xCenter+BxE-BWid,+yMax-BHgh,0.0,charLength};                             //p31
pB8  = newp; Point(pB8)  = {xCenter+BxE-BWid,-yMax+BHgh,0.0,charLength};                             //p32

// Vertical line dividing farfield from nearfield
pV1  = newp; Point(pV1)  = {5.0,yCenter+sqrt2Inv*radPbndS+0.00,0.0,charLength};                      //p33
pV2  = newp; Point(pV2)  = {5.0,yCenter-sqrt2Inv*radPbndS-0.00,0.0,charLength};                      //p34
pV3  = newp; Point(pV3)  = {5.0,+yMax,0.0,charLength};                                              //p35
pV4  = newp; Point(pV4)  = {5.0,-yMax,0.0,charLength};                                              //p36
pV5  = newp; Point(pV5)  = {5.0,yCenter+sqrt2Inv*radPbndS+0.55,0.0,charLength};                     //p37
pV6  = newp; Point(pV6)  = {5.0,yCenter-sqrt2Inv*radPbndS-0.55,0.0,charLength};                     //p38

// Vertical line at Xmax
pV7  = newp; Point(pV7)  = {xMax,yCenter+sqrt2Inv*radPbndS+0.15,0.0,charLength};                    //p39
pV8  = newp; Point(pV8)  = {xMax,yCenter-sqrt2Inv*radPbndS-0.15,0.0,charLength};                    //p40
pV9  = newp; Point(pV9)  = {xMax,yCenter+sqrt2Inv*radPbndS+0.80,0.0,charLength};                    //p41
pV10 = newp; Point(pV10) = {xMax,yCenter-sqrt2Inv*radPbndS-0.80,0.0,charLength};                    //p42

////////////////////////////////////   CONNECTORS   /////////////////////////////////////////////////////

// define Circle at radius = 0.5
lCirc1 = newl; Circle(lCirc1) = {pCirc1,pCenter,pCirc2};   //l1
lCirc2 = newl; Circle(lCirc2) = {pCirc2,pCenter,pCirc3};   //l2
lCirc3 = newl; Circle(lCirc3) = {pCirc3,pCenter,pCirc4};   //l3
lCirc4 = newl; Circle(lCirc4) = {pCirc4,pCenter,pCirc1};   //l4

//  Normal connectors
lBndC1 = newl; Line(lBndC1) = {pCirc1,pBndL1};             //l5
lBndC2 = newl; Line(lBndC2) = {pCirc2,pBndL2};             //l6
lBndC3 = newl; Line(lBndC3) = {pCirc3,pBndL3};             //l7
lBndC4 = newl; Line(lBndC4) = {pCirc4,pBndL4};             //l8

//  Concentric ``circle''
lBndL1 = newl; Circle(lBndL1) = {pBndL1,pCenter,pBndL2};   //l9
lBndL2 = newl; Circle(lBndL2) = {pBndL2,pCenter,pBndL3};   //l10
lBndL3 = newl; Circle(lBndL3) = {pBndL3,pCenterSW,pBndL4}; //l11
lBndL4 = newl; Circle(lBndL4) = {pBndL4,pCenter,pBndL1};   //l12

//  East Innermost shock circumferance
lShkM1 = newl; Circle(lShkM1) = {pBndL5,pCenterSE,pBndL6}; //l13
lShkM2 = newl; Circle(lShkM2) = {pBndL6,pCenterSE,pBndL7}; //l14
lShkM3 = newl; Circle(lShkM3) = {pBndL7,pCenterSE,pBndL8}; //l15
//  East Outermost shock circumferance
lShkW1 = newl; Circle(lShkW1) = {pBndL9,pCenterSM,pBndL10};  //l16
lShkW2 = newl; Circle(lShkW2) = {pBndL10,pCenterSM,pBndL11}; //l17
lShkW3 = newl; Circle(lShkW3) = {pBndL11,pCenterSM,pBndL12}; //l18

// Orthogonal connectors between W and M shock circumferences
lBndS1 = newl; Line(lBndS1) = {pBndL5,pBndL9};               //l19
lBndS2 = newl; Line(lBndS2) = {pBndL6,pBndL10};              //l20
lBndS3 = newl; Line(lBndS3) = {pBndL7,pBndL11};              //l21
lBndS4 = newl; Line(lBndS4) = {pBndL8,pBndL12};              //l22

// Orthogonal connectors between M shock and inner circle circumferences
lBndS5 = newl; Line(lBndS5) = {pBndL3,pBndL6};               //l23
lBndS6 = newl; Line(lBndS6) = {pBndL4,pBndL7};               //l24

// Bounce Region connectors
lBndB1 = newl; Line(lBndB1) = {pB2,pB4};               //l25
lBndB2 = newl; Line(lBndB2) = {pB6,pB8};               //l26
lBndB3 = newl; Line(lBndB3) = {pB2,pB6};               //l27
lBndB4 = newl; Line(lBndB4) = {pB4,pB8};               //l28

lBndB5 = newl; Line(lBndB5) = {pB1,pB3};               //l29
lBndB6 = newl; Line(lBndB6) = {pB5,pB7};               //l30
lBndB7 = newl; Line(lBndB7) = {pB1,pB5};               //l31
lBndB8 = newl; Line(lBndB8) = {pB3,pB7};               //l32

// Circumferencial shock to bounce connectors
// Top  W -> E

lShkBnc1 = newl; Line(lShkBnc1) = {pBndL9,pB1};          //l33
lShkBnc2 = newl; Line(lShkBnc2) = {pBndL5,pB5};          //l34
// Bottom  W -> E
lShkBnc3 = newl; Line(lShkBnc3) = {pBndL12,pB2};         //l35
lShkBnc4 = newl; Line(lShkBnc4) = {pBndL8, pB6};         //l36

// Bottom to top connectors (W -> E)  Bounce points to vertical middle divider
lBncMid1 = newl; Line(lBncMid1) = {pB4,pV6};             //l37
lBncMid2 = newl; Line(lBncMid2) = {pB8,pV2};             //l38
lBncMid3 = newl; Line(lBncMid3) = {pB7,pV1};             //l39
lBncMid4 = newl; Line(lBncMid4) = {pB3,pV5};             //l40

LmidVert1 = newl; Line(LmidVert1) = {pV4,pV6} ;          //l41
LmidVert2 = newl; Line(LmidVert2) = {pV6,pV2} ;          //l42
LmidVert3 = newl; Line(LmidVert3) = {pV2,pV1} ;          //l43
LmidVert4 = newl; Line(LmidVert4) = {pV1,pV5} ;          //l44
LmidVert5 = newl; Line(LmidVert5) = {pV5,pV3} ;          //l45

LBottom1  = newl; Line(LBottom1)  = {pFarF4,pB2} ;       //l46
LBottom2  = newl; Line(LBottom2)  = {pB4,pV4} ;          //l47
LBottom3  = newl; Line(LBottom3)  = {pV4,pFarF1} ;       //l48

LTop1  = newl; Line(LTop1)        = {pFarF3,pB1} ;       //l49
LTop2  = newl; Line(LTop2)        = {pB3,pV3} ;          //l50
LTop3  = newl; Line(LTop3)        = {pV3,pFarF2} ;       //l51

LFarVert1 = newl; Line(LFarVert1) = {pFarF1,pV10} ;      //l52
LFarVert2 = newl; Line(LFarVert2) = {pV10,pV8} ;         //l53
LFarVert3 = newl; Line(LFarVert3) = {pV8,pV7} ;          //l54
LFarVert4 = newl; Line(LFarVert4) = {pV7,pV9} ;          //l55
LFarVert5 = newl; Line(LFarVert5) = {pV9,pFarF2} ;       //l56

LNearVert1 = newl; Line(LNearVert1) = {pFarF4,pFarF3} ;  //l57

LHorizFar1 = newl; Line(LHorizFar1) = {pV6,pV10} ;       //l58
LHorizFar2 = newl; Line(LHorizFar2) = {pV2,pV8} ;        //l59
LHorizFar3 = newl; Line(LHorizFar3) = {pV1,pV7} ;        //l60
LHorizFar4 = newl; Line(LHorizFar4) = {pV5,pV9} ;        //l61

///////////////////////////////////  Surfaces /////////////////////////////////////////////

//  Boundary layer lines
Transfinite Line{-lBndC1,-lBndC2,-lBndC3,-lBndC4}  = nBL    Using Progression 0.90;

//  East nearfield Cylinder
Transfinite Line{-lCirc1,+lBndL1}                  = ncylE  Using Bump 0.80;
LL1 = newl; Line Loop(LL1) = {lBndL1, -lBndC2, -lCirc1, lBndC1}; 
S1 = news; Plane Surface(S1) = {LL1}; Transfinite Surface(S1) = {pCirc1,pBndL1,pBndL2,pCirc2}; Recombine Surface(S1);

//  North
Transfinite Line{-lCirc2,+lBndL2}                 = ncylNS  Using Bump 0.80;
LL2 = newl; Line Loop(LL2) = {lBndL2, -lBndC3, -lCirc2, lBndC2}; 
S2 = news; Plane Surface(S2) = {LL2}; Transfinite Surface(S2) = {pCirc2,pBndL2,pBndL3,pCirc3}; Recombine Surface(S2);

//  West
Transfinite Line{-lCirc3,+lBndL3}                 = ncylW   Using Bump 0.80;
LL3 = newl; Line Loop(LL3) = {lBndL3, -lBndC4, -lCirc3, lBndC3};
S3 = news; Plane Surface(S3) = {LL3}; Transfinite Surface(S3) = {pCirc3,pBndL3,pBndL4,pCirc4}; Recombine Surface(S3);

//  South 
Transfinite Line{-lCirc4,+lBndL4}                 = ncylNS  Using Bump 0.80;
LL4 = newl; Line Loop(LL4) = {lBndL4, -lBndC1, -lCirc4, lBndC4};
S4 = news; Plane Surface(S4) = {LL4}; Transfinite Surface(S4) = {pCirc4,pBndL4,pBndL1,pCirc1}; Recombine Surface(S4);

//  West-West between shock and cylinder
Transfinite Line{lShkM2}                          = ncylW  Using Bump 0.8;
Transfinite Line{-lBndS5,-lBndS6}                   = 9  Using Progression 0.9;
LL5 = newl; Line Loop(LL5) = {-lBndL3, +lBndS5, +lShkM2,-lBndS6};
S5 = news; Plane Surface(S5) = {LL5}; Transfinite Surface(S5) = {pBndL4,pBndL6,pBndL7,pBndL3}; Recombine Surface(S5);

Transfinite Line{lShkW2}                          = ncylW  Using Bump 0.8;
LL6 = newl; Line Loop(LL6) = {-lShkM2, +lBndS2, +lShkW2,-lBndS3};
S6 = news; Plane Surface(S6) = {LL6}; Transfinite Surface(S6) = {pBndL7,pBndL6,pBndL10,pBndL11}; Recombine Surface(S6);

//  Small shoulder regions in shock between bow and bounce 
Transfinite Line{lShkM1,lShkW1,-lShkM3,-lShkW3}   = 7  Using Progression 0.99;
LL7 = newl; Line Loop(LL7) = {-lShkM1, lBndS1, lShkW1, -lBndS2};
S7 = news; Plane Surface(S7) = {LL7}; Transfinite Surface(S7) = {pBndL6,pBndL5,pBndL9,pBndL10}; Recombine Surface(S7);
LL8 = newl; Line Loop(LL8) = {-lShkM3, lBndS3, lShkW3, -lBndS4};
S8 = news; Plane Surface(S8) = {LL8}; Transfinite Surface(S8) = {pBndL8,pBndL7,pBndL11,pBndL12}; Recombine Surface(S8);

//  Shock to Bounce surfaces
nShkBnc = 20;
Transfinite Line{lShkBnc1,lShkBnc2,lShkBnc3,lShkBnc4}   = nShkBnc     Using Bump 0.85;
LL9  = newl; Line Loop(LL9)  = {lShkBnc2 ,-lBndB7 , -lShkBnc1 ,-lBndS1};
S9 = news; Plane Surface(S9) = {LL9}; Transfinite Surface(S9) = {pBndL9,pBndL5,pB5,pB1}; Recombine Surface(S9);

LL10 = newl; Line Loop(LL10) = {-lShkBnc4 , lBndS4 , lShkBnc3 , lBndB3 };
S10 = news; Plane Surface(S10) = {LL10}; Transfinite Surface(S10) = {pBndL8, pBndL12, pB2 ,pB6} ; Recombine Surface(S10);

//  Bounce surfaces
nBounce = 10;
Transfinite Line{lBndB1,lBndB2,lBndB5,lBndB6}       = nBounce Using Bump 1.0 ;
LL11 = newl; Line Loop(LL11) = {-lBndB6, lBndB8, lBndB5, -lBndB7} ;
S11 = news; Plane Surface(S11) = {LL11}; Transfinite Surface(S11) = {pB5, pB7, pB3, pB1}; Recombine Surface(S11);
LL12 = newl; Line Loop(LL12) = {lBndB1, lBndB4, -lBndB2, -lBndB3} ;
S12 = news; Plane Surface(S12) = {LL12}; Transfinite Surface(S12) = {pB2, pB4, pB8, pB6} ; Recombine Surface(S12);

//  Bounce to Mid TFI surfaces
nBncMid = 40;
Transfinite Line{lBncMid1,lBncMid2,lBncMid3,lBncMid4}   = nBncMid Using Progression 1.00 ;
LL13 = newl; Line Loop(LL13) = {lBncMid1 , LmidVert2 , -lBncMid2 , -lBndB4 } ;
S13 = news; Plane Surface(S13) = {LL13}; Transfinite Surface(S13) = {pB4, pV6, pV2, pB8} ; Recombine Surface(S13);
LL14 = newl; Line Loop(LL14) = {lBncMid3 , LmidVert4 , -lBncMid4 , lBndB8} ;
S14 = news; Plane Surface(S14) = {LL14}; Transfinite Surface(S14) = {pB7, pV1, pV5, pB3} ; Recombine Surface(S14);

//  Farfield TFI surfaces
nVertTBFar = 10 ;
nVertMFar  = 17 ;
Transfinite Line{LBottom3,LHorizFar1,LHorizFar2,LHorizFar3,LHorizFar4,LTop3}   = nFarX       Using Progression 1.010  ;
Transfinite Line{LFarVert1,LmidVert1,-LFarVert5,-LmidVert5}                    = nVertTBFar  Using Progression 0.9 ;
Transfinite Line{LFarVert3,LmidVert3}                                          = nVertMFar   Using Bump 1.0 ;
LL15 = newl; Line Loop(LL15) = {LBottom3  ,LFarVert1,-LHorizFar1,-LmidVert1} ;
LL16 = newl; Line Loop(LL16) = {LHorizFar1,LFarVert2,-LHorizFar2,-LmidVert2} ;
LL17 = newl; Line Loop(LL17) = {LHorizFar2,LFarVert3,-LHorizFar3,-LmidVert3} ;
LL18 = newl; Line Loop(LL18) = {LHorizFar3,LFarVert4,-LHorizFar4,-LmidVert4} ;
LL19 = newl; Line Loop(LL19) = {LHorizFar4,LFarVert5,-LTop3,-LmidVert5} ;

S15 = news; Plane Surface(S15) = {LL15}; Transfinite Surface(S15) = {pV4,pFarF1,pV10,pV6} ; Recombine Surface(S15);
S16 = news; Plane Surface(S16) = {LL16}; Transfinite Surface(S16) = {pV6,pV10,pV8,pV2} ; Recombine Surface(S16);
S17 = news; Plane Surface(S17) = {LL17}; Transfinite Surface(S17) = {pV2,pV8,pV7,pV1} ; Recombine Surface(S17);
S18 = news; Plane Surface(S18) = {LL18}; Transfinite Surface(S18) = {pV1,pV7,pV9,pV5} ; Recombine Surface(S18);
S19 = news; Plane Surface(S19) = {LL19}; Transfinite Surface(S19) = {pV5,pV9,pFarF2,pV3} ; Recombine Surface(S19);

Transfinite Line{LFarVert4, LmidVert4, lBndB8, lBndB7, lBndS1, lBndS2, lBndS3, lBndS4, lBndB3, lBndB4, LmidVert2, LFarVert2} = nShock Using Bump 1.0 ;
nVerticalIn = 15 ;
nHorizonIn  = 20 ;
nHorizonMid = 20 ;
Transfinite Line{LNearVert1}                                               = nVerticalIn Using Bump 1.2 ;
Transfinite Line{-LBottom1,-LTop1}                                         = nHorizonIn  Using Progression 1.1 ;
Transfinite Line{+LBottom2,+LTop2}                                         = nHorizonMid Using Progression 1.05 ;

LL20 = newl; Line Loop(LL20) = {+LBottom1,-lShkBnc3,-lShkW3,-lShkW2,-lShkW1,+lShkBnc1,-LTop1,-LNearVert1} ;
S20 = news; Plane Surface(S20) = {LL20} ; Recombine Surface{S20};

LL21 = newl; Line Loop(LL21) = {+LBottom2,+LmidVert1,-lBncMid1};
S21 = news; Plane Surface(S21) = {LL21} ; Recombine Surface{S21};

LL22 = newl; Line Loop(LL22) = {+lBncMid4,+LmidVert5,-LTop2} ;
S22 = news; Plane Surface(S22) = {LL22} ; Recombine Surface{S22};

LL23 = newl; Line Loop(LL23) = {+lBndB2,+lBncMid2,+LmidVert3,-lBncMid3,-lBndB6,-lShkBnc2,+lShkM1,-lBndS5,-lBndL2,-lBndL1,-lBndL4,+lBndS6,+lShkM3,+lShkBnc4} ;
S23 = news; Plane Surface(S23) = {LL23} ; Recombine Surface{S23};

// Make a list of all surfaces which need to be extruded
allParts[] = {S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23} ;


// Extrude geometry and quadrilateral mesh in Z direction //
zdir[] = Extrude{0, 0, 1} { Surface{allParts[]}; Layers{1}; Recombine;};

// periodic_1 Back
Physical Surface(674) = {83, 105, 79, 89, 99, 75, 98, 65, 63, 67, 97, 71, 69, 73, 77, 96, 87, 81, 107, 101, 95, 103, 85};
// periodic_1 Front
Physical Surface(675) = {349, 601, 305, 415, 525, 261, 503, 151, 129, 173, 481, 217, 195, 239, 283, 459, 393, 327, 673, 567, 437, 584, 371};
// uniform_free_stream inflow
Physical Surface(676) = {566};
// uniform_free_stream Outflow
Physical Surface(678) = {516, 494, 472, 450, 428};
// symmetry_plane Top
Physical Surface(679) = {424, 575, 358, 538};
// symmetry_plane Bottom
Physical Surface(680) = {520, 600, 344, 562};
// Solid Surface
Physical Surface(681) = {124, 146, 190, 168};
// Volumes
Physical Volume(682) = {20, 6, 8, 5, 7, 3, 4, 2, 10, 9, 1, 12, 11, 23, 13, 21, 14, 22, 15, 16, 17, 18, 19};
