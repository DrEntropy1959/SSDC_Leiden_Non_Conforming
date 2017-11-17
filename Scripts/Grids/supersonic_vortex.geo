// parameters
charLength = 1.0;
xCenter    = -1.0;
yCenter    =  0.0;
radius1    = 1.0;
radius2    = 1.384;

sqrt2 = Sqrt(2);
//sqrt2Inv = 1.0/1.414213562;
sqrt2Inv = 1.0/sqrt2 ;

// define points
pCenter = newp; Point(pCenter) = {xCenter                ,yCenter                ,0.0,charLength};

pCirc1  = newp; Point(pCirc1)  = {xCenter+sqrt2Inv*radius1,yCenter-sqrt2Inv*radius1,0.0,charLength};
pCirc2  = newp; Point(pCirc2)  = {xCenter+sqrt2Inv*radius1,yCenter+sqrt2Inv*radius1,0.0,charLength};

pCirc3  = newp; Point(pCirc3)  = {xCenter+sqrt2Inv*radius2,yCenter-sqrt2Inv*radius2,0.0,charLength};
pCirc4  = newp; Point(pCirc4)  = {xCenter+sqrt2Inv*radius2,yCenter+sqrt2Inv*radius2,0.0,charLength};

// define lines
lCirc1 = newl; Circle(lCirc1) = {pCirc1,pCenter,pCirc2};
lCirc2 = newl; Circle(lCirc2) = {pCirc3,pCenter,pCirc4};

lLinC1 = newl; Line(lLinC1) = {pCirc2,pCirc4};
lLinC2 = newl; Line(lLinC2) = {pCirc1,pCirc3};

ll1 = newl; Line Loop(ll1) = {-lCirc1,lLinC2,lCirc2,-lLinC1};

// define surfaces
s1 = news; Plane Surface(s1) = {ll1};

Transfinite Surface{s1} = {pCirc2,pCirc1,pCirc4,pCirc3};


q = 2 ;
p = 0 ;
Cylinder = 1 + q*2*(2^p) ;
EndCaps  = 1 + q*1*(2^p) ;
Transfinite Line{lLinC1,lLinC2}  = EndCaps   Using Progression 1.0;
Transfinite Line{lCirc1,lCirc2}  = Cylinder  Using Progression 1.0;
//Transfinite Line{lCirc1}  = Cylinder+0  Using Progression 1.0;
//Transfinite Line{lCirc2}  = Cylinder+1  Using Progression 1.0;

// define recombine surfaces
Recombine Surface(s1);

//// Make a list of all surfaces which need to be extruded
allParts[] = {s1};

//// Extrude geometry and quadrilateral mesh in Z direction //
zdir[] = Extrude{0, 0, 1} { Surface{allParts[]}; Layers{1}; Recombine;};

// Radius inner and outer
Physical Surface(29) = {15, 23};
// inflow and outflow
Physical Surface(30) = {27, 19};
// periodic front
Physical Surface(31) = {28};
// periodic back
Physical Surface(32) = {6};
// Single Volume 
Physical Volume(33) = {1};
