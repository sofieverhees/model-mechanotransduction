// Gmsh project created on Thu May 25 11:54:54 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {24, 0, 0, 1.0};
//+
Point(3) = {0, 24, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Line(2) = {1, 2};
//+
Line(3) = {1, 3};
//+
Point(4) = {0, 0, 5.5, 0.01};
//+
Line(4) = {1, 4};
//+
Point(5) = {5, 0, 7.5, 0.01};
//+
Point(6) = {7.5, 0, 4.5, 0.01};
//+
Point(7) = {22, 0, 2, 0.01};
//+
Point(8) = {2, 0, 7.9, 0.01};
//+
Point(9) = {9.5, 0, 2.5, 0.01};
//+
BSpline(5) = {4, 8, 5, 6, 9, 7, 2};
//+
Point(10) = {0, 5, 7.5, 0.01};
//+
Point(11) = {0, 7.5, 4.5, 0.01};
//+
Point(12) = {0, 22, 2, 0.01};
//+
Point(13) = {0, 2, 7.9, 0.01};
//+
Point(14) = {0, 9.5, 2.5, 0.01};
//+
BSpline(6) = {4, 13, 10, 11, 14, 12, 3};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, 0.5*Pi} {Curve{5}; }
//+
Curve Loop(2) = {1,2,3};
//+
Curve Loop(3) = {5,2,4};
//+
Curve Loop(4) = {6,3,4};
//+ 
Surface(2) = {2};
//+ 
Surface(3) = {3};
//+
Surface(4) = {4};
//+
Coherence;
//+
Surface Loop(1) = {1,2,3,4};
//+ 
Volume(1) = {1};
//+
Physical Volume("cell",9) = {1};
//+
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 1.5;
