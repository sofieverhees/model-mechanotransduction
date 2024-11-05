// Gmsh project created on Thu May 25 11:54:54 2023
SetFactory("OpenCASCADE");
//+
Disk(1) = {0, 0, 0, 12, 12};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0, 0, 8, 1.0};
//+
Point(4) = {3, 0, 7.5, 1.0};
//+
Point(5) = {4.5, 0, 4.5, 1.0};
//+
Point(6) = {11, 0, 2, 1.0};
//+
Point(7) = {1, 0, 7.9, 1.0};
//+
Point(8) = {5.5, 0, 2.5, 1.0};
//+
BSpline(2) = {3, 7, 4, 5, 8, 6, 1};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} {
  Curve{2}; Layers{5}; Recombine;
}
//+
Point(10) = {0, 0, 5, 1.0};
//+
Point(11) = {0, 0, 7, 1.0};
//+
Point(12) = {0, 0, 3, 1.0};
//+
Point(13) = {2.5, 0, 5, 1.0};
//+
Ellipse(5) = {11,10,13};
//+
Ellipse(6) = {13,10,12};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} {
  Curve{5}; Layers{5}; Recombine;
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} {
  Curve{6}; Layers{5}; Recombine;
}
//+
Coherence;
//+
Surface Loop(10) = {1,2};
//+
Surface Loop(11) = {3,4};
//+
Volume(12) = {10,11};
//+
//Volume(13) = {11};
//+
//Physical Volume("cell", 14) = {12,13};
Physical Volume("cell", 14) = {12};
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 1.5;
