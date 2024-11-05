// Gmsh project created on Thu May 25 11:54:54 2023
SetFactory("OpenCASCADE");
//+
Disk(1) = {0, 0, 0, 12, 12};
//+
Point(2) = {0, 0, 0, 0.01};
//+
Point(3) = {0, 0, 8, 0.01};
//+
Point(4) = {3, 0, 7.5, 0.01};
//+
Point(5) = {4.5, 0, 4.5, 0.01};
//+
Point(6) = {11, 0, 2, 0.01};
//+
Point(7) = {1, 0, 7.9, 0.01};
//+
Point(8) = {5.5, 0, 2.5, 0.01};
//+
BSpline(2) = {3, 7, 4, 5, 8, 6, 1};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} {
  Curve{2}; Layers{5}; Recombine;
}
//+
Coherence;
//+
Surface Loop(10) = {1,2};
//+
Volume(11) = {10};
//+
Physical Volume("cell", 12) = {11};
Mesh.MeshSizeMin = 0.5;
Mesh.MeshSizeMax = 1.5;
//+
Hide "*";
//+
Show {
  Point{2}; Point{7}; Point{8}; Curve{3}; Surface{1}; Surface{2}; Volume{11}; 
}
//+
Hide "*";
//+
Show {
  Point{2}; Point{7}; Point{8}; Curve{3}; Surface{1}; Surface{2}; Volume{11}; 
}
