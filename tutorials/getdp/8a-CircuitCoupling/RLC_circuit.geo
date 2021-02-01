// Just 3 cubes

Include "RLC_circuit_common.pro";

lc = 0.25; // Characteristic length

Point(1) = {0, 0, 0, lc};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {0, 1, 0} {
  Line{1};
}
Extrude {0, 0, 1} {
  Surface{5};
}
Translate {2, 0, 0} {
  Duplicata { Volume{1}; }
}
Translate {4, 0, 0} {
  Duplicata { Volume{1}; }
}

Physical Surface(Cube1Top)      = {27};
Physical Surface(Cube1Bottom)   = {5};
Physical Surface(Cube2Top)      = {34};
Physical Surface(Cube2Bottom)   = {29};
Physical Surface(Cube3Top)      = {61};
Physical Surface(Cube3Bottom)   = {56};

Physical Volume(Cube1) = {1};
Physical Volume(Cube2) = {28};
Physical Volume(Cube3) = {55};
