/* -------------------------------------------------------------------
   Tutorial 8a : circuit coupling - resistive circuit

   Features:
   - Electrokinetic formulation coupled with resistive circuit
   - Definition of circuit elements (sources, resistances)
   - Implementation of a netlist

   To compute the solution in a terminal:
       getdp R_circuit -solve static -pos Cube_Top_Values_ASCII

   To compute the solution interactively from the Gmsh GUI:
       File > Open > R_circuit.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/*
  This example shows how to implement a circuit with discrete elements in a
  finite element simulation. Only lumped resistors, voltage- and current-sources
  are used. There are four separate cubes, each with a certain conductivity. The
  bottom of all cubes is kept at a constant voltage via the constraint
  Voltage_3D. The top of the cubes are connected with a circuit according the
  `R_circuit.jpg' schematic.
*/

Include "R_circuit_common.pro";

// FEM domain group data
Group {
  // Surfaces
  Top1    = Region[Cube1Top];
  Bottom1 = Region[Cube1Bottom];
  Top2    = Region[Cube2Top];
  Bottom2 = Region[Cube2Bottom];
  Top3    = Region[Cube3Top];
  Bottom3 = Region[Cube3Bottom];
  Top4    = Region[Cube4Top];
  Bottom4 = Region[Cube4Bottom];

  // Volumes
  Cube1 = Region[Cube1];
  Cube2 = Region[Cube2];
  Cube3 = Region[Cube3];
  Cube4 = Region[Cube4];

  // FEM Electrical regions
  // * all volumes:
  Vol_Ele = Region[{Cube1, Cube2, Cube3, Cube4}];
  // * all surfaces connected to the lumped element circuit
  Sur_Ele = Region[{Top1, Bottom1, Top2, Bottom2, Top3, Bottom3, Top4, Bottom4}];
  // * total electrical computation domain:
  Dom_Ele = Region[{Vol_Ele, Sur_Ele}];
}

// Circuit group data (fictive regions for the circuit elements)
Group {
  // Sources
  CurrentSource1 = Region[{3001}];
  VoltageSource1 = Region[{3002}];

  SourceI_Cir = Region[{CurrentSource1}]; // all current sources
  SourceV_Cir = Region[{VoltageSource1}]; // all voltage sources
  Sources_Cir = Region[{SourceI_Cir, SourceV_Cir}];

  // Resistors
  Res1 = Region[{4001}];
  Res2 = Region[{4002}];
  Res3 = Region[{4003}];
  Res4 = Region[{4004}];

  Resistance_Cir = Region[{Res1, Res2, Res3, Res4}]; // All resistors

  // Complete circuit domain containing all circuit elements
  Domain_Cir = Region[{Resistance_Cir, SourceI_Cir, SourceV_Cir}];
}

// FEM domain function data
Function {
  // Absolute voltage at the bottom of all cubes; set via constraint Voltage_3D,
  // not via a voltage source
  Vbottom = 4.0;

  // Geometry of a cube
  w = 1.0; // Width [m]
  l = 1.0; // Length [m]
  h = 1.0; // Height [m]

  // Resistance
  Rcube1 = 20.0; // Resistance of cube 1 [Ohm]
  Rcube2 = 40.0; // Resistance of cube 2 [Ohm]
  Rcube3 = 60.0; // Resistance of cube 3 [Ohm]
  Rcube4 = 10.0; // Resistance of cube 4 [Ohm]

  // Specific electrical conductivity [S*m/m^2]
  kappa[Cube1] = h / (w * l * Rcube1);
  kappa[Cube2] = h / (w * l * Rcube2);
  kappa[Cube3] = h / (w * l * Rcube3);
  kappa[Cube4] = h / (w * l * Rcube4);
}

// Circuit domain function data
Function {
  // External impedances functions
  R[Res1] = 10.0;
  R[Res2] = 20.0;
  R[Res3] = 15.0;
  R[Res4] = 10.0;
}

// FEM domain constraint data
Constraint {
  { Name Current_3D ; Type Assign ;
    Case {
    }
  }
  { Name Voltage_3D ; Type Assign ;
    Case {
      { Region Bottom1; Type Assign; Value Vbottom; }
      { Region Bottom2; Type Assign; Value Vbottom; }
      { Region Bottom3; Type Assign; Value Vbottom; }
      { Region Bottom4; Type Assign; Value Vbottom; }
    }
  }
}

// Circuit domain constraint data
Constraint {
  { Name Current_Cir ;
    Case {
      { Region CurrentSource1; Value 1.0; } // CurrentSource1 has 1.0 A
    }
  }
  { Name Voltage_Cir ;
    Case {
      { Region VoltageSource1; Value 10.0; } // VoltageSource1 has 10.0 V
    }
  }
  { Name ElectricalCircuit ; Type Network ;
    Case Circuit1 {
      // Here the circuit netlist is implemented. The important thing is that
      // each surface of the 3D geometry which is connected to the circuit
      // (Sur_Ele) has a single voltage and current, but in the netlist they
      // appear with two terminals (two nodes). The voltage between the two
      // nodes corresponds to the absolute voltage of the surface. The current
      // through the two nodes is the current flowing through the surface.  In
      // this simple case here we need the voltage drop across a cube (e.g. cube
      // 1) in the netlist rather than the absolute voltage. Therefore the
      // bottom-voltage has to be subtracted from the top-voltage. This is
      // accomplished by connecting both surfaces anti-serial (e.g. at node 4;
      // The voltage drop across cube 1 is then between node 2 and node 10).

      { Region VoltageSource1; Branch{1, 10}; }
      { Region Res1;           Branch{1,  2}; }
      // Voltage between node 2 and node 4 corresponds to the absolute voltage
      // of Top1:
      { Region Top1;           Branch{2,  4}; }
      // Note the reverse order of the nodes here - it's for subtracting the
      // bottom voltage:
      { Region Bottom1;        Branch{10, 4}; }
      // Voltage between node 3 and node 5 corresponds to the absolute voltage
      // of Top2:
      { Region Res2;           Branch{2,  3}; }
      { Region Top2;           Branch{3,  5}; }
      // Note the reverse order of the nodes here! It's for subtracting the
      // bottom voltage.:
      { Region Bottom2;        Branch{10, 5}; }

      { Region CurrentSource1; Branch{6, 10}; }
      { Region Res3;           Branch{6, 10}; }
      // Voltage between node 6 and node 8 corresponds to the absolute voltage
      // of Top3:
      { Region Top3;           Branch{6,  8}; }
      // Note the reverse order of the nodes here! It's for subtracting the
      // bottom voltage.:
      { Region Bottom3;        Branch{10, 8}; }
      { Region Res4;           Branch{6,  7}; }
      // Voltage between node 7 and node 9 corresponds to the absolute voltage
      // of Top4:
      { Region Top4;           Branch{7,  9}; }
      // Note the reverse order of the nodes here! It's for subtracting the
      // bottom voltage:
      { Region Bottom4;        Branch{10, 9}; }
    }
  }
}

Jacobian {
  { Name JVol ;
    Case {
      { Region Region[{Vol_Ele}]; Jacobian Vol; }
      { Region Region[{Sur_Ele}]; Jacobian Sur; }
    }
  }
}

Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints  1 ; }
          { GeoElement Line        ; NumberOfPoints  5 ; }
          { GeoElement Triangle    ; NumberOfPoints  6 ; }
          { GeoElement Quadrangle  ; NumberOfPoints  7 ; }
          { GeoElement Tetrahedron ; NumberOfPoints  15 ; }
          { GeoElement Hexahedron  ; NumberOfPoints  34 ; }
          { GeoElement Prism       ; NumberOfPoints  21 ; }
        }
      }
    }
  }
}


FunctionSpace {
  { Name Hgrad_V; Type Form0;
    BasisFunction {
      // All nodes but the grouped nodes of the surfaces for connecting the
      // circuitry
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Ele; Entity NodesOf[ All, Not Sur_Ele ]; }
      // Grouped nodes: Each surface of Sur_Ele has just one single voltage ->
      // each surface is just one node effectively
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes;
        Support Dom_Ele; Entity GroupsOfNodesOf[ Sur_Ele ]; }
    }
    GlobalQuantity {
      // Voltage of the surfaces of Sur_Ele
      { Name U; Type AliasOf; NameOfCoef vf; }
      // Current flowing through the surfaces
      { Name I; Type AssociatedWith; NameOfCoef vf; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf ; NameOfConstraint Voltage_3D; }
      { NameOfCoef U; EntityType GroupsOfNodesOf ; NameOfConstraint Voltage_3D; }
      { NameOfCoef I; EntityType GroupsOfNodesOf ; NameOfConstraint Current_3D; }
    }
  }

  { Name Hregion_Cir; Type Scalar;
    BasisFunction {
      { Name sn; NameOfCoef ir; Function BF_Region;
        Support Domain_Cir; Entity Domain_Cir; }
    }
    GlobalQuantity {
      { Name Iz; Type AliasOf; NameOfCoef ir; }
      { Name Uz; Type AssociatedWith; NameOfCoef ir; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }
}

Formulation {
  { Name static; Type FemEquation;
    Quantity {
      { Name V;  Type Local; NameOfSpace Hgrad_V; }
      { Name U;  Type Global; NameOfSpace Hgrad_V [U]; }
      { Name I;  Type Global; NameOfSpace Hgrad_V [I]; }
      { Name Uz; Type Global; NameOfSpace Hregion_Cir [Uz]; }
      { Name Iz; Type Global; NameOfSpace Hregion_Cir [Iz]; }
    }
    Equation {
      // FEM domain
      Integral   { [ kappa[] * Dof{d V}, {d V} ];
        In Vol_Ele; Integration I1; Jacobian JVol;  }
      GlobalTerm { [ Dof{I} , {U} ]; In Sur_Ele; }

      // Circuit related terms: Resistance equation
      GlobalTerm { [ Dof{Uz}, {Iz} ]; In Resistance_Cir; }
      GlobalTerm { [ R[] * Dof{Iz}, {Iz} ]; In Resistance_Cir; }

      // Here is the connection between the circuit and the finite element simulation
      GlobalEquation { Type Network; NameOfConstraint ElectricalCircuit;
        { Node {I};  Loop {U};  Equation {I};  In Sur_Ele; }
        { Node {Iz}; Loop {Uz}; Equation {Uz}; In Domain_Cir; }
      }
    }
  }
}

Resolution {
  { Name static;
    System {
      { Name A; NameOfFormulation static; }
    }
    Operation {
      Generate[A];
      Solve[A];
      SaveSolution[A];
    }
  }
}

PostProcessing {
  { Name Ele; NameOfFormulation static;
    Quantity {
      { Name V; Value{ Local{ [ {V} ]; In Vol_Ele; Jacobian JVol;} } }
      { Name Jv; Value{ Local{ [ -kappa[]*{d V} ]; In Vol_Ele; Jacobian JVol;} } }
      { Name J; Value{ Local{ [ Norm[kappa[]*{d V}] ];  In Vol_Ele; Jacobian JVol;} } }
      { Name U; Value { Term { [ {U} ]; In Sur_Ele; } } }
      { Name I; Value { Term { [ {I} ]; In Sur_Ele; } } }
    }
  }
}

PostOperation {
  // Absolute voltage [V]
  { Name V; NameOfPostProcessing Ele;
    Operation {
      Print[ V, OnElementsOf Vol_Ele, File "Result_V.pos"];
    }
  }
  // Current density vectors [A/m^2]
  { Name J_vectors; NameOfPostProcessing Ele;
    Operation {
      Print[ Jv, OnElementsOf Vol_Ele, File "Result_J_vectors.pos"];
    }
  }
  // Magnitude of current density [A/m^2]
  { Name J_magnitude; NameOfPostProcessing Ele;
    Operation {
      Print[ J, OnElementsOf Vol_Ele, File "Result_J.pos"];
    }
  }
  // Store the results in an ASCII file
  { Name Cube_Top_Values_ASCII; NameOfPostProcessing Ele;
    Operation {
      Print[U, OnRegion Top1, File   "Result_CubeTopValues.txt", Format SimpleTable];
      Print[U, OnRegion Top2, File > "Result_CubeTopValues.txt", Format SimpleTable];
      Print[U, OnRegion Top3, File > "Result_CubeTopValues.txt", Format SimpleTable];
      Print[U, OnRegion Top4, File > "Result_CubeTopValues.txt", Format SimpleTable];

      Print[I, OnRegion Top1, File > "Result_CubeTopValues.txt", Format SimpleTable];
      Print[I, OnRegion Top2, File > "Result_CubeTopValues.txt", Format SimpleTable];
      Print[I, OnRegion Top3, File > "Result_CubeTopValues.txt", Format SimpleTable];
      Print[I, OnRegion Top4, File > "Result_CubeTopValues.txt", Format SimpleTable];
    }
  }
}
