/* -------------------------------------------------------------------
   Tutorial 8b : circuit coupling - RLC circuit

   Features:
   - Electrokinetic formulation coupled with RLC circuit
   - Transient resolution

   To compute the solution in a terminal:
       getdp RLC_circuit -solve dynamic -pos Cube_Top_Values_ASCII

   To compute the solution interactively from the Gmsh GUI:
       File > Open > RLC_circuit.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/*
  This example shows how to implement a circuit with discrete resistors,
  capacitors, inductors, voltage- and current-sources in a transient finite
  element simulation. There are three separate cubes, each with a certain
  conductivity. The bottom of all cubes is kept at 0V via the constraint
  Voltage_3D. The cubes are connected according the `RLC_circuit.jpg` schematic.
*/

Include "RLC_circuit_common.pro";

// FEM domain group data
Group {
  // Surfaces
  Top1          = Region[Cube1Top];
  Bottom1       = Region[Cube1Bottom];
  Top2          = Region[Cube2Top];
  Bottom2       = Region[Cube2Bottom];
  Top3          = Region[Cube3Top];
  Bottom3       = Region[Cube3Bottom];

  // Volumes
  Cube1         = Region[Cube1];
  Cube2         = Region[Cube2];
  Cube3         = Region[Cube3];

  // FEM Electrical regions
  // * all volumes
  Vol_Ele   = Region[{Cube1, Cube2, Cube3}];
  // * all surfaces connected to the lumped element circuit
  Sur_Ele   = Region[{Top1, Bottom1, Top2, Bottom2, Top3, Bottom3}];
  // * total electrical computation domain
  Dom_Ele   = Region[{Vol_Ele, Sur_Ele}];
}

// Circuit group data - Fictive regions for the circuit elements
Group {

  // Sources
  CurrentSource1  = Region[{3001}];
  VoltageSource1  = Region[{3002}];

  // Resistors
  R1 = Region[{4001}];
  R2 = Region[{4002}];

  // Inductors
  L1 = Region[{4003}];

  // Capacitors
  C1 = Region[{4004}];

  Resistance_Cir  = Region[{R1, R2}]; // All resistors
  Inductance_Cir  = Region[{L1}]; // All inductors
  Capacitance_Cir = Region[{C1}]; // All capacitors
  SourceI_Cir     = Region[{CurrentSource1}]; // All current sources
  SourceV_Cir     = Region[{VoltageSource1}]; // All voltage sources

  // Complete circuit domain containing all circuit elements
  Domain_Cir = Region[{Resistance_Cir, Inductance_Cir, Capacitance_Cir,
                       SourceI_Cir, SourceV_Cir}];
}

Function {
  // Simulation parameters
  tStop  = 10.0e-6; // Stop time [s]
  tStep  = 100e-9;  // Time step [s]
  nStrobe= 1;       // Strobe periode for saving the result

  SaveFct[] = !($TimeStep % nStrobe); // Only each nStrobe-th time-step is saved

  // FEM domain function data
  // ------------------------

  Vbottom = 0.0; // Absolute voltage at the bottom of all cubes

  // Geometry of a cube
  w = 1.0; // Width [m]
  l = 1.0; // Length [m]
  h = 1.0; // Height [m]

  // Resistance
  Rcube1 = 10.0; // Resistance of cube 1 [Ohm]
  Rcube2 = 40.0; // Resistance of cube 2 [Ohm]
  Rcube3 = 20.0; // Resistance of cube 3 [Ohm]

  // Specific electrical conductivity [S*m/m^2]
  kappa[Cube1] = h / (w * l * Rcube1);
  kappa[Cube2] = h / (w * l * Rcube2);
  kappa[Cube3] = h / (w * l * Rcube3);

  // Circuit domain function data
  // ----------------------------

  // Resistors
  R[R1] = 30.0; // Resistance of R1 [Ohm]
  R[R2] = 40.0; // Resistance of R2 [Ohm]

  // Inductors
  L[L1] = 100.0e-6; // Inductance of L1 [H]

  // Capacitors
  C[C1] = 250.0e-9; // Capacitance of C1 [F]

  // Note: All voltages and currents are assigned / initialized in the
  // Constraint-block
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
    }
  }
}

// Circuit domain constraint data
Constraint {
  { Name Current_Cir ;
    Case {
      { Region CurrentSource1; Type Assign; Value 1.0; } // CurrentSource1 has 1.0 A
      { Region L1; Type Init; Value 1.0; } // Initial current of L1 is 1.0 A
    }
  }
  { Name Voltage_Cir ;
    Case {
      { Region VoltageSource1; Type Assign; Value 80.0; } // VoltageSource1 has 80.0 V
      { Region C1; Type Init; Value  0.0; } // Initial voltage of C1 = 0.0 V
    }
  }

  { Name ElectricalCircuit ; Type Network ;
    Case Circuit1 {
      // Circuit for cube 1:
      { Region VoltageSource1;  Branch{ 1,  10};}
      { Region R1;              Branch{ 1,  2};}
      // The voltage between nodes 2 and 3 corresponds to the absolute voltage
      // of the surface Top1:
      { Region Top1;            Branch{ 2,  3};}
      // Note the reverse order of the nodes in the next branch: it's for
      // subtracting the bottom voltage.
      { Region Bottom1;         Branch{ 10,  3};}
      // With the two lines above, the voltage between node 2 and 10 corresponds
      // to the voltage drop over the cube regardless what the absolute voltages
      // are.  The current between node 2 and 3 corresponds to the absolute
      // current flowing through Top1. The same is valid for Bottom1. Due to the
      // reverse order here the current has the opposite sign. This is correct
      // because the current flowing IN Top1 is flowing OUT of Bottom1. Note
      // that in this particular circuit it is actually not necessary to include
      // the bottom faces in the netlist, because all bottoms are kept at 0 V
      // via constraint Voltage_3D... But then all Bottom faces have to be
      // excluded also in the region Sur_Ele!

      // Circuit for cube 2:
      { Region L1;              Branch{ 4, 10};}
      { Region R2;              Branch{ 4, 10};}
      { Region Top2;            Branch{ 4,  5};}
      { Region Bottom2;         Branch{10,  5};}

      // Circuit for cube 3
      { Region CurrentSource1;  Branch{ 6, 10};}
      { Region C1;              Branch{ 6, 10};}
      { Region Top3;            Branch{ 6,  7};}
      { Region Bottom3;         Branch{10,  7};}
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
  // Function space for the FEM domain
  { Name Hgrad_V; Type Form0;
    BasisFunction {
       // All nodes but the grouped nodes of the surfaces for connecting the
       // circuitry
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Ele; Entity NodesOf[ All, Not Sur_Ele ]; }
      // Grouped nodes: Each surface of Sur_Ele has just one single voltage
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes;
        Support Dom_Ele; Entity GroupsOfNodesOf[ Sur_Ele ]; }
    }
    GlobalQuantity {
      { Name U; Type AliasOf; NameOfCoef vf; }
      { Name I; Type AssociatedWith;    NameOfCoef vf; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf ;         NameOfConstraint Voltage_3D; }
      { NameOfCoef U;  EntityType GroupsOfNodesOf ; NameOfConstraint Voltage_3D; }
      { NameOfCoef I;  EntityType GroupsOfNodesOf ; NameOfConstraint Current_3D; }
    }
  }

  // Function space for the circuit domain
  { Name Hregion_Cir; Type Scalar;
    BasisFunction {
      { Name sn; NameOfCoef ir; Function BF_Region;
        Support Domain_Cir; Entity Domain_Cir; }
    }
    GlobalQuantity {
      { Name Iz; Type AliasOf;          NameOfCoef ir; }
      { Name Uz; Type AssociatedWith;   NameOfCoef ir; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }

}


Formulation {
  { Name dynamic; Type FemEquation;
    Quantity {
      { Name V;  Type Local;    NameOfSpace Hgrad_V; }
      { Name U;  Type Global;   NameOfSpace Hgrad_V [U]; }
      { Name I;  Type Global;   NameOfSpace Hgrad_V [I]; }
      { Name Uz; Type Global;   NameOfSpace Hregion_Cir [Uz]; }
      { Name Iz; Type Global;   NameOfSpace Hregion_Cir [Iz]; }
    }
    Equation {
      // FEM domain
      Integral  { [ kappa[] * Dof{d V}, {d V} ]; In Vol_Ele;
        Integration I1; Jacobian JVol;  }
      GlobalTerm{ [ Dof{I},             {U} ];   In Sur_Ele; }

      // Circuit related terms
      // Resistance equation
      GlobalTerm{ [ Dof{Uz},            {Iz} ];  In Resistance_Cir; }
      GlobalTerm{ [ R[] * Dof{Iz},      {Iz} ];  In Resistance_Cir; }

      // Inductance equation
      GlobalTerm{ [ Dof{Uz},            {Iz} ];  In Inductance_Cir; }
      GlobalTerm{ DtDof[ L[] * Dof{Iz}, {Iz} ];  In Inductance_Cir; }

      // Capacitance equation
      GlobalTerm{ [ Dof{Iz},            {Iz} ]; In Capacitance_Cir; }
      GlobalTerm{ DtDof[ C[] * Dof{Uz}, {Iz} ]; In Capacitance_Cir; }

      // Inserting the network
      GlobalEquation { Type Network; NameOfConstraint ElectricalCircuit;
        { Node {I};     Loop {U};       Equation {I};   In Sur_Ele; }
        { Node {Iz};    Loop {Uz};      Equation {Uz};  In Domain_Cir; }
      }
    }
  }
}

Resolution {
  { Name dynamic;
    System {
      { Name A; NameOfFormulation dynamic; }
    }
    Operation {
      InitSolution[A];
      TimeLoopTheta[0.0, tStop, tStep, 1.0]{
        Generate[A];
        Solve[A];
        Test[SaveFct[]] { SaveSolution[A]; }
      }
    }
  }
}

PostProcessing {
  { Name Ele; NameOfFormulation dynamic;
    Quantity {
      { Name V;  Value{ Local{ [ {V} ];                 In Vol_Ele; Jacobian JVol;} } }
      { Name Jv; Value{ Local{ [ -kappa[]*{d V} ];      In Vol_Ele; Jacobian JVol;} } }
      { Name J;  Value{ Local{ [ Norm[kappa[]*{d V}] ]; In Vol_Ele; Jacobian JVol;} } }
      { Name U;  Value { Term { [ {U} ];                In Sur_Ele; } } }
      // The minus sign here is for getting positive currents at the top of the
      // cubes. This is because incomming currents are negative.
      { Name I;  Value { Term { [ -{I} ];               In Sur_Ele; } } }
    }
  }
}

PostOperation {
  // Absolute voltage everywhere [V]
  { Name V; NameOfPostProcessing Ele;
    Operation {
      Print[ V, OnElementsOf Vol_Ele, TimeLegend, File "Result_V.pos"];
    }
  }
  // Current through the top surfaces of the cubes [A]
  { Name I_Top; NameOfPostProcessing Ele;
    Operation {
      Print[ I, OnElementsOf Region[{Top1, Top2, Top3}], TimeLegend,
        File "Result_I_Top.pos"];
    }
  }
  // Current density vectors everywhere [A/m^2]
  { Name J_vectors; NameOfPostProcessing Ele;
    Operation {
      Print[ Jv, OnElementsOf Vol_Ele, TimeLegend, File "Result_J_vectors.pos"];
    }
  }
  // Magnitude of current density everywhere [A/m^2]
  { Name J_magnitude; NameOfPostProcessing Ele;
    Operation {
      Print[ J, OnElementsOf Vol_Ele, TimeLegend, File "Result_J.pos"];
    }
  }
  // Store the results in an ASCII file
  { Name Cube_Top_Values_ASCII; NameOfPostProcessing Ele;
    Operation {
      Print[U, OnRegion Top1, File "Result_Cube1TopValues.txt", Format TimeTable];
      Print[U, OnRegion Top2, File "Result_Cube2TopValues.txt", Format TimeTable];
      Print[U, OnRegion Top3, File "Result_Cube3TopValues.txt", Format TimeTable];

      Print[I, OnRegion Top1, File > "Result_Cube1TopValues.txt", Format TimeTable];
      Print[I, OnRegion Top2, File > "Result_Cube2TopValues.txt", Format TimeTable];
      Print[I, OnRegion Top3, File > "Result_Cube3TopValues.txt", Format TimeTable];
    }
  }

}
