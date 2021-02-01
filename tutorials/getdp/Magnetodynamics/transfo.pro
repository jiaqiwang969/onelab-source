/* -------------------------------------------------------------------
   Tutorial 7b : magnetodyamic model of a single-phase transformer

   Features:
   - Use of a generic template formulation library
   - Circuit coupling used as a black-box (see Tutorial 8 for details)

   To compute the solution in a terminal:
       getdp transfo -solve Magnetodynamics2D_av -pos dyn

   To compute the solution interactively from the Gmsh GUI:
       File > Open > transfo.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

Include "transfo_common.pro";

DefineConstant[
  // The "Massive" option is non-physical, but it's interesting to highlight the
  // fact that both massive and stranded conductors can be handled in the same
  // way as far as circuit-coupling is concerned
  ConductorType = {2, Choices{1 = "Massive", 2 = "Coil"}, Highlight "Blue",
    Name "Parameters/01Conductor type"}
  Freq = {1, Min 0, Max 1e3, Step 1,
    Name "Parameters/Frequency"}
  Flag_FrequencyDomain = {0, Choices{0, 1},
    Name "Parameters/Frequency-domain?"}
  mur_Core = {1000, Min 1, Max 10000, Step 1,
    Name "Parameters/Core relative permeability"}
];

Group {
  // Physical regions:
  Air = Region[{AIR_WINDOW, AIR_EXT}];
  Sur_Air_Ext = Region[SUR_AIR_EXT]; // exterior boundary
  Core = Region[CORE]; // magnetic core of the transformer, assumed non-conducting
  Coil_1_P = Region[COIL_1_PLUS]; // 1st coil, positive side
  Coil_1_M = Region[COIL_1_MINUS]; // 1st coil, negative side
  Coil_1 = Region[{Coil_1_P, Coil_1_M}];
  Coil_2_P = Region[COIL_2_PLUS]; // 2nd coil, positive side
  Coil_2_M = Region[COIL_2_MINUS]; // 2nd coil, negative side
  Coil_2 = Region[{Coil_2_P, Coil_2_M}];
  Coils = Region[{Coil_1, Coil_2}];

  // Abstract regions that will be used in the "Lib_Magnetodynamics2D_av_Cir.pro"
  // template file included below;
  Vol_Mag = Region[{Air, Core, Coils}]; // full magnetic domain
  If (ConductorType == 1)
    Vol_C_Mag = Region[{Coils}]; // massive conductors
  ElseIf (ConductorType == 2)
    Vol_S_Mag = Region[{Coils}]; // stranded conductors (coils)
  EndIf
}

Function {
  mu0 = 4e-7*Pi;
  mu[Air] = 1 * mu0;
  mu[Core] = mur_Core * mu0;
  mu[Coils] = 1 * mu0;
  nu[] = 1 / mu[];

  sigma[Coils] = 1e7;

  // To be defined separately for each coil portion, to fix the convention of
  // positive current (1: along Oz, -1: along -Oz)
  SignBranch[Coil_1_P] = 1;
  SignBranch[Coil_1_M] = -1;
  SignBranch[Coil_2_P] = 1;
  SignBranch[Coil_2_M] = -1;

  If(ConductorType == 2)
    // Number of turns in the coils (same for PLUS and MINUS portions) - half
    // values because half coils are defined geometrically due to symmetry!
    Ns[Coil_1] = 100;
    Ns[Coil_2] = 200;

    // To be defined separately for each coil portion:
    Sc[Coil_1_P] = SurfaceArea[];
    Sc[Coil_1_M] = SurfaceArea[];
    Sc[Coil_2_P] = SurfaceArea[];
    Sc[Coil_2_M] = SurfaceArea[];

    // Current density in each coil portion for a unit current (will be
    // multiplied by the actual total current in the coil), in the case of
    // stranded conductors
    js0[Coils] = Ns[] / Sc[] * Vector[0, 0, SignBranch[]];
  EndIf

  // For a correct definition of the voltage
  CoefGeos[Coils] = SignBranch[] * thickness_Core;
}

// We will use a circuit coupling to connect the PLUS and MINUS portions of the
// coil in series, for both the primary and secondary. We will then apply a
// voltage source (with a small resistance in series) on the resulting branch on
// the primary, and connect a resistive load on the resulting branch on the
// secondary.
Flag_CircuitCoupling = 1;

// Note that the voltage will not be equally distributed in the PLUS and MINUS
// parts, which is the reason why we must apply the total voltage through a
// circuit -- and we cannot simply use a current source like in Tutorial 7a.

// Here is the definition of the circuits on primary and secondary sides:
Group {
  // Empty Groups to be filled
  Resistance_Cir  = Region[{}]; // all resistances
  Inductance_Cir  = Region[{}] ; // all inductances
  Capacitance_Cir = Region[{}] ; // all capacitances
  SourceV_Cir = Region[{}]; // all voltage sources
  SourceI_Cir = Region[{}]; // all current sources

  // Primary side
  E_in = Region[10001]; // arbitrary region number (not linked to the mesh)
  SourceV_Cir += Region[{E_in}];
  R_in = Region[10002]; // arbitrary region number (not linked to the mesh)
  Resistance_Cir += Region[{R_in}];

  // Secondary side
  R_out = Region[10101]; // arbitrary region number (not linked to the mesh)
  Resistance_Cir += Region[{R_out}];
}

Function {
  deg = Pi/180;
  // Input RMS voltage (half of the voltage because of symmetry; half coils are
  // defined!)
  val_E_in = 1.;
  phase_E_in = 90 *deg; // Phase in radian (from phase in degree)

  // High value for an open-circuit test; Low value for a short-circuit test;
  // any value in-between for any charge
  Resistance[R_out] = 10;

  // End-winding primary winding resistance for more realistic primary coil
  // model
  Resistance[R_in] = 1e-3;
}

Constraint {
  { Name MagneticVectorPotential_2D;
    Case {
      { Region Sur_Air_Ext; Value 0; }
    }
  }
  { Name Current_2D;
    Case {
    }
  }
  { Name Voltage_2D;
    Case {
    }
  }
  { Name Current_Cir ;
    Case {
    }
  }
  { Name Voltage_Cir ;
    Case {
      { Region E_in; Value val_E_in;
        // F_Cos_wt_p[] is a built-in function with two parameters (w and p),
        // that can be used to evaluate cos(w * t + p) in both frequency- and
        // time-domain
        TimeFunction F_Cos_wt_p[]{2*Pi*Freq, phase_E_in}; }
    }
  }
  { Name ElectricalCircuit ; Type Network ;
    Case Circuit_1 {
      // PLUS and MINUS coil portions to be connected in series, together with
      // E_in; an additional resistor is defined in series to represent the
      // Coil_1 end-winding, which is not considered in the 2D model.
      { Region E_in; Branch {1,4}; }
      { Region R_in; Branch {4,2}; }

      { Region Coil_1_P; Branch {2,3} ; }
      { Region Coil_1_M; Branch {3,1} ; }
    }
    Case Circuit_2 {
      // PLUS and MINUS coil portions to be connected in series, together with
      // R_out (an additional resistor could be defined to represent the Coil_2
      // end-winding - but we can directly add it to R_out as well)
      { Region R_out; Branch {1,2}; }

      { Region Coil_2_P; Branch {2,3} ; }
      { Region Coil_2_M; Branch {3,1} ; }
    }
  }
}

Include "Lib_Magnetodynamics2D_av_Cir.pro";

PostOperation {
  { Name dyn; NameOfPostProcessing Magnetodynamics2D_av;
    Operation {
      Print[ j, OnElementsOf Region[{Vol_C_Mag, Vol_S_Mag}], Format Gmsh, File "j.pos" ];
      Print[ b, OnElementsOf Vol_Mag, Format Gmsh, File "b.pos" ];
      Print[ az, OnElementsOf Vol_Mag, Format Gmsh, File "az.pos" ];
      If (Flag_FrequencyDomain)
        // In text file UI.txt: voltage and current of the primary coil (from
        // E_in) (real and imaginary parts!)
        Echo[ "E_in", Format Table, File "UI.txt" ];
        Print[ U, OnRegion E_in, Format FrequencyTable, File > "UI.txt" ];
        Print[ I, OnRegion E_in, Format FrequencyTable, File > "UI.txt"];

        // In text file UI.txt: voltage and current of the secondary coil (from
        // R_out)
        Echo[ "R_out", Format Table, File > "UI.txt" ];
        Print[ U, OnRegion R_out, Format FrequencyTable, File > "UI.txt" ];
        Print[ I, OnRegion R_out, Format FrequencyTable, File > "UI.txt"];
      EndIf
    }
  }
}
