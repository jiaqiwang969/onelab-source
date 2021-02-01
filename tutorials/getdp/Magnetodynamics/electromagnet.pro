/* -------------------------------------------------------------------
   Tutorial 7a : magnetic fields of an electromagnet

   Features:
   - Use of a template formulation library
   - Identical to Tutorial 2 for a static current source
   - Frequency-domain or time-domain solution for a time-dependent
     current source

   To compute the static solution in a terminal:
       getdp electromagnet -solve Magnetostatics2D_a -pos sta

   To compute the frequency-domain solution in a terminal:
       getdp electromagnet -solve Magnetodynamics2D_av -pos dyn

   To compute the time-dependent dynamic solution in a terminal:
       getdp electromagnet -setnumber TimeDomain 1 -solve Magnetodynamics2D_av -pos dyn

   To compute the solution interactively from the Gmsh GUI:
       File > Open > electromagnet.pro
       You may choose the Resolution in the left panel:
       Magnetodynamics2D_av (default) or Magnetostatics2D_a
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

Include "electromagnet_common.pro";

Group {
  // Physical regions
  Air    = Region[ 101 ];   Core   = Region[ 102 ];
  Ind    = Region[ 103 ];   AirInf = Region[ 111 ];
  Surface_ht0 = Region[ 1100 ];
  Surface_bn0 = Region[ 1101 ];
  Surface_Inf = Region[ 1102 ];

  // Abstract regions used in the "Lib_Magnetodynamics2D_av_Cir.pro" template file
  // that is included below:
  Vol_Mag = Region[{Air, Core, Ind, AirInf}]; // full magnetic domain
  Vol_C_Mag = Region[Core]; // massive conductors
  Vol_S_Mag = Region[Ind]; // stranded conductors (coils)
  Vol_Inf_Mag = Region[AirInf]; // ring-shaped shell for infinite transformation
  Val_Rint = rInt; Val_Rext = rExt; // interior and exterior radii of ring
}

Function {
  DefineConstant[
    murCore = {100, Name "Model parameters/Mur core"},
    Current = {0.01, Name "Model parameters/Current"},
    TimeDomain = {0, Choices{0 = "Frequency-domain", 1 = "Time-domain"},
      Name "Model parameters/03Analysis type"}
    frequency = {50, Visible !TimeDomain,
      Name "Model parameters/Frequency" }
  ];

  If(TimeDomain)
    // Fix parameters from the "Lib_Magnetodynamics2D_av_Cir.pro" template:
    Flag_FrequencyDomain = 0;
    TimeInit = 0; // start simulation at time = 0s
    TimeFinal = 20e-3; // stop simulation at time = 20 ms
    DeltaTime = 1e-3; // use time steps equal to 1 ms
    // Define the time modulation of the current source, i.e. a linear ramp from
    // 0 to 1 until 10 ms, then a constant value of 1:
    myModulation[] = ($Time < TimeFinal / 2) ? (2 / TimeFinal * $Time) : 1;
  Else
    // Fix parameters from the "Lib_Magnetodynamics2D_av_Cir.pro" template:
    Flag_FrequencyDomain = 1;
    Freq = frequency;
  EndIf

  mu0 = 4.e-7 * Pi;
  nu[ Region[{Air, Ind, AirInf}] ]  = 1. / mu0;
  nu[ Core ]  = 1. / (murCore * mu0);

  sigma[ Core ] = 1e6 / 10;
  sigma[ Ind ] = 5e7;

  Ns[ Ind ] = 1000 ; // number of turns in coil
  Sc[ Ind ] = SurfaceArea[] ; // area of coil cross section

  // Current density in each coil portion for a unit current (will be multiplied
  // by the actual total current in the coil)
  js0[ Ind ] = Ns[] / Sc[] * Vector[0, 0, -1];

  // For a correct definition of the voltage:
  CoefGeos[] = 1; // planar model, 1 meter thick
}

Constraint {
  { Name MagneticVectorPotential_2D;
    Case {
      { Region Surface_bn0; Value 0; }
      { Region Surface_Inf; Value 0; }
    }
  }
  { Name Current_2D;
    Case {
      If(Flag_FrequencyDomain)
        // Amplitude of the phasor is set to "Current"
        { Region Ind; Value Current; }
      Else
        // Time-dependent value is set to "Current * myModulation[]"
        { Region Ind; Value Current; TimeFunction myModulation[]; }
      EndIf
    }
  }
  { Name Voltage_2D;
    Case {
      { Region Core; Value 0; }
    }
  }
}

Include "Lib_Magnetodynamics2D_av_Cir.pro";

PostOperation {
  { Name dyn; NameOfPostProcessing Magnetodynamics2D_av;
    Operation {
      Print[ a, OnElementsOf Vol_Mag, File "a.pos" ];
      Print[ b, OnElementsOf Vol_Mag, File "b.pos" ];
      Print[ j, OnElementsOf Vol_Mag, File "j.pos" ];
    }
  }
  { Name sta; NameOfPostProcessing Magnetostatics2D_a;
    Operation {
      Print[ a, OnElementsOf Vol_Mag, File "a.pos" ];
      Print[ b, OnElementsOf Vol_Mag, File "b.pos" ];
      Print[ j, OnElementsOf Vol_Mag, File "j.pos" ];
    }
  }
}
