/* -------------------------------------------------------------------
   Tutorial 2 : magnetostatic field of an electromagnet

   Features:
   - Infinite ring geometrical transformation
   - Parameters shared by Gmsh and GetDP, and ONELAB parameters
   - FunctionSpaces for the 2D vector potential formulation

   To compute the solution in a terminal:
       getdp electromagnet -solve MagSta_a
       getdp electromagnet -pos Map_a

   To compute the solution interactively from the Gmsh GUI:
       File > Open > electromagnet.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/* This model computes the static magnetic field produced by a DC current. This
   corresponds to a "magnetostatic" physical model, obtained by combining the
   time-invariant Maxwell-Ampere equation (curl h = js, with h the magnetic
   field and js the source current density) with Gauss' law (Div b = 0, with b
   the magnetic flux density) and the magnetic constitutive law (b = mu h, with
   mu the magnetic permeability).

   Since Div b = 0, b can be derived from a vector magnetic potential a, such
   that b = curl a. Plugging this potential in Maxwell-Ampere's law and using
   the constitutive law leads to a vector Poisson equation in terms of the
   magnetic vector potential: curl(nu curl a) = js, where nu = 1/mu is
   the reluctivity.

   Electromagnetic fields expand to infinity. The corresponding boundary
   condition can be imposed rigorously by means of a geometrical transformation
   that maps a ring (or shell) of finite elements to the complementary of its
   interior.  As this is a mere geometric transformation, it is enough in the
   model description to attribute a special Jacobian to the ring region
   ("AirInf") - see the "Jacobian" section below.  With this information, GetDP
   is able to correctly transform all quantities involved in the model.

   The special Jacobian "VolSphShell" takes 2 parameters in this case,
   "Val_Rint" and "Val_Rext", which represent the inner and outer radii of the
   transformed ring region and the value of which must match those used in the
   geometrical description of the model (.geo file).  This is a typical case
   where Gmsh and GetDP must consistently share parameter values.  To ensure
   consistency, common parameters are defined is a specific file
   "electromagnet_common.pro", which is included in both the .geo and .pro file
   of the model.

   Besides sharing parameters between Gmsh and GetDP, it is also useful to share
   some parameters (not all) with the user of the model, i.e., to make them
   editable in the GUI before running the model.  Such variables are called
   ONELAB variables (because the sharing mechanism between the model and the GUI
   uses the ONELAB interface).  ONELAB parameters are defined with a
   "DefineNumber" statement, which can be invoked in the .geo and .pro files. */

Group {
  // Physical regions:
  Air    = Region[ 101 ];   Core   = Region[ 102 ];
  Ind    = Region[ 103 ];   AirInf = Region[ 111 ];

  Surface_ht0 = Region[ 1100 ];
  Surface_bn0 = Region[ 1101 ];
  Surface_Inf = Region[ 1102 ];

  /* Abstract regions :
     The purpose of abstract regions is to allow a generic definition of the
     FunctionSpace, Formulation and PostProcessing fields with no reference to
     model-specific Physical regions.  We will show in a later tutorial how
     abstract formulations can then be isolated in geometry independent template
     files, thanks to an appropriate declaration mechanism (using
     DefineConstant[], DefineGroup[] and DefineFunction[]).

     The abstract regions in this model have the following interpretation:
     - Vol_Mag     : full volume domain
     - Vol_S_Mag   : region with imposed current source js
     - Vol_Inf_Mag : region where the infinite ring geometric transformation applies
     - Sur_Dir_Mag : part of the boundary with homogenous Dirichlet conditions
     - Sur_Neu_Mag : part of the boundary with non-homogeneous Neumann conditions
  */
  Vol_Mag     = Region[ {Air, AirInf, Core, Ind} ];
  Vol_S_Mag   = Region[ Ind ];
  Vol_Inf_Mag = Region[ AirInf ];
  Sur_Dir_Mag = Region[ {Surface_bn0, Surface_Inf} ];
  Sur_Neu_Mag = Region[ {} ]; // empty
}

/* The weak formulation for this problem is derived in a similar way to the
   electrostatic weak formulation from Tutorial 1. The main differences are that
   the fields are now vector-valued, and that we have one linear (source) term
   in addition to the bilinear term. The weak formulation reads: find a such
   that

   (curl(nu curl a), a')_Vol_Mag = (js, a')_Vol_S_Mag

   holds for all test functions a'. After integration by parts it reads: find a
   such that

   (nu curl a, curl a')_Vol_Mag + (n x (nu curl a), a')_Bnd_Vol_Mag = (js, a')_Vol_S_Mag

   In this electromagnet model the second (boundary) term vanishes, as there is
   either no test function a' (on the Dirichlet boundary), or "n x (nu curl a) =
   n x h" is zero (on the homogeneous Neumann boundary). We are thus eventually
   looking for functions a such that

   (nu curl a, curl a')_Vol_Mag = (js, a')_Vol_S_Mag

   holds for all a'. */

Function {
  mu0 = 4.e-7 * Pi;

  /* New ONELAB variables can then be defined using DefineNumber, e.g.: */
  murCore = DefineNumber[100, Name "Model parameters/Mur core",
			 Help "Magnetic relative permeability of Core"];

  /* When the script is run, if the parameter named "Model parameters/Mur core"
     has not been previously defined, it takes the value (100) provided in
     DefineNumber and is sent to the ONELAB server. The "/" character in the
     variable name is interpreted as a path separator, and results in the
     creation of a sub-tree in the graphical user interface. If the script is
     re-run later, the value will be updated using the value from the server
     (unless it is labeled as "ReadOnly").

     Various attributes can be specified for ONELAB parameters: see
     https://gitlab.onelab.info/doc/tutorials/wikis/ONELAB-syntax-for-Gmsh-and-GetDP */

  nu [ Region[{Air, Ind, AirInf}] ]  = 1. / mu0;
  nu [ Core ]  = 1. / (murCore * mu0);

  Current = DefineNumber[0.01, Name "Model parameters/Current",
			 Help "Current injected in coil [A]"];

  NbTurns = 1000 ; // number of turns in the coil
  js_fct[ Ind ] = -NbTurns*Current/SurfaceArea[];
  /* The minus sign is to have the current in -e_z direction,
     so that the magnetic induction field is in +e_y direction */
}

Constraint {
  { Name Dirichlet_a_Mag;
    Case {
      { Region Sur_Dir_Mag ; Value 0.; }
    }
  }
  { Name SourceCurrentDensityZ;
    Case {
      { Region Vol_S_Mag ; Value js_fct[]; }
    }
  }
}

/* In the 2D approximation, the magnetic vector potential a and the current
   density js are vectors with only the z-component, i.e.:

   a  = Vector [ 0, 0, az(x,y) ]
   js = Vector [ 0, 0, jsz(x,y) ]

   These vector fields behave differently under derivation and geometrical
   transformation though, and GetDP needs this information to perform these
   operations correctly. This is reflected in the Type, BasisFunction and Entity
   specified in the "Hcurl_a_Mag_2D" FunctionSpace for a ("Perpendicular 1-form"
   with "BF_PerpendicularEdge" basis functions associated to nodes of the mesh)
   and in the "Hregion_j_Mag_2D" FunctionSpace for js ("Vector" with
   "BF_RegionZ" basis functions associated with the region Vol_S_Mag). Without
   providing explicit details, a is thus node-based, whereas js is region-wise
   constant. */

Group {
  Dom_Hcurl_a_Mag_2D = Region[ {Vol_Mag, Sur_Neu_Mag} ];
}

FunctionSpace {
  { Name Hcurl_a_Mag_2D; Type Form1P; // Magnetic vector potential a
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Dom_Hcurl_a_Mag_2D ; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType NodesOf;
        NameOfConstraint Dirichlet_a_Mag; }
    }
  }

  { Name Hregion_j_Mag_2D; Type Vector; // Electric current density js
    BasisFunction {
      { Name sr; NameOfCoef jsr; Function BF_RegionZ;
        Support Vol_S_Mag; Entity Vol_S_Mag; }
    }
    Constraint {
      { NameOfCoef jsr; EntityType Region;
        NameOfConstraint SourceCurrentDensityZ; }
    }
  }

}

Include "electromagnet_common.pro";
Val_Rint = rInt; Val_Rext = rExt;

Jacobian {
  { Name Vol ;
    Case { { Region Vol_Inf_Mag ;
             Jacobian VolSphShell {Val_Rint, Val_Rext} ; }
           { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name Int ;
    Case { { Type Gauss ;
             Case { { GeoElement Triangle    ; NumberOfPoints  4 ; }
                    { GeoElement Quadrangle  ; NumberOfPoints  4 ; }
	}
      }
    }
  }
}

/* The function space "Hregion_j_Mag_2D" provides one basis function, and hence
   one degree of freedom, per physical region in the abstract region
   "Vol_S_Mag".  The constraint "SourceCurrentDensityZ" fixes all these dofs (a
   single one!), so the FunctionSpace "Hregion_j_Mag_2D" is fully fixed and has
   no FE unknowns.  One could thus have replaced it by a simple function and the
   Integral term below in the weak formulation could have been written as

   Integral { [ Vector[ 0,0,-js_fct[] ] , {a} ];
     In Vol_S_Mag; Jacobian Vol; Integration Int; }

   instead of

   Integral { [ - Dof{js} , {a} ];
     In Vol_S_Mag; Jacobian Vol; Integration Int; }

   The chosen implementation below is however more efficient as it avoids
   evaluating repeatedly the function js_fct[] during assembly.
*/

Formulation {
  { Name Magnetostatics_a_2D; Type FemEquation;
    Quantity {
      { Name a ; Type Local; NameOfSpace Hcurl_a_Mag_2D; }
      { Name js; Type Local; NameOfSpace Hregion_j_Mag_2D; }
    }
    Equation {
      // all terms on the left-hand side (hence the "-" sign in front of
      // Dof{js}):
      Integral { [ nu[] * Dof{d a} , {d a} ];
        In Vol_Mag; Jacobian Vol; Integration Int; }
      Integral { [ -Dof{js} , {a} ];
        In Vol_S_Mag; Jacobian Vol; Integration Int; }
    }
  }
}

Resolution {
  { Name MagSta_a;
    System {
      { Name Sys_Mag; NameOfFormulation Magnetostatics_a_2D; }
    }
    Operation {
      Generate[Sys_Mag]; Solve[Sys_Mag]; SaveSolution[Sys_Mag];
    }
  }
}

PostProcessing {
  { Name MagSta_a_2D; NameOfFormulation Magnetostatics_a_2D;
    Quantity {
      { Name a;
        Value {
          Term { [ {a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name az;
        Value {
          Term { [ CompZ[{a}] ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name b;
        Value {
          Term { [ {d a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name h;
        Value {
          Term { [ nu[] * {d a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name js;
        Value {
          Term { [ {js} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
    }
  }
}

e = 1.e-5;
h = 0.02;
p1 = {e,h,0};
p2 = {0.25-e,h,0}; // horizontal cut through model, just above x-axis.

PostOperation {

  { Name Map_a; NameOfPostProcessing MagSta_a_2D;
    Operation {
      Echo[ Str["l=PostProcessing.NbViews-1;",
		"View[l].IntervalsType = 1;",
		"View[l].NbIso = 40;"],
	    File "tmp.geo", LastTimeStepOnly] ;
      Print[ a, OnElementsOf Dom_Hcurl_a_Mag_2D, File "a.pos" ];
      Print[ js, OnElementsOf Dom_Hcurl_a_Mag_2D, File "js.pos" ];
      Print[ az, OnElementsOf Dom_Hcurl_a_Mag_2D, File "az.pos" ];
      Print[ b, OnElementsOf Dom_Hcurl_a_Mag_2D, File "b.pos" ];
      Print[ b, OnLine{{List[p1]}{List[p2]}} {50}, File "by.pos" ];
    }
  }
}
