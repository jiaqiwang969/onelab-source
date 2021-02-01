

// Define first the candidate design variables of the problem
// with appropriate bounds and a default first guess value. 
// Explain in the "Help" option their meaning. 
// These definitions must mandatorily be enclosed 
// in a DefineConstant[] statement.

DefineConstant[
 D = DefineNumber[0.01, Name "Optimization/Parameters/D", 
		  Min 0.001, Max 0.050, Help "Air gap length"]
 E = DefineNumber[0.1 , Name "Optimization/Parameters/E", 
		  Min 0.01 , Max 0.15, Help "C-core thickness"]
 F = DefineNumber[0.01, Name "Optimization/Parameters/F", 
		  Min 0.001, Max 0.090, Help "Core-coil gap"]
];

// "OnelabOptimize.py" needs to know the Onelab name of the design variables.
// They are communicated as the "Choices" of a Onelab string variable 
// with the reserved name "Optimization/Design parameters"
DesignParameters =
  DefineString["", Name "Optimization/Design parameters", Visible 0,
	       Choices {"Optimization/Parameters/D",
			"Optimization/Parameters/E",
			"Optimization/Parameters/F"}];

// "OnelabOptimize.py" can work in two modes:
// - parameter sweep and optimization over one individual design parameter
// - optimization over a multidimensional design space
// The mode is selected by the user with the parameter "SweepOnParameter".

// If "SweepOnParameter" is positive, the objective function and the sensibility 
// of the selected parameter are evaluated for "StepsInRange" points 
// distributed evenly over the parameter range.
// When completed, python prints in the terminal the results of the sweep
// with in particular a finite difference approximation of the sensitivity,
// in order to validate the sensitivity analysis for that parameter.
// The parameter sweep is followed by a minimization over the selected parameter,
// and skipped if "StepsInRange" is set to zero. 

// If "SweepOnParameter" is zero, the optimization is performed
// over the design space determined by the indices given in the "Choices"
// of the variable "DesignSpace". 
// The indices of the design parameter correspond to their order
// in the Choices list of "DesignParameters". 

SweepOnParameter = 
  DefineNumber[0, Name "Optimization/Sweep on parameter", Min 1, Max 3,
              Choices {0="Optimize on selected design space",
		       1="D",
		       2="E",
		       3="F"}];
StepsInRange = 
  DefineNumber[20, Name "Optimization/Steps in range", Visible SweepOnParameter];

DesignSpace =
  DefineNumber[0, Name "Optimization/Design space", Visible 0,
	       Choices {1,2}, Min 1, Max 3];

// available optimization algorithms
Optimizer = 
  DefineNumber[0, Name "Optimization/Algorithm",
               Choices {0="scipy", 
                        1="conveks-MMA",
                        2="conveks-GCMMA"}];
w = 
  DefineNumber[0, Name "Optimization/Results/w", 
	       Graph "030000000000000000000000000000000000"];


DefineConstant[
	     // all dimensions in meter
  L=1        // sidelength of the air box square

  CoreX=0.3  // X-position of C-core's bottom left corner
  CoreY=0.3  // Y-position of C-core's bottom left corner
  A=0.4      // C-core length along X-axis
  B=0.4      // C-core length along Y-axis
  D=0.01     // air gap length 
  E=0.1      // C-core thickness
  F=0.01     // gap between C-core and coil
  G=0.05     // coil width along X-axis

  CoilSection_ref = G*(B-2*E-2*F)

  R=1        // refinement factor: R=5 far raw mesh, R=1 normal mesh
];

By_target = 1.;
Px_target = CoreX + A - E/2 ;
Py_target = CoreY + B/2 ;
