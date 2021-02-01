# Open this file with Gmsh (interactively with File->Open->topo.py, or on the
# command line with 'gmsh topo.py')
#

import numpy as np
import onelab
import conveks

c = onelab.client(__file__)

# get paths to the gmsh and getdp executable from Gmsh options
mygmsh = c.getString('General.ExecutableFileName')
mygetdp = ''
for s in range(9):
   n = c.getString('Solver.Name' + str(s))
   if(n == 'GetDP'):
      mygetdp = c.getString('Solver.Executable' + str(s))
      break
if(not len(mygetdp)):
   c.sendError('This appears to be the first time you are trying to run GetDP')
   c.sendError('Please run a GetDP model interactively once with Gmsh to ' +
               'initialize the solver location')
   exit(0)

c.sendInfo('Will use gmsh={0} and getdp={1}'.format(mygmsh, mygetdp))

# append correct path to model file names
topo_geo = c.getPath('topo.geo')
topo_msh = c.getPath('topo.msh')
topo_pro = c.getPath('topo.pro')

# build command strings with appropriate parameters and options
mygetdp += ' -p 0 ' + topo_pro + ' -msh ' + topo_msh + ' '
mygmsh += ' -p 0 -v 2 ' + topo_geo + ' '

# load geometry in GUI to have something to look at
c.openProject(topo_geo)

# dry getdp run (without -solve or -pos option) to get model parameters in the GUI
c.runSubClient('myGetDP', mygetdp)

# get now optimization parameters
# some of them as Onelab parameter, to be editable in the GUI
Optimizer = c.getNumber('Optimization/00Optimizer')
maxIter = c.getNumber('Optimization/01Max iterations')
maxChange = c.getNumber('Optimization/02Max change')
filtering = c.getNumber('Optimization/3Density/1filter')

# end of check phase. We're done if we do not start the optimization
if c.action == 'check': exit(0)

# mesh the geometry to create the density (background) mesh
c.runSubClient('myGmsh', mygmsh + '-2')

# initializes the density 'xe' field
# This is done with a GetDP problem named 'DensityField'
# that solves no FE equations
# but only fills up the function space with an initial value.
c.runSubClient('myGetDP', mygetdp + '-solve DensityField')

# read in the 'xe' field from GetDP into a python variable 'd'
d = c.getNumbers('Optimization/Results/designVariable')

# the variable d contains the following serialized data: [N, ele1, xe1, ..., eleN, xeN]
# split it into homogeneous vectors
numElements = d[0] # N
elements = np.array(d[1::2]) # [ele1, ..., eleN ]
#x = np.array(d[2::2]) # [xe1, ..., xeN]
x = 0.5*np.ones(numElements)

print('number of elements in topology optimization:', numElements)

# This function serializes the python variables 'ele' and 'var' into
# an ElementTable variable [N, ele1, xe1, ..., eleN, xeN]
# and set it as a Onelab parameter named `var'
def setElementTable(var, ele, val):
   data = np.ravel(np.column_stack((ele, val)))
   data = np.insert(data, 0, numElements)
   c.setNumber(var, values=data, visible=0)

def get_objective(xFromOpti, flag_grad=1):
    if filtering:
        setElementTable('Optimization/Results/filterInput', elements, xFromOpti)
        c.runSubClient('myGetDP', mygetdp + '-solve FilterDensity')
        xFromOpti = np.array(c.getNumbers('Optimization/Results/filterOutput')[2::2])

    setElementTable('Optimization/Results/designVariable', elements, xFromOpti)

    c.runSubClient('myGetDP', mygetdp + '-solve ComputePerformancesAndGradient '\
        + '-setnumber Flag_getGradient '+str(flag_grad))

    # get the value of the objective function and of the constraints
    # as well as their sensitivity with respect to design variables at `x'
    objective = c.getNumber('Optimization/Results/pnorm of Von-Mises')
    constraints = np.array([c.getNumber('Optimization/Results/VolumeDomain')\
        / (0.5*c.getNumber('Optimization/Results/VolumeDomainFull')) - 1.0])

    return objective, constraints

def get_grad(xFromOpti):
    grad_objective = np.asarray(c.getNumbers('Optimization/Results/GradOfVonMises')[2::2])
    grad_constraints = np.asarray(c.getNumbers('Optimization/Results/GradOfVolume')[2::2])\
        /(0.5*c.getNumber('Optimization/Results/VolumeDomainFull'))
    if filtering:
        # apply the chain rule for the sensitivity of objective
        setElementTable('Optimization/Results/filterInput', elements, grad_objective)
        c.runSubClient('myGetDP', mygetdp + '-solve FilterDensity')
        grad_objective = np.array(c.getNumbers('Optimization/Results/filterOutput')[2::2])

        # apply the chain rule for the sensitivity of constraints
        d_con_dx = grad_constraints.reshape((1, numVariables))
        for k,dfe in enumerate(d_con_dx):
            setElementTable('Optimization/Results/filterInput', elements, dfe)
            c.runSubClient('myGetDP', mygetdp + '-solve FilterDensity')
            d_con_dx[k] = np.array(c.getNumbers('Optimization/Results/filterOutput')[2::2])
        grad_constraints = d_con_dx.reshape((1 * numVariables,))
    return grad_objective, grad_constraints

# number of design variables and number of constraints
numVariables = len(x)

# lower bound for design variables (here all set to 0.001)
lowerBound = 0.001*np.ones(numVariables)

# upper bound for design variables (here all set to 1)
upperBound = np.ones(numVariables)

# Initialize conveks
conveks.initialize()

# Initialize the MMA optimizer
conveks.mma.initialize(x, lowerBound, upperBound)

# Set appropriate options for MMA
conveks.mma.option.setNumber('General.Verbosity', 0)
# if Optimizer==0:conveks.mma.option.setNumber('SubProblem.move', 0.1)
# conveks.mma.option.setNumber('SubProblem.asymptotesRmax', 10)
# conveks.mma.option.setNumber('SubProblem.asymptotesRmin', 0.01)
# conveks.mma.option.setNumber('SubProblem.adaptSubproblem', 1)
if Optimizer != 1:
    conveks.mma.option.setNumber('SubProblem.isContinuation', 1)
    conveks.mma.option.setNumber('SubProblem.cInitial'  , 1e2)
    conveks.mma.option.setNumber('SubProblem.cFinal'    , 1e5)
    conveks.mma.option.setNumber('SubProblem.cMultiply' , 10)
    conveks.mma.option.setNumber('SubProblem.cIteration', 50)

# Get iteration count (here it will be 1 - could be different in case of restart)
it = conveks.mma.getOuterIteration()
innerit = 0; inneritTot = 0; change = 1.0; kkt_norm=1.0
all_objectives = []

# get (copy of) current point
x = conveks.mma.getCurrentPoint()

# compute objective function and constraints
# as well as their sensitivity with respect to design variables at `x'
objective,constraints = get_objective(x) 
grad_objective,grad_constraints = get_grad(x)

# send data to onelab server
c.setNumber('Optimization/Results/objective', value=objective) 
c.addNumberChoice('Optimization/Results/objective', value=objective) 
c.setNumber('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 
c.addNumberChoice('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 
c.setNumber('Optimization/03Inner Iteration', value=0) 
show_inner = 1

# show current iteration
# print 'iter.  inner iter.   obj.   max(constr.)   L2-norm(kkt)   point'
# print '%3i  %3i   %2.4e   %2.4e   %2.4e'%(it-1, innerit, objective, np.max(np.abs(constraints)), kkt_norm),
# print x[0:5]
print 'iter.  inner iter.   obj.   max(constr.)   L2-norm(kkt)   point'
print '%3i  %3i  %.6e  %.6e  %.6e  %3i'%(it-1, innerit, objective, np.max(np.abs(constraints)), kkt_norm, 0),
print '[',
for xk in x[0:5]:
    print '%.6e'%(float(xk)),
print ']'

while change > maxChange and it <= maxIter and c.getString('topo/Action') != 'stop':

    # Solve MMA
    conveks.mma.updateCurrentPoint(objective,grad_objective,constraints,grad_constraints)

    # and get the total number of nonlinear loops
    nlloops = conveks.mma.countNonlinearIterations()

    # Let us check now if the MMA approximation is conservative 
    # at the new point; and adapt the latter if it is not the case.
    if Optimizer == 1:
        # Evaluate the objective at the new point and check if 
        # the approximation is conservative at the current point
        obj_new,c_new = get_objective(conveks.mma.getCurrentPoint())
        conserv = conveks.mma.isConservativeSubProblem(obj_new, c_new) 

        innerit = 0
        # print '\t it.  new-point (xn)       f(xc)       f(xn)      max(c(xc))     max(c(xn))'
        # print '\t%3i'%(innerit),
        # print conveks.mma.getCurrentPoint()[0:5],
        # print '  %.4e  %.4e  %.4e  %.4e'%(objective,obj_new,np.max(constraints), np.max(c_new))
    
        while innerit <= 15 and conserv==0:
            # Do not show fields (density, von-Mises, ...) during inner iterations
            c.setNumber('Optimization/03Inner Iteration', value=1)
            innerit += 1

            # adapt the subproblem to make it conservative and generate a new point
            conveks.mma.updateCurrentPoint(objective,grad_objective,
                objectiveAtUpdatedPoint=obj_new,constraintsAtUpdatedPoint=c_new)

            # and get the total number of nonlinear loops
            nlloops = conveks.mma.countNonlinearIterations()

            # evaluate the objective and the constraints at the new point
            # and check if the subproblem is conservative at this point
            obj_new,c_new = get_objective(conveks.mma.getCurrentPoint())
            conserv = conveks.mma.isConservativeSubProblem(obj_new, c_new) 

            # show inner iteration
            # print '\t%3i'%(innerit),
            # print conveks.mma.getCurrentPoint()[0:5],
            # print '  %.4e  %.4e  %.4e  %.4e  %3i'%(objective,obj_new,np.max(constraints), np.max(c_new), nlloops)
            if show_inner:
                print '%3i  %3i '%(it, innerit),
                print '%.6e  %.6e  %.6e  %3i'%(obj_new,np.max(np.abs(c_new)), kkt_norm, nlloops),
                print '[',
                for xk in conveks.mma.getCurrentPoint()[0:5]:
                    print '%.6e'%(float(xk)),
                print ']'

        inneritTot += innerit
        c.setNumber('Optimization/03Inner Iteration', value=0) 

    # get (copy of) current point
    x = conveks.mma.getCurrentPoint()

    # compute objective function and constraints
    # as well as their sensitivity with respect to design variables at `x'
    objective,constraints = get_objective(x) 
    grad_objective,grad_constraints = get_grad(x)

    # evaluate the L2-norm of KKT 
    kkt_norm = conveks.mma.getKKTNorm(constraints,grad_objective,grad_constraints)
    
    # show the current state
    # print '%3i  %3i  %2.4e  %2.4e  %2.4e  %4i'%(it, innerit, objective, np.max(np.abs(constraints)), kkt_norm, nlloops),
    # print x[0:5]
    print '%3i  %3i  %.6e  %.6e  %.6e  %3i'%(it, innerit, objective, np.max(np.abs(constraints)), kkt_norm, nlloops),
    print '[',
    for xk in x[0:5]:
        print '%.6e'%(float(xk)),
    print ']'

    # get the outer iteration count here as well as the 
    # change between two successive points 
    it = conveks.mma.getOuterIteration()
    change = conveks.mma.getDesignChange()

    # send data to the onelab server
    c.setNumber('Optimization/Results/01Current iteration', value=it, readOnly=1)
    c.setNumber('Optimization/Results/02Current change', value=change, readOnly=1)
    c.setNumber('Optimization/Results/03KKT L2-norm', value=kkt_norm, readOnly=1)
    c.setNumber('Optimization/Results/objective', value=objective) 
    c.addNumberChoice('Optimization/Results/objective', value=objective) 
    c.setNumber('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 
    c.addNumberChoice('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 

# show the problem at the converged state 
print '*'*80
print 'Optimization has converged after %i iterations.'%(it-1)
print '  > point           :',x[0:5]
print '  > objective       : %.4e'%(objective)
print '  > max(constraints): %.4e'%(np.max(np.abs(constraints)))
print '  > L2-norm of kkt  : %.4e'%(kkt_norm)
print '  > inner iterations: %i'%(inneritTot)
print '*'*80

# This should be called at the end of MMA
conveks.mma.finalize()

# This should be called at the end
conveks.finalize()

print all_objectives