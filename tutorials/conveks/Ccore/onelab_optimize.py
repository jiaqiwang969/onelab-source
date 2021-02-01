import numpy as np
import conveks
import onelab
import scipy.optimize as spo
from shutil import copyfile

c = onelab.client(__file__)
x = {}

def setNodeTable(var, nodes, val):
    data = np.ravel(np.column_stack((nodes, val)))
    data = np.insert(data, 0, float(len(nodes)))
    c.setNumber(var, values=data, visible=0)

def getVelocityField(xvar, perturb=1e-7):
    for id, var in xvar.iteritems():
        c.runNonBlockingSubClient('myGmsh', mygmsh \
          + ' -setnumber PerturbMesh 1'\
          + ' -setnumber Perturbation '+str(perturb)\
          + ' -setnumber VelocityTag '+str(id)\
          + ' -setnumber '+var[-1]+' '+str(var[1]+perturb)\
          + ' -run')
    c.waitOnSubClients()
    # send the x,y,z components of each velocity field
    for label, var in x.iteritems():
        d = c.getNumbers('Optimization/Results/velocity_'+str(label))
        if label==0:
            nodes = np.array(d[1::4])
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_1',nodes,np.array(d[2::4]))
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_2',nodes,np.array(d[3::4]))
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_3',nodes,np.array(d[4::4]))
        c.setNumber(var[0],value=x[label][1])

def fobj(xFromOpti):
    # send the current point to GetDP model
    for label, var in x.iteritems():
        x[label][1] = xFromOpti[label]
        c.setNumber(var[0],value=x[label][1])

    # reload the geometry in GUI to see the geometrical changes
    c.openProject(file_geo)

    c.runSubClient('myGmsh', mygmsh + ' -2')
    getVelocityField(x)
    c.runSubClient('myGetDP', mygetdp + '-solve GetPerformances' \
                   + ' - setnumber OuterIteration ' + str(0)) # it !!
    return c.getNumber('Optimization/Results/w')


def grad(xFromOpti):
    # send the current point to GetDP model
    for label, var in x.iteritems():
        x[label][1] = xFromOpti[label]
        c.setNumber(var[0],value=x[label][1])
        
    # as well as their sensitivity with respect to design variables at `x'
    for dv, var in x.iteritems():
        c.runNonBlockingSubClient('myGetDP', mygetdp \
          + ' -setnumber VelocityTag '+str(dv)\
          + ' -solve GetGradient_wrt_dv')
    c.waitOnSubClients()
    grads=[]
    for dv, var in x.iteritems():
      grads.append(c.getNumber('Optimization/Results/dwdtau_'+str(dv)))
    return np.asarray(grads)

def grad_fd(xFromOpti, obj, basefile, perturb):
    # generate the perturbed meshes (for each design variable)
    getVelocityField(x,perturb) 
        
    # run the simulation for the perturbed mesh
    grads = []
    for dv, var in x.iteritems():
        copyfile(basefile+'Perturb_'+str(dv)+'.msh', basefile+'.msh')
        c.runSubClient('myGetDP', mygetdp + '-solve GetPerformances') 
        obj_perturb = c.getNumber('Optimization/Results/w')
        grads.append((obj_perturb-obj)/perturb)
    return np.asarray(grads)

def test(n):
    if n==0: return
    perturb = 1e-5
    tab = np.zeros((n+1,4))
    for step in range(n+1):
        xstep = x[0][2] + (x[0][3]-x[0][2])/float(n)*step
        f,g = fobj([xstep]), grad([xstep])[0]
        g_fd = grad_fd([xstep], f, modelName, perturb)[0]
        tab[step] = [xstep, f, g, g_fd]
        print "%4d" % step,
        for item in tab[step] : print "%3.6f" % item,
        print
        if c.getString(c.name+'/Action') == 'stop': exit(1)

def callbackF(Xi):
    print Xi
    if (c.getString(c.name+'/Action') == 'stop'):
        raise StopIteration 
        
def Optimize(Optimizer):
    print Optimizer
    print x
    if Optimizer==0: # 'scipy'
        values = [var[1] for i,var in x.iteritems()]
        bounds = [(var[2], var[3]) for i,var in x.iteritems()]
        print values, bounds
        try:
            res = spo.minimize(fobj, values, jac=grad, bounds=bounds,\
                        callback=callbackF,\
                        method = 'L-BFGS-B', tol=None,\
                        options={'ftol': 1e-6, 'gtol': 1e-3, 'disp':True} )
        except StopIteration:
            print 'Stopped by Gmsh'
            return
        print res
    elif Optimizer==1 or Optimizer==2: # 'mma'
        numVariables = len(x)
        initialPoint = np.zeros(numVariables)
        lowerBound = np.zeros(numVariables)
        upperBound = np.zeros(numVariables)
        for label, var in x.iteritems():
            x[label] = var
            initialPoint[label] = var[1]
            lowerBound[label] = var[2]
            upperBound[label] = var[3]

        # Initialize conveks
        conveks.initialize()

        # Initialize the MMA optimizer
        conveks.mma.initialize(initialPoint, lowerBound, upperBound)

        # Set some options for MMA
        if Optimizer == 1:
            conveks.mma.option.setNumber('General.Verbosity', 0) # def 4
            conveks.mma.option.setNumber('SubProblem.isContinuation', 1)
            conveks.mma.option.setNumber('SubProblem.cInitial'  , 1e2)
            conveks.mma.option.setNumber('SubProblem.cFinal'    , 1e5)
            conveks.mma.option.setNumber('SubProblem.cMultiply' , 10)
            conveks.mma.option.setNumber('SubProblem.cIteration', 50)
        
        # Get iteration count (here it will be 1 - could be different in case of restart)
        it = conveks.mma.getOuterIteration()
        innerit = 0; inneritTot = 0; change_mma = 1.0
        all_objectives = []
        print 'outer iter.   inner iter.   change   obj.   constr.'

        while it <= maxIter :
            c.setNumber('Optimization/01Current iteration', value=it, readOnly=1,
                         attributes={'Highlight':'LightYellow'})

            xFromMMA = conveks.mma.getCurrentPoint()

            if (c.getString(c.name+'/Action') == 'stop'): break
            obj = fobj(xFromMMA)
            if (c.getString(c.name+'/Action') == 'stop'): break
            grad_obj = grad(xFromMMA)
            constraints = np.array([np.sum(xFromMMA)/np.sum(upperBound)-1.0])
            grad_constraints = np.ones(numVariables)/np.sum(upperBound)
            all_objectives.append(obj)
            #print '%3i  %3i  %.4e %.4e  %.4e'%(it,innerit,change_mma,obj,constraints[0])
            #print xFromMMA

            if it == 1:
                fscale = 1.0
            #    fscale = 1.0 / np.abs(obj)
                print 'Initial: %4d %3.6f' %(it, obj),
            if it > 1:
                change_rel = np.abs(obj-obj_prec)*fscale
                print '%4d %3.6f %3.6f %3.6f' %(it, obj, change_rel, change_mma),
                if change_rel < maxChange: break
            print xFromMMA, grad_obj
            
            obj_prec = obj
            objective = obj*fscale
            grad_objective = grad_obj*fscale

            c.sendInfo('Optimization: it. {}, obj. {}'.format(it,objective))

            # call MMA and update the current point
            #conveks.mma.setMoveLimits(lowerBound, upperBound, 0.01) # parametre interessant
            #conveks.mma.updateCurrentPoint(grad_objective)
            conveks.mma.updateCurrentPoint(objective,grad_objective,constraints,grad_constraints)

            if Optimizer==2: # we use GCMMA, so update the approximation to be conservative
                # Evaluate the objective at the new point
                xFromMMA = conveks.mma.getCurrentPoint()
                objective_new = fobj(xFromMMA)
                constraints_new = np.array([np.sum(xFromMMA)/np.sum(upperBound)-1.0])

                # check if the approximation is conservative at the current point
                conserv = conveks.mma.isConservativeSubProblem(objective_new,constraints_new) 

                innerit = 0
                print '\t',innerit,xFromMMA,objective_new,objective
                while innerit <= 15 and conserv==0:
                    innerit += 1

                    # adapt the subproblem to make it conservative and generate a new point
                    conveks.mma.updateCurrentPoint(objective,grad_objective,constraints,grad_constraints,
                        objectiveAtUpdatedPoint=objective_new,constraintsAtUpdatedPoint=constraints_new)

                    # evaluate the objective and the constraints at the new point
                    xFromMMA = conveks.mma.getCurrentPoint()
                    objective_new = fobj(xFromMMA)
                    constraints_new = np.array([np.sum(xFromMMA)/np.sum(upperBound)-1.0])
                    print '\t',innerit,xFromMMA,objective_new,objective,constraints_new,constraints

                    # and check if the subproblem is conservative at this point
                    conserv = conveks.mma.isConservativeSubProblem(objective_new,constraints_new) 
                inneritTot += innerit

            change_mma = conveks.mma.getDesignChange()
            if change_mma < maxChange: break
            it = conveks.mma.getOuterIteration()

        print
        print 'Converged: %4d %3.6f %3.6f %3.6f' %(it, obj, change_rel, change_mma),
        print xFromMMA, grad_obj
        print 'objective'
        for i in all_objectives:print i
        
        # This should be called at the end of MMA
        conveks.mma.finalize()

        # This should be called at the end 
        conveks.finalize()

def minimize(name):
    global x, modelName, file_geo, file_pro, file_msh, mygmsh, mygetdp, maxIter, maxChange
    modelName = name
    # get gmsh and getdp locations from Gmsh options
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
    file_geo = c.getPath(modelName+'.geo')
    file_msh = c.getPath(modelName+'.msh')
    file_pro = c.getPath(modelName+'.pro')

    # build command strings with appropriate parameters and options
    mygmsh += ' -p 0 -v 2 -parametric ' + file_geo + ' '
    mygetdp += ' -p 0 ' + file_pro + ' -msh ' + file_msh + ' '


    # load geometry in GUI to have something to look at
    c.openProject(file_geo)

    c.sendCommand('Solver.AutoMesh = -1;')

    # dry getdp run (without -solve or -pos option) to get model parameters in the GUI
    c.runSubClient('myGetDP', mygetdp)

    # get design space from Onelab
    sweepOnParameter = int(c.getNumber('Optimization/Sweep on parameter'))
    if sweepOnParameter:
        designSpace = [ sweepOnParameter ]
    else:
        designSpace = c.getNumberChoices('Optimization/Design space')
    designParams = c.getStringChoices('Optimization/Design parameters')
        
    designParameters = {}
    iParam = 0
    for i in designSpace:
        label = designParams[int(i)-1]
        myrange = c.getNumberRange(label)
        param = [label]
        param.extend(myrange)
        param.extend([label.split('/')[-1]])
        designParameters[iParam]=param
        iParam +=1

    print
    print 'Working with design space of dimension', len(designParameters), ' :'
    for i, val in designParameters.iteritems():
        print i, val
    
    maxIter = c.defineNumber('Optimization/00Max iterations', value=30, visible=1)
    maxChange = c.defineNumber('Optimization/01Max change', value=1e-4,  visible=1)

    # end of check phase. We're done if we do not start the optimization
    if c.action == 'check' :
        exit(0)

    Optimizer = c.getNumber('Optimization/Algorithm')
    x=designParameters # x is global
    if sweepOnParameter:
        test(int(c.getNumber('Optimization/Steps in range')))
    Optimize(Optimizer)
    
    
## modelName = c.defineString('Optimization/Model name', value='')
## if modelName=='':
##     exit(1)

## minimize()
