from EXOSIMS.Prototypes.SurveySimulation import SurveySimulation
import astropy.units as u
import numpy as np
from ortools.linear_solver import pywraplp
from scipy.optimize import minimize,minimize_scalar
import os
try:
   import cPickle as pickle
except:
   import pickle
from astropy.time import Time
import pdb

class SLSQPScheduler(SurveySimulation):
    """SLSQPScheduler
    
    This class implements a continuous optimization of integration times
    using the scipy minimize function with method SLSQP.  ortools with the CBC 
    linear solver is used to find an initial solution consistent with the constraints.
    For details see Savransky et al. 2017 (SPIE).

    Args:         
        \*\*specs:
            user specified values

    Notes:
        Due to the time costs of the current comp_per_inttime calculation in GarrettCompleteness
        this should be used with BrownCompleteness.

        Requires ortools
    
    """

    def __init__(self, cacheOptTimes=False, staticOptTimes=False, selectionMetric='maxC', Izod='current',
        maxiter=60, ftol=1e-3, **specs): #fZminObs=False,
        
        #initialize the prototype survey
        SurveySimulation.__init__(self, **specs)

        #Calculate fZmax
        self.valfZmax, self.absTimefZmax = self.ZodiacalLight.calcfZmax(np.arange(self.TargetList.nStars), self.Observatory, self.TargetList,
            self.TimeKeeping, filter(lambda mode: mode['detectionMode'] == True, self.OpticalSystem.observingModes)[0], self.cachefname)

        assert isinstance(staticOptTimes, bool), 'staticOptTimes must be boolean.'
        self.staticOptTimes = staticOptTimes
        self._outspec['staticOptTimes'] = self.staticOptTimes

        assert isinstance(cacheOptTimes, bool), 'cacheOptTimes must be boolean.'
        self._outspec['cacheOptTimes'] = cacheOptTimes

        assert selectionMetric in ['maxC','Izod-Izodmin','Izod-Izodmax',
            '(Izod-Izodmin)/(Izodmax-Izodmin)',
            '(Izod-Izodmin)/(Izodmax-Izodmin)/CIzod', #(Izod-Izodmin)/(Izodmax-Izodmin)/CIzodmin is simply this but with Izod='fZmin'
            'TauIzod/CIzod', #TauIzodmin/CIzodmin is simply this but with Izod='fZmin'
            'random',
            'priorityObs'], 'selectionMetric not valid input' # Informs what selection metric to use
        self.selectionMetric = selectionMetric
        self._outspec['selectionMetric'] = self.selectionMetric

        assert Izod in ['fZmin','fZ0','fZmax','current'], 'Izod not valid input' # Informs what Izod to optimize integration times for [fZmin, fZmin+45d, fZ0, fZmax, current]
        self.Izod = Izod
        self._outspec['Izod'] = self.Izod

        assert isinstance(maxiter, int), 'maxiter is not an int' # maximum number of iterations to optimize integration times for
        assert maxiter >= 1, 'maxiter must be positive real'
        self.maxiter = maxiter
        self._outspec['maxiter'] = self.maxiter

        # assert isinstance(fZminObs, bool), 'fZminObs must be boolean' # True means Observations will occur at fZmin of targets
        # self.fZminObs = fZminObs

        assert isinstance(ftol, float), 'ftol must be boolean' # tolerance to place on optimization
        assert ftol > 0, 'ftol must be positive real'
        self.ftol = ftol
        self._outspec['ftol'] = self.ftol


        #some global defs
        self.detmode = filter(lambda mode: mode['detectionMode'] == True, self.OpticalSystem.observingModes)[0]
        self.ohTimeTot = self.Observatory.settlingTime + self.detmode['syst']['ohTime'] # total overhead time per observation
        self.maxTime = self.TimeKeeping.missionLife*self.TimeKeeping.missionPortion # total mission time

        self.constraints = {'type':'ineq',
                            'fun': lambda x: self.maxTime.to(u.d).value - np.sum(x[x*u.d > 0.1*u.s]) - #maxTime less sum of intTimes
                                             np.sum(x*u.d > 0.1*u.s).astype(float)*self.ohTimeTot.to(u.d).value, # sum of True -> goes to 1 x OHTime
                            'jac':lambda x: np.ones(len(x))*-1.}

        self.t0 = None
        if cacheOptTimes:
            #Generate cache Name########################################################################
            cachefname = self.cachefname + 't0'
            
            if os.path.isfile(cachefname):
                self.vprint("Loading cached t0 from %s"%cachefname)
                with open(cachefname, 'rb') as f:
                    self.t0 = pickle.load(f)
                sInds = np.arange(self.TargetList.nStars)
                fZ = np.array([self.ZodiacalLight.fZ0.value]*len(sInds))*self.ZodiacalLight.fZ0.unit
                self.scomp0 = -self.objfun(self.t0.to(u.d).value,sInds,fZ)


        if self.t0 is None:
            #1. find nominal background counts for all targets in list
            dMagint = 25.0 # this works fine for WFIRST
            self.vprint('dMagint: ' + str(self.dMagint))
            self.vprint('WAint: ' + str(self.WAint))
            _, Cbs, Csps = self.OpticalSystem.Cp_Cb_Csp(self.TargetList, range(self.TargetList.nStars),  
                    self.ZodiacalLight.fZ0, self.ZodiacalLight.fEZ0, dMagint, self.WAint, self.detmode)

            #find baseline solution with dMagLim-based integration times
            #DELETE self.vprint('Finding baseline fixed-time optimal target set.')
            #3.
            t0 = self.OpticalSystem.calc_intTime(self.TargetList, range(self.TargetList.nStars),  
                    self.ZodiacalLight.fZ0, self.ZodiacalLight.fEZ0, self.dMagint, self.WAint, self.detmode)
            #DELETE self.vprint('Calculated t0')
            #4.
            comp0 = self.Completeness.comp_per_intTime(t0, self.TargetList, range(self.TargetList.nStars), 
                    self.ZodiacalLight.fZ0, self.ZodiacalLight.fEZ0, self.WAint, self.detmode, C_b=Cbs, C_sp=Csps)
            
            #### 5. Formulating MIP to filter out stars we can't or don't want to reasonably observe
            #DELETE self.vprint('Instantiating Solver')
            solver = pywraplp.Solver('SolveIntegerProblem',pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING) # create solver instance
            xs = [ solver.IntVar(0.0,1.0, 'x'+str(j)) for j in range(len(comp0)) ] # define x_i variables for each star either 0 or 1
            self.vprint('Finding baseline fixed-time optimal target set.')

            #constraint is x_i*t_i < maxtime
            constraint = solver.Constraint(-solver.infinity(),self.maxTime.to(u.day).value) #hmmm I wonder if we could set this to 0,maxTime
            for j,x in enumerate(xs):
                constraint.SetCoefficient(x, t0[j].to(u.day).value + self.ohTimeTot.to(u.day).value) # this forms x_i*(t_0i+OH) for all i
            #DELETE self.vprint('Done defining constraint')

            #objective is max x_i*comp_i
            objective = solver.Objective()
            for j,x in enumerate(xs):
                objective.SetCoefficient(x, comp0[j])
            objective.SetMaximization()
            #DELETE self.vprint('Done defining objective')

            #solver.EnableOutput()# this line enables output of the CBC MIXED INTEGER PROGRAM (Was hard to find don't delete)
            solver.SetTimeLimit(5*60*1000)#time limit for solver in milliseconds
            cpres = solver.Solve() # actually solve MIP
            x0 = np.array([x.solution_value() for x in xs]) # convert output solutions

            self.scomp0 = np.sum(comp0*x0) # calculate sum Comp from MIP
            #DELETE self.vprint('number x0 set to 0: ' + str(len(x0[x0==0.])))
            self.t0 = t0 # assign calculated t0

            #Observation num x0=0 @ dMagint=25 is 1501
            #Observation num x0=0 @ dMagint=30 is 1501...
            #### Finished Running MIP for quick initial solution

            #now find the optimal eps baseline and use whichever gives you the highest starting completeness
            self.vprint('Finding baseline fixed-eps optimal target set.')
            def totCompfeps(eps):
                compstars,tstars,x = self.inttimesfeps(eps, Cbs.to('1/d').value, Csps.to('1/d').value)
                self.vprint('totCompfeps: ' + str(len(x[x==0.])) + ' len tstars: ' + str(len(tstars[tstars>1e-10])))
                return -np.sum(compstars*x)
            #print(saltyburrito)
            #Note: There is no way to seed an initial solution to minimize scalar 
            #0 and 1 are supposed to be the bounds on epsres. I could define upper bound to be 0.01, However defining the bounds to be 5 lets the solver converge
            epsres = minimize_scalar(totCompfeps,method='bounded',bounds=[0,7], options={'disp': 3, 'xatol':self.ftol, 'maxiter': self.maxiter})  #adding ftol for initial seed. could be different ftol
                #https://docs.scipy.org/doc/scipy/reference/optimize.minimize_scalar-bounded.html#optimize-minimize-scalar-bounded
            #DELETE self.vprint('Calculating inttimesfeps')
            comp_epsmax,t_epsmax,x_epsmax = self.inttimesfeps(epsres['x'],Cbs.to('1/d').value, Csps.to('1/d').value)
            if np.sum(comp_epsmax*x_epsmax) > self.scomp0:
                x0 = x_epsmax
                self.scomp0 = np.sum(comp_epsmax*x_epsmax) 
                self.t0 = t_epsmax*x_epsmax*u.day

            ##### Optimize the baseline solution
            self.vprint('Optimizing baseline integration times.')
            sInds = np.arange(self.TargetList.nStars)
            if self.Izod == 'fZ0': # Use fZ0 to calculate integration times
                fZ = np.array([self.ZodiacalLight.fZ0.value]*len(sInds))*self.ZodiacalLight.fZ0.unit
            elif self.Izod == 'fZmin': # Use fZmin to calculate integration times
                fZ = self.valfZmin
            elif self.Izod == 'fZmax': # Use fZmax to calculate integration times
                fZ = self.valfZmax
            elif self.Izod == 'current': # Use current fZ to calculate integration times
                fZ = self.ZodiacalLight.fZ(self.Observatory, self.TargetList, sInds, self.TimeKeeping.currentTimeAbs.copy()+np.zeros(self.TargetList.nStars)*u.d, self.detmode)

            maxIntTimeOBendTime, maxIntTimeExoplanetObsTime, maxIntTimeMissionLife = self.TimeKeeping.get_ObsDetectionMaxIntTime(self.Observatory, self.detmode, self.TimeKeeping.currentTimeNorm.copy())
            maxIntTime   = min(maxIntTimeOBendTime, maxIntTimeExoplanetObsTime, maxIntTimeMissionLife) # Maximum intTime allowed
            bounds = [(0,maxIntTime.to(u.d).value) for i in range(len(sInds))]
            #bounds = [(0,self.maxTime.to(u.d).value) for i in range(len(sInds))] # original method of defining bounds
            initguess = x0*self.t0.to(u.d).value
            self.save_initguess = initguess
            #DELETE self.vprint('x0')
            #DELETE self.vprint(x0)
            #DELETE self.vprint(len(x0[x0==1.]))
            #self.vprint([i for i.value in x0 if np.isinf(i)])
            #self.vprint([i for i in x0 if np.isnan(i)])
            #self.vprint([i for i in x0 if np.isinf(i)])
            #DELETE self.vprint('t0')
            #DELETE self.vprint(self.t0)
            #DELETE self.vprint([i for i in self.t0.value if np.isnan(i)])
            #DELETE self.vprint([i for i in self.t0.value if np.isinf(i)])
            #DELETE self.vprint([i for i in self.t0.value if i<0.])
            #DELETE self.vprint([i for i in self.t0.value if i==0.])
            #DELETE self.vprint('initguess: ' + str(initguess))
            #I specify ftol to be .5% of the scomp0 NEED TO IMPROVE THIS
            ires = minimize(self.objfun, initguess, jac=self.objfun_deriv, args=(sInds,fZ), 
                    constraints=self.constraints, method='SLSQP', bounds=bounds, options={'maxiter':self.maxiter, 'ftol':self.ftol, 'disp': True}) #original method
            #tol=self.scomp0*0.0001,

            assert ires['success'], "Initial time optimization failed."

            self.t0 = ires['x']*u.d
            self.scomp0 = -ires['fun']

            if cacheOptTimes:
                with open(cachefname,'wb') as f:
                    pickle.dump(self.t0, f)
                self.vprint("Saved cached optimized t0 to %s"%cachefname)

        #Redefine filter inds
        self.intTimeFilterInds = np.where((self.t0 > 0)*(self.t0 <= self.OpticalSystem.intCutoff) > 0)[0] # These indices are acceptable for use simulating    


    def inttimesfeps(self,eps,Cb,Csp):
        """
        Compute the optimal subset of targets for a given epsilon value
        where epsilon is the maximum completeness gradient.

        Everything is in units of days
        """

        tstars = (-Cb*eps*np.sqrt(np.log(10.)) + np.sqrt((Cb*eps)**2.*np.log(10.) + 
                   5.*Cb*Csp**2.*eps))/(2.0*Csp**2.*eps*np.log(10.)) # calculating Tau to achieve dC/dT #double check
        #if len(tstars[tstars<= 0.]) >= 1.:
        #    self.vprint('At least 1 tstar leq 0')
        #    print saltyburrito
        #tstars[tstars <= 0.] = self.maxTime.to(u.d).value # Setting to maxTime so cost of x being 1 is very high #NOT HOW THIS WORKS need to set to zero because some stars are returning negative
        #tstars[self.t0 == 0.] = self.maxTime.to(u.d).value # Setting to maxTime so cost of x being 1 is very high #added later to ensure previously filtered stars
        
        compstars = self.Completeness.comp_per_intTime(tstars*u.day, self.TargetList, 
                np.arange(self.TargetList.nStars), self.ZodiacalLight.fZ0, 
                self.ZodiacalLight.fEZ0, self.WAint, self.detmode, C_b=Cb/u.d, C_sp=Csp/u.d)

        
        solver = pywraplp.Solver('SolveIntegerProblem',pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
        xs = [ solver.IntVar(0.0,1.0, 'x'+str(j)) for j in range(len(compstars)) ]
        constraint = solver.Constraint(-solver.infinity(), self.maxTime.to(u.d).value)

        for j,x in enumerate(xs):
            constraint.SetCoefficient(x, tstars[j] + self.ohTimeTot.to(u.day).value)

        objective = solver.Objective()
        for j,x in enumerate(xs):
            objective.SetCoefficient(x, compstars[j])
        objective.SetMaximization()
        #solver.EnableOutput() # this line enables output of the CBC MIXED INTEGER PROGRAM (Was hard to find don't delete)
        solver.SetTimeLimit(5*60*1000)#time limit for solver in milliseconds


        cpres = solver.Solve()
        #self.vprint(solver.result_status())


        x = np.array([x.solution_value() for x in xs])
        self.vprint('Solver is FEASIBLE: ' + str(solver.FEASIBLE))
        self.vprint('Solver is OPTIMAL: ' + str(solver.OPTIMAL))
        self.vprint('Solver is BASIC: ' + str(solver.BASIC))
        #self.vprint('Solver is iterations: ' + str(solver.iterations))
        #self.vprint('Solver is nodes: ' + str(solver.nodes))

        #solver.reset()

        return compstars,tstars,x


    def objfun(self,t,sInds,fZ):
        """
        Objective Function for SLSQP minimization. Purpose is to maximize summed completeness

        Args:
            t (ndarray):
                Integration times in days. NB: NOT an astropy quantity.
            sInds (ndarray):
                Target star indices (of same size as t)
            fZ (astropy Quantity):
                Surface brightness of local zodiacal light in units of 1/arcsec2
                Same size as t

        """
        good = t*u.d >= 0.1*u.s # inds that were not downselected by initial MIP

        comp = self.Completeness.comp_per_intTime(t[good]*u.d, self.TargetList, sInds[good], fZ[good], 
                self.ZodiacalLight.fEZ0, self.WAint[sInds][good], self.detmode)
        self.vprint(-comp.sum())
        return -comp.sum()


    def objfun_deriv(self,t,sInds,fZ):
        """
        Jacobian of objective Function for SLSQP minimization. 

        Args:
            t (astropy Quantity):
                Integration times in days. NB: NOT an astropy quantity.
            sInds (ndarray):
                Target star indices (of same size as t)
            fZ (astropy Quantity):
                Surface brightness of local zodiacal light in units of 1/arcsec2
                Same size as t

        """
        good = t*u.d >= 0.1*u.s # inds that were not downselected by initial MIP

        tmp = self.Completeness.dcomp_dt(t[good]*u.d, self.TargetList, sInds[good], fZ[good], 
                self.ZodiacalLight.fEZ0, self.WAint[sInds][good], self.detmode).to("1/d").value

        jac = np.zeros(len(t))
        jac[good] = tmp
        return -jac



    def calc_targ_intTime(self, sInds, startTimes, mode):
        """
        Given a subset of targets, calculate their integration times given the
        start of observation time.

        This implementation updates the optimized times based on current conditions and 
        mission time left.

        Note: next_target filter will discard targets with zero integration times.
        
        Args:
            sInds (integer array):
                Indices of available targets
            startTimes (astropy quantity array):
                absolute start times of observations.  
                must be of the same size as sInds 
            mode (dict):
                Selected observing mode for detection

        Returns:
            intTimes (astropy Quantity array):
                Integration times for detection 
                same dimension as sInds
        """
 
        if self.staticOptTimes:
            intTimes = self.t0[sInds]
        else:
            # assumed values for detection
            if self.Izod == 'fZ0': # Use fZ0 to calculate integration times
                fZ = np.array([self.ZodiacalLight.fZ0.value]*len(sInds))*self.ZodiacalLight.fZ0.unit
            elif self.Izod == 'fZmin': # Use fZmin to calculate integration times
                fZ = self.valfZmin
            elif self.Izod == 'fZmax': # Use fZmax to calculate integration times
                fZ = self.valfZmax
            elif self.Izod == 'current': # Use current fZ to calculate integration times
                fZ = self.ZodiacalLight.fZ(self.Observatory, self.TargetList, sInds, startTimes, mode)

            #### instead of actual time left, try bounding by maxTime - detection time used
            #need to update time used in choose_next_target
            
            timeLeft = (self.TimeKeeping.missionLife.copy() - self.TimeKeeping.currentTimeNorm.copy())*self.TimeKeeping.missionPortion.copy()
            bounds = [(0,timeLeft.to(u.d).value) for i in range(len(sInds))]

            initguess = self.t0[sInds].to(u.d).value
            ires = minimize(self.objfun, initguess, jac=self.objfun_deriv, args=(sInds,fZ), constraints=self.constraints,
                    method='SLSQP', bounds=bounds, options={'disp':True,'maxiter':self.maxiter,'ftol':self.ftol})
            
            #update default times for these targets
            self.t0[sInds] = ires['x']*u.d

            intTimes = ires['x']*u.d
            
        intTimes[intTimes < 0.1*u.s] = 0.0*u.d
            
        return intTimes

    def choose_next_target(self, old_sInd, sInds, slewTimes, intTimes):
        """
        
        Given a subset of targets (pre-filtered by method next_target or some 
        other means), select the best next one. 

        Args:
            old_sInd (integer):
                Index of the previous target star
            sInds (integer array):
                Indices of available targets
            slewTimes (astropy quantity array):
                slew times to all stars (must be indexed by sInds)
            intTimes (astropy Quantity array):
                Integration times for detection in units of day
        
        Returns:
            sInd (integer):
                Index of next target star
            waitTime (astropy Quantity):
                the amount of time to wait (this method returns None)
        
        """
        #Do Checking to Ensure There are Targetswith Positive Nonzero Integration Time
        tmpsInds = sInds
        sInds = sInds[np.where(intTimes.value > 1e-10)]#filter out any intTimes that are essentially 0
        intTimes = intTimes[intTimes.value > 1e-10]
        #self.vprint(len(sInds))
        if len(sInds) == 0:#If there are no stars... arbitrarily assign 1 day for observation length otherwise this time would be wasted
            return None, None
            self.vprint('len sInds is 0')
            sInds = tmpsInds #revert to the saved sInds
            intTimes = (np.zeros(len(sInds)) + 1.)*u.d  

        # calcualte completeness values for current intTimes
        if self.Izod == 'fZ0': # Use fZ0 to calculate integration times
            fZ = np.array([self.ZodiacalLight.fZ0.value]*len(sInds))*self.ZodiacalLight.fZ0.unit
        elif self.Izod == 'fZmin': # Use fZmin to calculate integration times
            fZ = self.valfZmin[sInds]
        elif self.Izod == 'fZmax': # Use fZmax to calculate integration times
            fZ = self.valfZmax[sInds]
        elif self.Izod == 'current': # Use current fZ to calculate integration times
            fZ = self.ZodiacalLight.fZ(self.Observatory, self.TargetList, sInds,  
                self.TimeKeeping.currentTimeAbs.copy() + slewTimes[sInds], self.detmode)
        comps = self.Completeness.comp_per_intTime(intTimes, self.TargetList, sInds, fZ, 
                self.ZodiacalLight.fEZ0, self.WAint[sInds], self.detmode)

        #### Selection Metric Type
        valfZmax = self.valfZmax[sInds]
        valfZmin = self.valfZmin[sInds]
        if self.selectionMetric == 'maxC': #A choose target with maximum completeness
            sInd = np.random.choice(sInds[comps == max(comps)])
        elif self.selectionMetric == 'Izod-Izodmin': #B choose target closest to its fZmin
            selectInd = np.argmin(fZ - valfZmin)
            sInd = sInds[selectInd]
        elif self.selectionMetric == 'Izod-Izodmax': #C choose target furthest from fZmax
            selectInd = np.argmin(fZ - valfZmax)#this is most negative when fZ is smallest 
            sInd = sInds[selectInd]
        elif self.selectionMetric == '(Izod-Izodmin)/(Izodmax-Izodmin)': #D choose target closest to fZmin with largest fZmin-fZmax variation
            selectInd = np.argmin((fZ - valfZmin)/(valfZmin - valfZmax))#this is most negative when fZ is smallest 
            sInd = sInds[selectInd]
        elif self.selectionMetric == '(Izod-Izodmin)/(Izodmax-Izodmin)/CIzod': #E = D + current completeness at intTime optimized at 
            selectInd = np.argmin((fZ - valfZmin)/(valfZmin - valfZmax)*(1./comps))
            sInd = sInds[selectInd]
        #F is simply E but where comp is calculated sing fZmin
        # elif self.selectionMetric == '(Izod-Izodmin)/(Izodmax-Izodmin)/CIzodmin': #F = D + current completeness at Izodmin and intTime
        #     selectInd = np.argmin((fZ - valfZmin)/(valfZmin - valfZmax)*(1./comps))
        #     sInd = sInds[selectInd]
        elif self.selectionMetric == 'TauIzod/CIzod': #G maximum C/T
            selectInd = np.argmin(intTimes/comps)
            sInd = sInds[selectInd]
        elif self.selectionMetric == 'random': #I random selection of available
            sInd = np.random.choice(sInds)
        elif self.selectionMetric == 'priorityObs': # Advances time to 
            valfZmax = self.valfZmax[sInds].copy()
            valfZmin = self.valfZmin[sInds].copy()
            TK = self.TimeKeeping

            #Time relative to now where fZmin occurs
            timeWherefZminOccursRelativeToNow = self.absTimefZmin.copy().value - TK.currentTimeAbs.copy().value #of all targets
            indsLessThan0 = np.where((timeWherefZminOccursRelativeToNow < 0))[0] # find all inds that are less than 0
            cnt = 0.
            while len(indsLessThan0) > 0: #iterate until we know the next time in the future where fZmin occurs for all targets
                cnt += 1.
                timeWherefZminOccursRelativeToNow[indsLessThan0] = self.absTimefZmin.copy().value[indsLessThan0]\
                    - TK.currentTimeAbs.copy().value + cnt*365.25 #take original and add 365.25 until we get the right number of years to add
                indsLessThan0 = np.where((timeWherefZminOccursRelativeToNow < 0))[0]
            timeToStartfZmins = timeWherefZminOccursRelativeToNow#contains all "next occuring" fZmins in absolute time

            timefZminAfterNow = [timeToStartfZmins[i] for i in sInds]#range(len(timeToStartfZmins)) if i in sInds]#filter by times in future and times not filtered
            timeToAdvance = min(np.asarray(timefZminAfterNow))#find the minimum time
            #self.vprint(len(timeToStartfZmins))
            #self.vprint(len(sInds))
            sInd = np.where((timeToStartfZmins == timeToAdvance))[0][0]#find the index of the minimum time and return that sInd
            del timefZminAfterNow
            #if len(self.DRM) > 0:
            #    print 'prev Arr: ' + str(self.DRM[-1]['arrival_time'].value) + ' prev DetTime' + str(self.DRM[-1]['det_time'].value)
            #pdb.set_trace()
            #Advance To fZmin of Target
            success = self.TimeKeeping.advanceToAbsTime(Time(timeToAdvance+TK.currentTimeAbs.copy().value, format='mjd', scale='tai'), False)
            #print 'advc: ' + str(np.round(timeToAdvance,2)) + '  '
            #DELETE self.vprint(success)
            #self.vprint('advancedToAbsTime: ' + str(TK.currentTimeAbs.copy()))
            waitTime = None

            fZ = self.ZodiacalLight.fZ(self.Observatory, self.TargetList, sInds,  
                self.TimeKeeping.currentTimeAbs.copy() + slewTimes[sInds], self.detmode)
            selectInd = np.argmin(np.abs(fZ - valfZmin))#this is most negative when fZ is smallest 
            sInd = sInds[selectInd]

            # #Check if exoplanetObsTime would be exceeded
            # OS = self.OpticalSystem
            # Comp = self.Completeness
            # TL = self.TargetList
            # Obs = self.Observatory
            # TK = self.TimeKeeping
            # allModes = OS.observingModes
            # mode = filter(lambda mode: mode['detectionMode'] == True, allModes)[0]
            # maxIntTimeOBendTime, maxIntTimeExoplanetObsTime, maxIntTimeMissionLife = TK.get_ObsDetectionMaxIntTime(Obs, mode)
            # maxIntTime = min(maxIntTimeOBendTime, maxIntTimeExoplanetObsTime, maxIntTimeMissionLife)#Maximum intTime allowed
            # intTimes2 = self.calc_targ_intTime(sInd, TK.currentTimeAbs.copy(), mode)
            # if intTimes2 > maxIntTime: # check if max allowed integration time would be exceeded
            #     self.vprint('max allowed integration time would be exceeded')
            #     sInd = None
            #     waitTime = 1.*u.d
        #H is simply G but where comp and intTime are calculated using fZmin
        #elif self.selectionMetric == 'TauIzodmin/CIzodmin': #H maximum C at fZmin / T at fZmin
        #else:
        # if self.t0[sInd].value == 0:
        #     print saltyburrito


        if not sInd == None:
            if self.t0[sInd] < 1.0*u.s: # We assume any observation with integration time of less than 1 second is not a valid integration time
                sInd = None
        
        return sInd, None

    def arbitrary_time_advancement(self,dt):
        """ Handles fully dynamically scheduled case where OBduration is infinite and
        missionPortion is less than 1.
        Input dt is the total amount of time, including all overheads and extras
        used for the previous observation.
        """
        if self.selectionMetric == 'priorityObs':
            pass
        else:
            self.TimeKeeping.allocate_time( dt*(1 - self.TimeKeeping.missionPortion)/self.TimeKeeping.missionPortion,\
                addExoplanetObsTime=False )


