"""
Purpose: To Plot C_0 vs T_0 and C_actual vs T_actual
Written by: Dean Keithly on 5/17/2018
"""
"""Example 1
I have 1000 pkl files in /home/dean/Documents/SIOSlab/Dean17Apr18RS01C01fZ01OB01PP01SU01/run146279583107.pkl and
1qty outspec file in /home/dean/Documents/SIOSlab/Dean17Apr18RS01C01fZ01OB01PP01SU01/outspec.json

To generate timelines for these run the following code from an ipython session
from ipython
%run PlotC0vsT0andCvsT.py '/home/dean/Documents/SIOSlab/Dean17Apr18RS01C01fZ01OB01PP01SU01/run146279583107.pkl' \
'/home/dean/Documents/SIOSlab/Dean17Apr18RS01C01fZ01OB01PP01SU01/outspec.json'
"""
"""Example 2
I have several folders with foldernames /home/dean/Documents/SIOSlab/*fZ*OB*PP*SU*/
each containing ~1000 pkl files and 1 outspec.json file

To plot a random Timeline from each folder, from ipython
%run PlotC0vsT0andCvsT.py '/home/dean/Documents/SIOSlab/' None
"""
#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/SIOSlab/Dean6May18RS09CXXfZ01OB09PP01SU01.json/run95764934358.pkl' '/home/dean/Documents/SIOSlab/Dean6May18RS09CXXfZ01OB09PP01SU01.json/outspec.json'
#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/SIOSlab/Dean6May18RS09CXXfZ01OB13PP01SU01/run295219944902.pkl' '/home/dean/Documents/SIOSlab/Dean6May18RS09CXXfZ01OB13PP01SU01/outspec.json'
#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/SIOSlab/Dean21May18RS09CXXfZ01OB01PP01SU01/run6012655441614.pkl' '/home/dean/Documents/SIOSlab/Dean21May18RS09CXXfZ01OB01PP01SU01/outspec.json'
#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean21May18RS09CXXfZ01OB01PP01SU01/run3492624809.pkl' '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean21May18RS09CXXfZ01OB01PP01SU01/outspec.json'
#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean19May18RS09CXXfZ01OB56PP01SU01/run1636735874.pkl' '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean19May18RS09CXXfZ01OB56PP01SU01/outspec.json'

#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean6June18RS09CXXfZ01OB56PP01SU01/run5442111239.pkl' '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean6June18RS09CXXfZ01OB56PP01SU01/outspec.json'
#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean6June18RS09CXXfZ01OB01PP01SU01/run7000640433.pkl' '/home/dean/Documents/exosims/EXOSIMS/EXOSIMS/Scripts/Dean6June18RS09CXXfZ01OB01PP01SU01/outspec.json'



#%run PlotC0vsT0andCvsT.py '/home/dean/Documents/SIOSlab/Dean6June18RS09CXXfZ01OB56PP01SU01/run254150360189.pkl' '/home/dean/Documents/SIOSlab/Dean6June18RS09CXXfZ01OB56PP01SU01/outspec.json'
#Dean6June18RS09CXXfZ01OB56PP01SU01.json  run245043802546.pkl
#Dean6June18RS09CXXfZ01OB01PP01SU01.json

try:
    import cPickle as pickle
except:
    import pickle
import os
if not 'DISPLAY' in os.environ.keys(): #Check environment for keys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt 
else:
    import matplotlib.pyplot as plt 
import numpy as np
from numpy import nan
import argparse
import json
import sys, os.path, EXOSIMS, EXOSIMS.MissionSim
import astropy.units as u
import copy
import random
import datetime
import re
from EXOSIMS.util.vprint import vprint
from scipy.optimize import minimize,minimize_scalar

class plotC0vsT0andCvsT(object):
    """Designed to plot Planned Completeness and Observed Completeness
    """
    _modtype = 'util'

    def __init__(self, args=None):
        """
        Args:
            args (dict) - 'file' keyword specifies specific pkl file to use
        """
        self.args = args
        pass

    def singleRunPostProcessing(self, PPoutpath, folder):
        """Generates a single yield histogram for the run_type
        Args:
            PPoutpath (string) - output path to place data in
            folder (string) - full filepath to folder containing runs
        """
        #Get name of pkl file
        if isinstance(self.args,dict):
            if 'file' in self.args.keys():
                file = self.args['file']
        else:
            file = self.pickPKL(folder)
        fullPathPKL = os.path.join(folder,file) # create full file path
        if not os.path.exists(fullPathPKL):
            raise ValueError('%s not found'%fullPathPKL)

        #Load pkl and outspec files
        try:
            with open(fullPathPKL, 'rb') as f:#load from cache
                DRM = pickle.load(f)
        except:
            vprint('Failed to open fullPathPKL %s'%fullPathPKL)
            pass
        outspecPath = os.path.join(folder,'outspec.json')
        try:
            with open(outspecPath, 'rb') as g:
                outspec = json.load(g)
        except:
            vprint('Failed to open outspecfile %s'%outspecPath)
            pass

        #Create Simulation Object
        sim = EXOSIMS.MissionSim.MissionSim(scriptfile=None, nopar=True, **outspec)
        SS = sim.SurveySimulation
        ZL = SS.ZodiacalLight
        COMP = SS.Completeness
        OS = SS.OpticalSystem
        Obs = SS.Observatory
        TL = SS.TargetList
        TK = SS.TimeKeeping

        #IF SurveySimulation module is SLSQPScheduler
        initt0 = None
        comp0 = None
        if 'SLSQPScheduler' in outspec['modules']['SurveySimulation']:
            #Extract Initial det_time and scomp0
            initt0 = sim.SurveySimulation.t0#These are the optmial times generated by SLSQP
            numObs0 = initt0[initt0.value>1e-10].shape[0]
            timeConservationCheck = numObs0*(outspec['settlingTime'] + outspec['starlightSuppressionSystems'][0]['ohTime'].value) + sum(initt0).value # This assumes a specific instrument for ohTime
            #assert abs(timeConservationCheck-outspec['missionLife']*outspec['missionPortion']*365.25) < 0.1, 'total instrument time not consistent with initial calculation'
            if not abs(timeConservationCheck-outspec['missionLife']*outspec['missionPortion']*365.25) < 0.1:
                vprint('total instrument time used is not within total allowed time with 0.1d')
            assert abs(timeConservationCheck-outspec['missionLife']*outspec['missionPortion']*365.25) < 0.5, 'total instrument time not consistent with initial calculation'
            #THIS IS JUST SUMCOMP initscomp0 = sim.SurveySimulation.scomp0

            _, Cbs, Csps = OS.Cp_Cb_Csp(TL, range(TL.nStars), ZL.fZ0, ZL.fEZ0, 25.0, SS.WAint, SS.detmode)

            #find baseline solution with dMagLim-based integration times
            #self.vprint('Finding baseline fixed-time optimal target set.')
            # t0 = OS.calc_intTime(TL, range(TL.nStars),  
            #         ZL.fZ0, ZL.fEZ0, SS.dMagint, SS.WAint, SS.detmode)
            comp0 = COMP.comp_per_intTime(initt0, TL, range(TL.nStars), 
                    ZL.fZ0, ZL.fEZ0, SS.WAint, SS.detmode, C_b=Cbs, C_sp=Csps)#Integration time at the initially calculated t0
            sumComp0 = sum(comp0)

            #Plot t0 vs c0
            plt.figure()
            plt.rc('axes',linewidth=2)
            plt.rc('lines',linewidth=2)
            #rcParams['axes.linewidth']=2
            plt.rc('font',weight='bold')
            #scatter(initt0.value, comp0, label='SLSQP $C_0$ ALL')
            plt.scatter(initt0[initt0.value > 1e-10].value, comp0[initt0.value > 1e-10], label=r'SLSQP $C_0$, $\sum C_0$' + "=%0.2f"%sumComp0, alpha=0.5, color='blue')


            #This is a calculation check to ensure the targets at less than 1e-10 d are trash
            sIndsLT1us = np.arange(TL.nStars)[initt0.value < 1e-10]
            t0LT1us = initt0[initt0.value < 1e-10].value + 0.1
            comp02 = COMP.comp_per_intTime(t0LT1us*u.d, TL, sIndsLT1us.tolist(), 
                    ZL.fZ0, ZL.fEZ0, SS.WAint[sIndsLT1us], SS.detmode, C_b=Cbs[sIndsLT1us], C_sp=Csps[sIndsLT1us])

            #Overwrite DRM with DRM just calculated
            res = sim.run_sim()
            DRM['DRM'] = sim.SurveySimulation.DRM



        #extract mission information from DRM
        arrival_times = [DRM['DRM'][i]['arrival_time'].value for i in np.arange(len(DRM['DRM']))]
        star_inds = [DRM['DRM'][i]['star_ind'] for i in np.arange(len(DRM['DRM']))]
        sumOHTIME = outspec['settlingTime'] + outspec['starlightSuppressionSystems'][0]['ohTime'].value
        raw_det_time = [DRM['DRM'][i]['det_time'].value for i in np.arange(len(DRM['DRM']))]#DOES NOT INCLUDE overhead time
        det_times = [DRM['DRM'][i]['det_time'].value+sumOHTIME for i in np.arange(len(DRM['DRM']))]#includes overhead time
        det_timesROUNDED = [round(DRM['DRM'][i]['det_time'].value+sumOHTIME,1) for i in np.arange(len(DRM['DRM']))]
        ObsNums = [DRM['DRM'][i]['ObsNum'] for i in np.arange(len(DRM['DRM']))]
        y_vals = np.zeros(len(det_times)).tolist()
        char_times = [DRM['DRM'][i]['char_time'].value*(1.+outspec['charMargin'])+sumOHTIME for i in np.arange(len(DRM['DRM']))]
        OBdurations = np.asarray(outspec['OBendTimes'])-np.asarray(outspec['OBstartTimes'])
        #sumOHTIME = [1 for i in np.arange(len(DRM['DRM']))]
        vprint(sum(det_times))
        vprint(sum(char_times))



        #calculate completeness at the time of each star observation
        slewTimes = np.zeros(len(star_inds))
        fZ = ZL.fZ(Obs, TL, star_inds, TK.missionStart + (arrival_times + slewTimes)*u.d, SS.detmode)
        comps = COMP.comp_per_intTime(raw_det_time*u.d, TL, star_inds, fZ, 
                ZL.fEZ0, SS.WAint[star_inds], SS.detmode)
        sumComps = sum(comps)


        if not plt.get_fignums(): # there is no figure open
            plt.figure()
        plt.rc('axes',linewidth=2)
        plt.rc('lines',linewidth=2)
        plt.rcParams['axes.linewidth']=2
        plt.rc('font',weight='bold')
        plt.scatter(raw_det_time, comps, label=r'SLSQP $C_{t_{Obs}}$, $\sum C_{t_{Obs}}$' + "=%0.2f"%sumComps, alpha=0.5, color='black')
        plt.xlim([0, 1.1*max(raw_det_time)])
        plt.ylim([0, 1.1*max(comps)])
        plt.xlabel(r'Integration Time, $\tau_i$, in (days)',weight='bold')
        plt.ylabel(r'Target Completeness, $C_i$',weight='bold')
        legend_properties = {'weight':'bold'}
        plt.legend(prop=legend_properties)

        #Done plotting Comp vs intTime of Observations
        date = unicode(datetime.datetime.now())
        date = ''.join(c + '_' for c in re.split('-|:| ',date)[0:-1])#Removes seconds from date
        fname = 'C0vsT0andCvsT_' + folder.split('/')[-1] + '_' + date
        plt.savefig(os.path.join(PPoutpath, fname + '.png'))
        plt.savefig(os.path.join(PPoutpath, fname + '.svg'))
        plt.savefig(os.path.join(PPoutpath, fname + '.eps'))


        #Manually Calculate the difference to veryify all det_times are the same
        tmpdiff = np.asarray(initt0[star_inds]) - np.asarray(raw_det_time)
        vprint(max(tmpdiff))

        vprint(-2.5*np.log10(ZL.fZ0.value)) # This is 23
        vprint(-2.5*np.log10(np.mean(fZ).value))

        self.writeDATAtoFile(initt0, numObs0, sumOHTIME, raw_det_time, PPoutpath, folder, date)
        self.plotCvsTlines(TL, Obs, TK, OS, SS, ZL, sim, COMP, PPoutpath, folder, date)

    def writeDATAtoFile(self, initt0, numObs0, sumOHTIME, raw_det_time, PPoutpath, folder, date):
        ############################################
        #### Calculate Lines for Data Output
        lines = []
        lines.append('Sum Planned Integration Time: ' + str(sum(initt0[initt0.value>1e-10])))
        lines.append('Number Planned Observations: ' + str(numObs0))
        lines.append('Planned Tsettling+Toh: ' + str(numObs0*sumOHTIME))
        RDT = [rdt for rdt in raw_det_time if rdt>1e-10]
        sumrdt = sum(RDT)
        lines.append('Sum Obs Integration Time: ' + str(sumrdt))
        lines.append('Number Obs Made: ' + str(len(RDT)))
        lines.append('Obs Tsettling+Toh: ' + str(len(RDT)*sumOHTIME))

        #### Save Data File
        fname = 'C0vsT0andCvsTDATA_' + folder.split('/')[-1] + '_' + date
        with open(os.path.join(PPoutpath, fname + '.txt'), 'w') as g:
            g.write("\n".join(lines))

    def plotCvsTlines(self, TL, Obs, TK, OS, SS, ZL, sim, COMP, PPoutpath, folder, date):
        """ Plots CvsT with Lines
        #From starkAYO_staticSchedule_withPlotting_copy_Feb6_2018.py
        #Lines 1246-1313, 1490-1502
        """
        sInds = np.arange(TL.nStars)
        mode = filter(lambda mode: mode['detectionMode'] == True, OS.observingModes)[0]
        fZ, fZabsTime = ZL.calcfZmin(sInds, Obs, TL, TK, mode, SS.cachefname)
        fEZ = ZL.fEZ0
        WA = OS.WA0
        dmag = np.linspace(1, COMP.dMagLim, num=1500,endpoint=True)
        Cp = np.zeros([sInds.shape[0],dmag.shape[0]])
        Cb = np.zeros(sInds.shape[0])
        Csp = np.zeros(sInds.shape[0])
        for i in xrange(dmag.shape[0]):
            Cp[:,i], Cb[:], Csp[:] = OS.Cp_Cb_Csp(TL, sInds, fZ, fEZ, dmag[i], WA, mode)
        Cb = Cb[:]#Cb[:,0]/u.s#note all Cb are the same for different dmags. They are just star dependent
        Csp = Csp[:]#Csp[:,0]/u.s#note all Csp are the same for different dmags. They are just star dependent
        #self.Cp = Cp[:,:] #This one is dependent upon dmag and each star
        
        TaylorCDFofGaussianFigsTau = plt.figure(703)
        cmap = plt.cm.get_cmap('autumn_r')
        intTimes = np.logspace(-6,3,num=400,base=10.0)#define integration times we will evaluate at
        actualComp = np.zeros([sInds.shape[0],intTimes.shape[0]])
        for j in np.arange(intTimes.shape[0]):
            actualComp[:,j] = COMP.comp_per_intTime((intTimes[j]+np.zeros([sInds.shape[0]]))*u.d, TL, sInds, fZ, fEZ, WA, mode, Cb/u.s, Csp/u.s)
        
        #Plot Top 10 black Lines
        compObs = COMP.comp_per_intTime(sim.SurveySimulation.t0, TL, sInds, fZ, fEZ, WA, mode, Cb/u.s, Csp/u.s)
        compObs2 = np.asarray([gg for gg in compObs if gg > 0.])
        tmpI = np.asarray([gg for gg in sInds if compObs[gg] > 0.])
        maxCI = np.argmax(compObs) # should return ind of max C0
        minCI = tmpI[np.argmin(compObs2)] # should return ind of min C0
        tmpI2 = np.argsort(compObs)[-10:]
        middleCI = compObs.tolist().index(np.percentile(compObs2,50,interpolation='nearest'))

        for l in np.arange(10):
            plt.plot(intTimes,actualComp[tmpI2[l],:],color='k',zorder=1)
        plt.plot(intTimes,actualComp[middleCI,:],color='k',zorder=1)
        plt.plot(intTimes,actualComp[minCI,:],color='k',zorder=1)
        ###############################

        plt.xscale('log')
        plt.rcParams['axes.linewidth']=2
        plt.rc('font',weight='bold') 
        plt.title('Generic Title I Forgot to Update',weight='bold')
        plt.xlabel(r'Integration Time, $\tau$ (days)',weight='bold',fontsize=14)
        plt.ylabel('Completeness',weight='bold',fontsize=14)
        plt.rc('axes',linewidth=2)
        plt.rc('lines',linewidth=2)
        #Plot Colorbar
        cmap = plt.cm.get_cmap('autumn_r')

        compatt0 = np.zeros([sInds.shape[0]])
        for j in np.arange(sInds.shape[0]):
            compatt0[j] = COMP.comp_per_intTime(sim.SurveySimulation.t0[j], TL, sInds[j], fZ[j], fEZ, WA, mode, Cb[j]/u.s, Csp[j]/u.s)
        plt.scatter(sim.SurveySimulation.t0,compatt0,color='k',marker='o',zorder=3,label=r'$C_{i}(\tau_{0})$')


        def plotSpecialPoints(ind, TL, OS, fZ, fEZ, COMP, WA, mode, sim):
            #### Plot Top Performer at dMagLim, max(C/t)
            tCp, tCb, tCsp = OS.Cp_Cb_Csp(TL, ind, fZ[ind], fEZ, COMP.dMagLim, WA, mode)
            Cdmaglim = COMP.comp_per_intTime(sim.SurveySimulation.t0[ind], TL, ind, fZ[ind], fEZ, WA, mode, tCb[0], tCsp[0])
            plt.scatter(sim.SurveySimulation.t0[ind],Cdmaglim,marker='x',color='red',zorder=3)

            def objfun(t, TL, tmpI, fZ, fEZ, WA, mode, OS):
                dmag = OS.calc_dMag_per_intTime(t*u.d, TL, tmpI, fZ, fEZ, WA, mode)#We must calculate a different dmag for each integraiton time
                Cp, Cb, Csp = OS.Cp_Cb_Csp(TL, tmpI, fZ, fEZ, dmag, WA, mode)#We must recalculate Cb and Csp at each dmag
                return -COMP.comp_per_intTime(t*u.d, TL, tmpI, fZ, fEZ, WA, mode, Cb, Csp)/t

            out = minimize_scalar(objfun,method='bounded',bounds=[0,10**3.], args=(TL, ind, fZ[ind], fEZ, WA, mode, OS))#, options={'disp': 3, 'xatol':self.ftol, 'maxiter': self.maxiter}) 
            tMaxCbyT = out['x']
            CtMaxCbyT = COMP.comp_per_intTime(tMaxCbyT*u.d, TL, ind, fZ[ind], fEZ, WA, mode, tCb[0], tCsp[0])
            plt.scatter(tMaxCbyT,CtMaxCbyT,marker='D',color='blue',zorder=3)
        plt.scatter(10**0,-1.,marker='x',color='red',zorder=3,label=r'$C_{\Delta mag_{lim}}$')
        plt.scatter(10**0,-1.,marker='D',color='blue',zorder=3,label=r'Max $C/\tau$')

        plotSpecialPoints(maxCI, TL, OS, fZ, fEZ, COMP, WA, mode, sim)
        plotSpecialPoints(middleCI, TL, OS, fZ, fEZ, COMP, WA, mode, sim)
        plotSpecialPoints(minCI, TL, OS, fZ, fEZ, COMP, WA, mode, sim)

        plt.plot([1e-5,1e-5],[0,0],color='k',label=r'Numerical $C_{i}(\tau)$',zorder=1)
        plt.legend(loc=2)
        plt.xlim([1e-6,10.*max(sim.SurveySimulation.t0.value)])
        plt.ylim([0,1.1*max(compatt0)])
        plt.show(block=False)

        fname = 'CvsTlines_' + folder.split('/')[-1] + '_' + date
        plt.savefig(os.path.join(PPoutpath, fname + '.png'))
        plt.savefig(os.path.join(PPoutpath, fname + '.svg'))
        plt.savefig(os.path.join(PPoutpath, fname + '.eps'))

    def multiRunPostProcessing(self, PPoutpath, folders):
        """Does Nothing
        Args:
            PPoutpath (string) - output path to place data in
            folders (string) - full filepaths to folders containing runs of each run_type
        """
        pass

    def pickPKL(self,folder):
        """Picks a PKL file from the provided folder
        """
        assert os.path.isdir(folder), 'The provided folder %s is not a folder'
        files = os.listdir(folder) # get files located in the provided folder
        assert len(files) > 0, 'There are no files in %s' %(folder)
        assert any('.pkl' in mystring for mystring in files), 'no files in folder are .pkl'
        return random.choice([file for file in files if '.pkl' in file])