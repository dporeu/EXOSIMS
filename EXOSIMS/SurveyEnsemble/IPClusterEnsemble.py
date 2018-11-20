from __future__ import print_function

from ipyparallel import Client
from EXOSIMS.Prototypes.SurveyEnsemble import SurveyEnsemble 
from EXOSIMS.util.get_module import get_module
import time
from IPython.core.display import clear_output
import sys
import json
import os
import numpy as np
import EXOSIMS
import EXOSIMS.MissionSim
import os
import os.path
import cPickle
import random
import traceback
import subprocess


class IPClusterEnsemble(SurveyEnsemble):
    """Parallelized suvey ensemble based on IPython parallel (ipcluster)
    
    """

    def __init__(self, **specs):
        
        SurveyEnsemble.__init__(self, **specs)

        self.verb = specs.get('verbose', True)
        
        # access the cluster
        self.rc = Client()
        self.dview = self.rc[:]
        self.dview.block = True
        with self.dview.sync_imports(): import EXOSIMS, EXOSIMS.util.get_module, \
                os, os.path, time, random, cPickle, traceback
        if specs.has_key('logger'):
            specs.pop('logger')
        if specs.has_key('seed'):
            specs.pop('seed')
        self.dview.push(dict(specs=specs))
        res = self.dview.execute("SS = EXOSIMS.util.get_module.get_module(specs['modules'] \
                ['SurveySimulation'], 'SurveySimulation')(**specs)")

        res2 = self.dview.execute("SS.reset_sim()")

        self.vprint("Created SurveySimulation objects on %d engines."%len(self.rc.ids))
        #for row in res.stdout:
        #    self.vprint(row)

        self.lview = self.rc.load_balanced_view()

        self.maxNumEngines = len(self.rc.ids)

    def run_ensemble(self, sim, nb_run_sim, run_one=None, genNewPlanets=True,
        rewindPlanets=True, kwargs={}):
        """
        Args:
            sim:

        """
        t1 = time.time()
        async_res = []
        for j in range(nb_run_sim):
            ar = self.lview.apply_async(run_one, genNewPlanets=genNewPlanets,
                    rewindPlanets=rewindPlanets, **kwargs)
            async_res.append(ar)
        
        print("Submitted %d tasks."%len(async_res))
        

        runStartTime = time.time()#create job starting time
        avg_time_per_run = 0.
        tmplenoutstandingset = nb_run_sim
        tLastRunFinished = time.time()
        ar= self.rc._asyncresult_from_jobs(async_res)
        pids = ar.get_dict()
        while not ar.ready():
            ar.wait(10.)
            clear_output(wait=True)
            if ar.progress > 0:
                timeleft = ar.elapsed/ar.progress * (nb_run_sim - ar.progress)
                if timeleft > 3600.:
                    timeleftstr = "%2.2f hours"%(timeleft/3600.)
                elif timeleft > 60.:
                    timeleftstr = "%2.2f minutes"%(timeleft/60.)
                else:
                    timeleftstr = "%2.2f seconds"%timeleft
            else:
                timeleftstr = "who knows"

            #Terminate hanging runs
            outstandingset = self.rc.outstanding#a set of msg_ids that have been submitted but resunts have not been received
            if len(outstandingset) > 0 and len(outstandingset) < nb_run_sim:#there is at least 1 run still going and we have not just started
                avg_time_per_run = (time.time() - runStartTime)/float(nb_run_sim - len(outstandingset))#compute average amount of time per run
                if len(outstandingset) < tmplenoutstandingset:#The scheduler has finished a run
                    tmplenoutstandingset = len(outstandingset)#update this. should decrease by ~1 or number of cores...
                    tLastRunFinished = time.time()#update tLastRunFinished to the last time a simulation finished (right now)
                    #self.vprint("tmplenoutstandingset %d, tLastRunFinished %0.6f"%(tmplenoutstandingset,tLastRunFinished))
                if time.time() - tLastRunFinished > avg_time_per_run*(1. + self.maxNumEngines*2.)*4.:
                    #nb_run_sim = len(self.rc.outstanding)
                    #restartRuns = True
                    self.vprint('Aborting ' + str(len(self.rc.outstanding)) + 'qty outstandingset jobs')
                    runningPIDS = os.listdir('/proc') # get all running pids
                    for pid in pids:
                        if pid in runningPIDS:
                            os.kill(pid,9) # send kill command to stop this worker
                    stopIPClusterCommand = subprocess.Popen(['ipcluster','stop'])
                    time.sleep(20) # doing this instead of waiting for ipcluster to terminate
                    #stopIPClusterCommand.wait() # waits for process to terminate
                    #call(["ipcluster","stop"]) # send command to stop ipcluster
                    #self.rc.abort(jobs=self.rc.outstanding.copy().pop())
                    #self.rc.abort()#by default should abort all outstanding jobs... #it is possible that this will not stop the jobs running
                    #ar.wait(100)
                    #self.rc.purge_everything() # purge all results if outstanding *because rc.abort() didn't seem to do the job right
                    tLastRunFinished = time.time()#update tLastRunFinished to the last time a simulation was restarted (right now)

            print("%4i/%i tasks finished after %4i s. About %s to go." % (ar.progress, nb_run_sim, ar.elapsed, timeleftstr), end="")
            sys.stdout.flush()
        #numRunStarts += 1 # increment number of run restarts



        t2 = time.time()
        print("\nCompleted in %d sec" % (t2 - t1))
        
        res = [ar.get() for ar in async_res]
        
        return res
