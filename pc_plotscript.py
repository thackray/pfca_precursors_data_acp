import numpy as np
import sys
import os
import pylab as pl
from datetime import datetime, timedelta

from boxmodel import BoxModel
import mypcm as pcm
from mypc import PCE
from pc_analysis import show_param_fracs, top_coeffs, get_vars, get_skew, gauss
from calculate_rates_yields import calculate_timescale as t_y
from pytools import pload, pdump, savefig

pfoa, pfna = 'C7F15C(O)OH','C8F17C(O)OH'
schain = 'C8F17O'
spc_interest = [pfoa, pfna, schain]

rundir = sys.argv[1]
n = 39; o = 2

if rundir.endswith('arctic/'):
    distlims = {'C8F17O':(0.001,1.0),
                'C8F17C(O)OH':(0.001,1.0),
                'C7F15C(O)OH':(0.001,1.0)}
elif rundir.endswith('urban/'):
    distlims = {'C8F17O':(0.001,1.0),
                'C8F17C(O)OH':(0.001,1),
                'C7F15C(O)OH':(0.001,1)}
elif rundir.endswith('ocean/'):
    distlims = {'C8F17O':(0.001,1.),
                'C8F17C(O)OH':(0.001,1.0),
                'C7F15C(O)OH':(0.001,1.0)}
    
    
for spc in spc_interest:

    Pt = pload(os.path.join(rundir, 'timing_poly_%s.pkl'%spc))
    Py = pload(os.path.join(rundir, 'yield_poly_%s.pkl'%spc))

    #print Pt.coeffs                                                           
    #    for lab, cof in zip( Py.labels, Py.coeffs):
    #    print lab, cof

    # yield distr.                                                             
    distyes = True
    if distyes:
 
        yield_unc = get_vars(Py)
        medy = Py.get_median()
        print medy
        stdy = yield_unc**0.5
        print stdy

        pl.figure(figsize=(16,9))
        lims = distlims[spc]
        x = np.linspace(np.log10(lims[0]),np.log10(lims[1]),1000)
        yyy = gauss(x, medy, stdy)
        pl.plot(10**x,yyy, lw=3)
        #pl.semilogx()
        pl.title('%s Yield Uncertainty'%spc, fontsize=30)
        savefig(os.path.join(rundir,'yield_dist_%s'%spc))
        pdump({'x':10**x,'y':yyy},os.path.join(rundir,'yield_dist_%s.pkl'%spc))
        # yield distr MC                                                            

        # timing distr.                                                             
        time_unc = get_vars(Pt)
        medt = Pt.get_median()
        stdt = time_unc**0.5

        pl.figure(figsize=(16,9))
        x = np.linspace(0,365,366)
        ttt = gauss(x, medt, stdt)
        pl.plot(x,ttt, lw=3)
        pl.title('%s Timing Uncertainty'%spc, fontsize=30)
        savefig(os.path.join(rundir,'timing_dist_%s'%spc))
        pdump({'x':x,'t':ttt},os.path.join(rundir,'timing_dist_%s.pkl'%spc))
    # timing distr MC                                                           

    fracyes = True
    if fracyes:
        # yield contribs.                                                           
        f = show_param_fracs(Py, range(n))
        pl.title('Yield %s'%spc, fontsize=20)
        savefig(os.path.join(rundir,'yield_frac_%s'%spc))
        pdump(f,os.path.join(rundir,'yield_frac_%s.pkl'%spc))
        # timing contribs.                                                          
        f = show_param_fracs(Pt, range(n))
        pl.title('Timing %s'%spc, fontsize=20)
        savefig(os.path.join(rundir,'timing_frac_%s'%spc))
        pdump(f,os.path.join(rundir,'timing_frac_%s.pkl'%spc))

    specyes = True
    if specyes:
        # spectrum of coefs                                                         
        top_coeffs(Py, ncoeffs = 15)
        pl.title('Yield %s'%spc, fontsize=20)
        savefig(os.path.join(rundir,'yield_coeffs_%s'%spc))

        top_coeffs(Pt, ncoeffs = 15)
        pl.title('Timing %s'%spc, fontsize=20)
        savefig(os.path.join(rundir,'timing_coeffs_%s'%spc))

pl.show()


