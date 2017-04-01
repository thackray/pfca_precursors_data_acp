"""
Boxmodel built around a chembox unit.
Will read in config and inputs yaml files, and do a model run accordingly.
Will need to specify reactions, emissions, output locations, etc.

Interface:
boxmodel(configfile,inputsfile)
output = boxmodel.run(verbose=False)
"""
from datetime import timedelta
from chembox import initialize, ddt_emissions, dC_from_ddt
import reactions as rxnsmod
import numpy as np

def make_perturbation(old_value, u, x):
    return old_value + u*x

class BoxModel(object):
    
    """A Chembox-based single-box model.

    Public methods:
    run -- run model with options pre-defined and return requested output

    Attributes:
    chembox -- Chembox object with meta and input info for this box model
    """

    def __init__(self, configfilename, inputsfilename):
        """Construct BoxModel through Chembox object.

        Arguments:
        configfilename -- path to config file
        inputsfilename -- path to inputs file
        """

        self.chembox = initialize(configfilename, input_file=inputsfilename)
        return

    def rate_perturbations(self, pert):
        """Perturb the reaction rates according to the dict 'pert'."""
#        print self.chembox.meta.reactions
        for ind in pert:
#            print ind
            rxn = self.chembox.meta.reactions[ind]
            old_k, u = rxn['k'], rxn['u']
#            print old_k
            new_k = old_k
            new_k[0] = make_perturbation( old_k[0], u, pert[ind] )
            self.chembox.meta.reactions[ind]['k'] = new_k
            #            print new_k

    def cond_perturbations(self, pert):
        """Perturb the conditions according to the dict 'pert'."""
        for spc in pert:
            (dist, params), p = pert[spc]
            self.chembox.set_concentration(spc, 
                                           10**dist(params[0], params[1], p))


    def run(self,):
        """Run the BoxModel according to the already defined config options.

        Returns:
        dictionary with output concentrations in subdict under 'concs' and 
        output times under 't'
        """
    
        if self.chembox.inputs.emissions:
            dext = ddt_emissions(self.chembox)
            self.chembox.take_external(dext)
        for tout in self.chembox.meta.get_tout():
            dtime = tout-self.chembox.t
            dt = float(dtime.total_seconds())
            self.chembox.do_internal(dt)
        return self.chembox.give_output()

    def run_to_eq(self,species):
        """Run the BoxModel until the species of interest are steady."""

        self.fluxes = {'species':self.chembox.meta.get_species(),
                       'conspecies':self.chembox.meta.get_constant_species(),
                       'flux':[]}
        if self.chembox.inputs.emissions:
            dext = ddt_emissions(self.chembox)
            self.chembox.take_external(dext)
        stop = False
        for tout in self.chembox.meta.get_tout():
            if not stop:
                temp = True
                dtime = tout-self.chembox.t
                dt = float(dtime.total_seconds())
                self.chembox.do_internal(dt)
                netflux = getflux(self.chembox.passer,
                                  self.chembox.inputs.get_condition('T'),
                                  np.array([float(self.chembox.concs[spc])
                                            for spc in 
                                            self.chembox.meta.get_species()]))
                self.fluxes['flux'].append(np.abs(netflux))
                try:
                    for spc in species:
                        s1 = self.chembox.output[spc][-1]
                        s2 = self.chembox.output[spc][-2]
                        s3 = self.chembox.output[spc][-10]
                        temp *= abs(s1-s2)/s1 < 0.00001
                        temp *= abs(s1-s3)/s1 < 0.00001
                    temp *= self.chembox.output_t[-1]-self.chembox.output_t[0]>\
                            timedelta(days=10)
                    stop = temp
                except IndexError:
                    pass
        return self.chembox.give_output()

def getflux(p, T, cs):
#    print cs
    csi = np.where(np.isnan(cs),0.,cs)
    nu = p.stoich
    ks = p.ks
 #   print nu.shape
  #  print csi.shape
    k = np.array([rxnsmod.get_rate_coeff(ki,T) for ki in ks])

    prod = csi**np.where(nu<0,-1*nu,0.)
    #prod = csi**nu
  #  print prod.shape
    R = np.prod(prod,axis=1)
    dCs = np.dot(k*R, nu)

    return dCs


if __name__=='__main__':
    B = BoxModel('testing/config_BOXTEST.yaml', 'testing/inputs_BOXTEST.yaml')
    pert = {'1':1.0}
    B.rate_perturbations(pert)
    o = B.run()
    print(o['t'])
    for spc in o['concs']:
        print(o['concs'][spc])
        
