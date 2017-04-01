"""
Chembox is the unit that does the single-box chemistry and physics to be
used in various applications.  
Interface from meta.py, and Euler/Lagrange/grid
and interface to outputs.py

Interface:
chembox.do_internal(delta_t)
chembox.take_external(delta_concs)
chembox.set_parameter(param, value)
chembox.initialize()
requested_output = chembox.give_output(args?)
"""

import os
import yaml
import operator
import numpy as np

from datetime import datetime, timedelta
from scipy.integrate import odeint

import reactions as rxnsmod
from meta import make_meta, make_inputs, load_yaml, Meta, Inputs
from error import ExternalForcingError


class Chembox(object):

    """An object to act as a box model or a base for larger undertakings.

    Public methods:
    do_internal -- use loaded modules to run chemistry
    take_external -- change concentrations based on method input
    scale_concentrations -- multiplicative scaling of concentrations
    give_output -- return model output concentrations
    give_meta_output -- return text about this chembox instance

    Attributes:
    meta -- Meta object containing info for this chembox instance
    inputs -- Inputs object containing info for this chembox instance
    concs -- dictionary of current concentrations for all species
    t -- current time in model world
    """

    def __init__(self, meta, inputs=None):
        """Construct Chembox according to meta (and input) info objects.

        Arguments:
        meta -- Meta object to inform Chembox initialization
        
        Keyword arguments:
        inputs -- Inputs object to further inform Chembox initialization
        """
        self.meta = meta
        self.inputs = inputs # T, p, emis, met, env types, etc.
        #unpack everything needed to init
        #make attr declarations
        self.concs = {}
        self.external = {}
        self.output = {}
        tstart = meta.get_tstart()
        self.t = tstart
        self.output_t = [tstart]
        self.output_diag = {'kR':[],'nu':[]}
        outspecies = self.meta.get_output_species()
        try:
            ic = self.inputs.get_initial_concentrations()
        except AttributeError:
            ic = {}
        for spc in self.meta.get_species():
            self.concs[spc] = ic.get(spc, 0.)
            if spc in outspecies:
                self.output[spc] = [self.concs[spc]]
        return

    def do_internal(self, delta_t):
        """Run the model for time step delta_t (seconds).

        Arguments:
        delta_t -- time step to advance chemistry by (float) in seconds
        """
        # do processes
        passer = Passer(self.meta,self.inputs)

        # new for matrix chem
        stoc = []
        ks = []
        RR = self.meta.reactions
        for rxn in RR:
            temp = []
            for spc in self.meta.get_species():
                if spc in RR[rxn]['R']:
                    temp.append(-1.)
                elif spc in RR[rxn]['P']:
                    temp.append(1.)
                else:
                    temp.append(0.)
            stoc.append(temp)
            ks.append(RR[rxn]['k'])
        stoc = np.array(stoc)
        passer.stoich = stoc
        passer.ks = ks
        # end of new

        assert isinstance(delta_t, float), 'delta_t must be float in sec'
        cs = [float(self.concs[spc]) for spc in self.meta.get_species()]
        kR_diag, nu_diag = ddt_chem_diag(cs, passer)
        self.output_diag['kR'].append(kR_diag)
        self.output_diag['nu'].append(nu_diag)
        wsol = odeint( ddt_total, cs, [0., delta_t], 
#                       args=({'meta':self.meta, 'inputs':self.inputs},) ) 
                       args=(passer,) )
        self.t += timedelta(seconds=delta_t)
        if True in [compare(self.t, to, eps=timedelta(seconds=delta_t/2.)) \
                    for to in self.meta.get_tout()]:
            self.output_t.append(self.t)
            for i, spc in enumerate(self.meta.get_species()):
                self.concs[spc] = wsol[-1,i]
                if spc in self.output:
                    self.output[spc].append(wsol[-1,i])        
        self.passer = passer
        return

    def take_external(self, delta_concs):
        """Use the given dictionary of delta_concs to add to model concs.

        Arguments:
        delta_concs -- dictionary of concentration values to add to Chembox
        """
        # unpack delta_conc dict
        # make changes to self.concs
        for spc in delta_concs:
            self.external[spc] = delta_concs[spc] + self.external.get(spc,0.)
        self.inputs.set_condition('external', self.external)
        return

    def scale_concentration(self, spc, scale_factor):
        """Scale concentration of given species by multiplicative factor.

        Arguments:
        spc -- the species to scale
        scale_factor -- the factor by which to scale conc of spc
        """
        self.concs[spc] = self.concs[spc]*scale_factor
        return

    def set_concentration(self, spc, new_value):
        """Set the concentration of spc to the given new value.

        Arguments:
        spc -- the species whose concentration to set
        new_value -- the new value at which to set the species concentration
        """
        self.concs[spc] = new_value
        return

    def give_output(self):
        """Give the requested outputs for this Chembox.

        Returns:
        dictionary with concentrations under 'concs' and times under 't'

        The output species and times are defined in Meta.
        """
        return {'concs':self.output, 't':self.output_t, 'meta':self.meta, 
                'diag':self.output_diag, }

    def give_meta_output(self):
        """Give the meta output.

        Returns:
        output string about this Chembox instance
        """
        out = ''
        out += self.meta.output()
        if self.inputs:
            out += self.inputs.output()
        return out


class Passer(object):
    def __init__(self, meta, inputs):
        self.meta = meta
        self.inputs = inputs
        assert isinstance(self.meta, Meta), "meta must be Meta object."
        assert isinstance(self.inputs, Inputs), "inputs must be Inputs object."
    def give(self):
        return self.meta, self.inputs


def mult(somelist):
    """Product of everything in somelist"""
    return reduce(operator.mul, somelist, 1)


def compare(x,y, eps=1e-6):
    """Return boolean of 'are x and y the same, within eps?'."""
    return abs(x-y) < eps


def ddt_total(cs, t, p):
    """Aggregate all the processes that are turned on, return dC array

    Arguments:
    cs -- array of concs aligned with species                                  
    t -- dummy time variable                                                   
    p -- Passer object                                          
                                                                             
    Returns:                                                                   
    array of dCs aligned with species 
    """

    dCs = np.zeros_like(cs)
    meta, inputs = p.give()
    if "DRY_DEPOSITION" in meta.process_list:
        dCs += ddt_drydep(cs, p)
    if "WET_DEPOSITION" in meta.process_list:
        dCs += ddt_wetdep(cs, p)
    if "CHEMISTRY" in meta.process_list:
        dCs += ddt_chem(cs, p)
    if "PHOTOLYSIS" in meta.process_list:
        print ("PHOTOLYSIS NO LONGER IMPLEMENTED")
    #        dCs += ddt_photolysis(cs, p)
    try:
        dCs += ddt_external(cs, p)
    except ExternalForcingError as e:
        pass
        
    return dCs
        

def ddt_external(cs, p):
    """Calculate change in concentrations due to external forcing.

    Arguments:
    cs -- array of concs aligned with species
    p -- Passer object

    Returns:
    array of dCs aligned with species
    """
    meta, inputs = p.give()
    conspecies = meta.get_constant_species()
    reactions = meta.get_reactions()
    species = meta.get_species()
    f = []
    try:
        external = inputs.get_condition('external')
    except KeyError as e:
        raise ExternalForcingError
    for i,spc in enumerate(species):
        dC = 0.
        if spc not in conspecies:
            dC += external.get(spc,0.)
        f.append(dC)
    return np.array(f)

            
def ddt_chem(cs, p):
    """Calculate change in concentrations due to chemistry.

    Arguments:                                                                 
    cs -- array of concs aligned with species                              
    p -- Passer object           

    Returns:                                                                   
    array of dCs aligned with species                                      
    """
    # To use with greater than 0-D, will have to flatten spatial array         
    """
    meta, inputs = p.give()
    conspecies = meta.get_constant_species()
    reactions = meta.get_reactions()
    species = meta.get_species()
    T = inputs.get_condition('T')
    f = []
    for i,spc in enumerate(species):
        #f.append(cs[i])                                                       
        dC = 0.
        if spc not in conspecies:
            for rxn in reactions:
                    dC += do_reaction(cs, reactions[rxn], spc, T, species)
        f.append(dC)
    return np.array(f)
    """
    csi = np.where(np.isnan(cs),0.,cs)
    meta, inputs = p.give()
    conspecies = meta.get_constant_species()
    reactions = meta.get_reactions()
    species = meta.get_species()
    T = inputs.get_condition('T')
    nu = p.stoich
    ks = p.ks
    k = np.array([rxnsmod.get_rate_coeff(ki,T) for ki in ks])

    prod = csi**np.where(nu<0,-1*nu,0.)
    R = np.prod(prod,axis=1)
    dCs = np.dot(k*R, nu)

    for i,spc in enumerate(species):
        if spc in conspecies:
            dCs[i] = 0.
    return dCs

def ddt_chem_diag(cs, p):
    """Calculate change in concentrations due to chemistry.

    Arguments:                                                                 
    cs -- array of concs aligned with species                              
    p -- Passer object           

    Returns:                                                                   
    array of dCs aligned with species                                      
    """
    # To use with greater than 0-D, will have to flatten spatial array         
    """
    meta, inputs = p.give()
    conspecies = meta.get_constant_species()
    reactions = meta.get_reactions()
    species = meta.get_species()
    T = inputs.get_condition('T')
    f = []
    for i,spc in enumerate(species):
        #f.append(cs[i])                                                       
        dC = 0.
        if spc not in conspecies:
            for rxn in reactions:
                    dC += do_reaction(cs, reactions[rxn], spc, T, species)
        f.append(dC)
    return np.array(f)
    """
    csi = np.where(np.isnan(cs),0.,cs)
    meta, inputs = p.give()
    conspecies = meta.get_constant_species()
    reactions = meta.get_reactions()
    species = meta.get_species()
    T = inputs.get_condition('T')
    nu = p.stoich
    ks = p.ks
    k = np.array([rxnsmod.get_rate_coeff(ki,T) for ki in ks])

    prod = csi**np.where(nu<0,-1*nu,0.)
    R = np.prod(prod,axis=1)
    return k*R, nu
                          

def do_reaction(cs, this_reaction, spc, T, species):
    """calculate dC resulting from this_reaction."""
    dC = 0.
    if this_reaction['k'][2]:
        if spc in this_reaction['R']:
            dC -= mult([cs[species.index(r)] for r in this_reaction['R']]) * \
                  rxnsmod.get_rate_coeff(this_reaction['k'], T)
        if spc in this_reaction['P']:
            dC += mult([cs[species.index(r)] for r in this_reaction['R']]) * \
                  rxnsmod.get_rate_coeff(this_reaction['k'], T)
    return dC


def ddt_drydep(cs, p):
    """Calculate dC due to dry deposition.

    Arguments:
    cs -- concentrations array aligned with species
    p -- Passer object

    Returns:
    array of dCs aligned with species
    """
    meta, inputs = p.give()
    try:
        delta_z = meta.get_size('z')
    except KeyError as e:
        print ("No delta z given, can't calculate drydep")
        raise e
    conspecies = meta.get_constant_species()
    species = meta.get_species()
    f = []
    for i,spc in enumerate(species):
        dC = 0.
        if spc not in conspecies:
            dC += -meta.get_parameter('vdep_%s'%spc)*cs[i]/delta_z
        f.append(dC)
    return np.array(f)


def ddt_wetdep(cs, p):
    """Calculate dC due to wet deposition.

    Arguments:
    cs -- concentrations vector aligned with species
    p -- Passer object

    Returns:
    array of dCs aligned with species
    """

    meta, inputs = p.give()
    Kaw = []
    species = meta.get_species()
    conspecies = meta.get_constant_species()
    for spc in species:
        Kaw.append(meta.get_parameter('Kaw_%s'%spc))
    Kaw = np.array(Kaw)
    R = inputs.get_condition('rainrate')
    zx = meta.get_size('z')
    E = 1.33
    z = inputs.get_condition('cloud_char_height')
    phi = inputs.get_condition('lw_air_frac')

    kwd = (R * E * np.exp(-z/zx)) / (zx * (Kaw+phi))    
    dC = - kwd * cs
    dC = np.where(compare(Kaw,-999),0.,dC)
    return dC


def ddt_photolysis(cs, p):
    """Calculate dC due to photolysis.

    Arguments:
    cs -- concentrations vector aligned with species
    p -- Passer object

    Returns:
    array of dCs aligned with species
    """
    meta, inputs = p.give()
    photo_yields = []
    for spc in meta.get_species():
        photo_yields.append(meta.get_parameter('photo_yield_%s'%spc))
    photo_yields = np.array(photo_yields)
    dC = -cs*photo_yields*inputs.get_condition('actinic_flux')
    return dC


def make_chembox(meta, inputs=None):
    """Wrapper for instantiating a Chembox from meta (and inputs, if given)."""
    chembox = Chembox(meta, inputs=inputs)
    return chembox


def initialize(config_file, input_file=None):
    """Instantiate a Chembox.

    If input_file not given, instantiates without inputs, will need
    to have them defined later.
    """

    config_options = load_yaml(config_file)
    meta = make_meta(config_options)
    if input_file:
        input_options = load_yaml(input_file)
        inputs = make_inputs(input_options)
        chembox = make_chembox(meta, inputs=inputs)
    else:
        chembox = make_chembox(meta)
    return chembox


def ddt_emissions(chembox):
    """Do timestep's emissions to chembox.

    Arguments:
    chembox -- Chembox object to be emitted into

    Returns:
    Chembox object with emissions done
    """

    dz = chembox.meta.get_size('z')
    dCs = {}
    species = chembox.meta.get_species()
    conspecies = chembox.meta.get_constant_species()
    for spc in species:
        if spc not in conspecies:
            E = chembox.inputs.get_emissions(spc)
            if E:
                M = chembox.meta.get_parameter('M',spc)
                dCs[spc] = E*6.022141e23 / (M*dz)
    return dCs


def dC_from_ddt(ddt, dt):
    """Converts rates to absolute changes over given time interval."""
    dC = {}
    for spc in ddt:
        dC[spc] = ddt[spc]*dt
    return dC


if __name__=="__main__":
    
    # python chembox.py at the command line will run this simple test suite.
    # These tests are meant to be inspectable by eye to determine that
    # the chembox interface is working properly, and to identify where
    # problems might be arising.  A more thorough automated suite of tests
    # will be added via other scripts.

    import testing


