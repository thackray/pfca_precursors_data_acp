"""module defining Meta and Inputs objects used to specify Chembox instances.
"""

import os
import yaml
from datetime import datetime, timedelta

import reactions as rxnsmod

from error import YamlBooleanError

VERSION = 'Chembox version 0.1'
DEFAULTSTARTTIME = datetime(1987, 6, 7)

def add_output_block(string, name, listy, breakcount = 10):
    """Add a section of output to the output string.

    Arguments:
    string -- existing output string to be added to
    name -- name for this block of output
    listy -- list or dictionary to be summarized in output

    Keyword arguments:
    breakcount -- the number of items from listy to display on a line

    returns:
    string with new output block added to it
    """
    string += name + os.linesep
    for count, elem in enumerate(listy):
        if not count%breakcount: string += os.linesep
        if isinstance(listy, dict):
            string += str(elem) + ': ' + str(listy[elem]) + ' , '
        else:
            string += str(elem) + ' , '
    string += os.linesep + '-'*25 + os.linesep
    return string


class Meta(object):
    
    """Container object for chembox metadata. 

    Model-level metadata: 
    parameters -- dictionary of parameter values
    process_list -- list of process to model (dep, chem, etc.)
    reactions -- reactions dictionary created my reactions.py
    tout -- times to collect output at
    tstart -- reference start time
    constant_species -- species to hold constant in run
    output_species -- species to collect output for
    species -- list of all species in model (automatically generated)
    size -- spatial size parameters of box (height, horizontal distances, etc.)

    Public methods:
    set_size
    get_size
    set_parameter
    get_parameter
    add_process
    remove_process
    set_reactions
    get_reactions
    set_tout
    get_tout
    set_tstart
    get_tstart
    add_constant_species
    get_constant_species
    add_output_species
    get_output_species
    add_species
    get_species
    output
    """

    def __init__(self,):
        self.parameters = {}
        self.process_list = []
        self.reactions = None
        self.tout = None
        self.tstart = None
        self.constant_species = []
        self.output_species = []
        self.species = []
        self.size = {}
        return

    def set_size(self, dimension, value):
        self.size[dimension] = value
        return

    def get_size(self, dimension):
        return self.size[dimension]

    #def calculate_sizey_things(self, thing):
        # calculate area or whatever
        #self.size[thing] = calculated

    def set_parameter(self, name, value):
        self.parameters[name] = value
        return

    def get_parameter(self, name, spc=None):
        if spc:
            return self.parameters['%s_%s'%(name,spc)]
        else:
            return self.parameters[name]

    def add_process(self, process):
        self.process_list.append(process)
        return

    def remove_process(self, process):
        self.process_list.pop(process_list.index(process))
        return

    def set_reactions(self, reactions_dict):
        for rxn in reactions_dict:
            reactions_dict[rxn]['k'][0] = float(reactions_dict[rxn]['k'][0])
        self.reactions = reactions_dict
        return

    def get_reactions(self,):
        return self.reactions

    def set_tout(self, tout):
        self.tout = []
        assert isinstance(tout, list), 'OUTPUT TIME MUST BE LIST'
        for ti in tout:
            if isinstance(ti, str):
                try:
                    self.tout.append(datetime.strptime(ti, '%Y %m %d %H:%M'))
                except ValueError:
                    self.tout.append(datetime.strptime(ti, '%Y %m %d'))
            if isinstance(ti, float):
                self.tout.append(DEFAULTSTARTTIME+timedelta(seconds=ti))
            if isinstance(ti, datetime):
                self.tout.append(ti)
        return

    def get_tout(self, ):
        return self.tout

    def set_tstart(self, t):
        if not isinstance(t, datetime):
            t = DEFAULTSTARTTIME+timedelta(seconds=t)
        self.tstart = t

    def get_tstart(self):
        return self.tstart

    def add_constant_species(self, spc):
        self.constant_species.append(spc)
        return

    def get_constant_species(self,):
        return self.constant_species

    def add_output_species(self, spc):
        # make sure not to duplicate
        self.output_species.append(spc)
        return

    def get_output_species(self,):
        return self.output_species

    def add_species(self, spc):
        self.species.append(spc)
        return

    def get_species(self,):
        return self.species

    def output(self):
        """Return metadata output for informative purposes"""
        so = ''
        so += VERSION + os.linesep
        so += str(datetime.now()) + os.linesep
        so += '-'*25 + os.linesep

        so = add_output_block(so, 'SPECIES IN MODEL', self.species)
        so = add_output_block(so, 'SPECIES HELD CONSTANT', 
                              self.constant_species)
        so = add_output_block(so, 'SPECIES TO OUTPUT', 
                              self.output_species)
        so = add_output_block(so, 'INCLUDED PROCESSES', self.process_list)
        so = add_output_block(so, 'MODELED REACTIONS', self.reactions)
        so = add_output_block(so, 'OUTPUT TIMES', self.tout)
        so = add_output_block(so, 'MODEL PARAMETERS', self.parameters)
        return so


class Inputs(object):

    """Container object for model inputs.  

    Inputs herein:
    conditions -- dictionary of atmospheric conditions (T, rainrate, etc.)
    emissions -- dictionary of emission values for each emitted species
    initial_concentrations -- initial values for existant species at start

    Public methods:
    set_condition
    get_condition
    set_initial_concentration
    get_initial_concentrations
    set_emissions
    get_emissions
    output
    """

    def __init__(self, ):
        self.conditions = {}
        self.met = {}
        self.env = {}
        self.emissions = {} # by species
        self.initial_concentrations = {} # by species

        return

    def set_condition(self, condition, value):
        self.conditions[condition] = value
        return

    def set_initial_concentration(self, spc, value):
        self.initial_concentrations[spc] = value

    def get_initial_concentrations(self, ):
        return self.initial_concentrations

    def get_condition(self, condition):
        return self.conditions[condition]

    def set_emissions(self, spc, value):
        self.emissions[spc] = value
        return

    def get_emissions(self, spc):
        return self.emissions.get(spc, 0.)

    def output(self):
        so = ''
        so = add_output_block(so, 'CONDITIONS', self.conditions)
        so = add_output_block(so, 'EMISSIONS', self.emissions)
        so = add_output_block(so, 'METEOROLOGY', self.met)
        so = add_output_block(so, 'ENVIRONMENT', self.env)
        so = add_output_block(so, 'INITIAL CONCENTRATIONS', 
                              self.initial_concentrations)
        return so


def load_yaml(yaml_file):
    """Load a dictionary from a yaml-formatted text file."""
    with open(yaml_file, 'rb') as f:
        yam = yaml.load(f)
    return yam


def make_meta(config):
    """Using config options, set up a Meta object to feed to a Chembox."""
    meta = Meta()
    # Get required fields from config file
    if 'all' in config['OUTPUT_SPECIES']:
        ALLOUT = True
    else:
        ALLOUT = False
        for spc in config['OUTPUT_SPECIES']:
            meta.add_output_species(spc)
    for spc in config['CONSTANT_SPECIES']:
        if type(spc) is bool:
            raise YamlBooleanError('Constant species: %s'%spc)
        meta.add_constant_species(spc)
    for proc in config['PROCESSES']:
        if config['PROCESSES'][proc]:
            meta.add_process(proc)
    for param in config['PARAMETERS']:
        meta.set_parameter(param, config['PARAMETERS'][param])
    rxns = rxnsmod.create_reaction_dict(config['REACTION_FILE'])
    for rxn in rxns:
        if type(rxns[rxn]['k']) == str:
            try:
                rxns[rxn]['k'] = meta.get_parameter(rxns[rxn]['k'])
            except KeyError:
                raise( KeyError, 'rate const. not defined for %s'%rxns[rxn]['k'])
    meta.set_reactions(rxns)
    spcs = rxnsmod.get_species_list(rxns)
    for spc in spcs:
        meta.add_species(spc)
    if ALLOUT:
        for spc in meta.get_species():
            meta.add_output_species(spc)
    # look for non-required fields in config file
    try:
        meta.set_tout(config['TIME_OUT'])
    except KeyError as e:
        print( "make_meta: No TIME_OUT defined in config")
    try:
        tstart = config['START TIME']
    except:
        tstart = 0.
    meta.set_tstart(tstart)
    try:
        for siz in config['BOX_SIZE']:
            meta.set_size(siz, config['BOX_SIZE'][siz])
    except KeyError as e:
        print ("No box size has been defined")
    return meta


def make_inputs(config):
    """Using input options, make an Inputs object for use in a Chembox."""
    inputs = Inputs()

    # Get required input fields for inputs file
    for cond in config['CONDITIONS']:
        inputs.set_condition(cond, config['CONDITIONS'][cond])
    
    # Look for other input fields in inputs file
    try:
        emissions = load_yaml(config['EMISSIONS FILE'])
    except KeyError as e:
        emissions = {}
    for spc in emissions:
        inputs.set_emissions(spc, emissions[spc])
    try:
        ic = config['INITIAL CONCENTRATIONS']
    except:
        ic = {}
    for spc in ic:
        if type(spc) is bool:
            raise YamlBooleanError('Initial condition: %s'%spc)
        inputs.set_initial_concentration(spc, ic[spc])

    return inputs


