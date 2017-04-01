# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 21:12:42 2013

@author: thackray

Module containing the reaction parsing tools
"""

import numpy as np
import re
import os

def read_reaction_withk(textr):
    """ Convert text form of a reaction to dictionary for model use 
    Argument: string form of reaction: e.g.
          "C8F17CH2CH2OH + OH -> C8F17CH2CHO     k = 3.2e-11exp(-1000/T)"
    returns: dictionary of {'R': reactants, 'P': products, 'k': rate coeff}
          k elements are [base k, T dependence, type flag]                
    """
    text = re.split('(?i)k', textr)
    R = [ss.strip(' ') for ss in re.split('<*-*>',text[0])[0].split(' + ')]
    P = [ss.strip(' ') for ss in re.split('<*-*>',text[0])[1].split(' + ')]
    k = re.split('(?i)exp\(',text[1])
    if re.match('.*<-*>.*',textr):
        k = [float(text[1].strip(' = ')),0.,0] # partition k
    elif len(k)==1:
        k = k[0].split('(')
        temp = k[1].split('/')
        if temp[0].startswith('T'):
            k = [float(k[0].strip(' = ')), float(temp[1].replace(')','')), 2]
        else:
            k = [float(k[0].strip(' = ')), float(temp[0]), -2]
    else:
        temp = k[1].split('/')
        if temp[0].startswith('T'):
            k = [float(k[0].strip(' = ')), float(temp[1].replace(')','')), 1]
        else:
            k = [float(k[0].strip(' = ')), float(temp[0]), -1]
    return {'R':R, 'P':P, 'k':k}


def read_reaction(textr):
    """ Convert text form of a reaction to dictionary for model use
    Argument: string form of reaction: e.g.                               
    "C8F17CH2CH2OH + OH -> C8F17CH2CHO"           
    returns: dictionary of {'R': reactants, 'P': products, 'k': rate coef name}
    """
    if textr.startswith('#') or textr.startswith('\n') \
       or textr.startswith('\r'):
        return None
    texta = re.split(',',textr)
    text = texta[0]
    R = [ss.strip(' ') for ss in re.split('<*-*>',text)[0].split(' + ')]
    P = [ss.rstrip('\n').rstrip('\r').strip() for ss in \
         re.split('<*-*>',text)[1].split(' + ')]
    if len(texta) == 1:
        if ('<' in text) and ('>' in text):
            k = 'K_' + '+'.join(R) + '<>' + '+'.join(P)
        else:
            k = 'k_' + '+'.join(R) + '>' + '+'.join(P)
    else:
        ktext = texta[1]
        k = [float(x) for x in ktext.split(':')]
        if len(k) == 1:
            k = k+[0.,-1]
        if len(texta) > 2:
            try:
                utext = texta[2]
                u = float(utext)
            except ValueError:
                u = 0.
        else:
            u = 0.
    return {'R':R, 'P':P, 'k':k, 'u':u}



def create_reaction_dict(filename):
    """
    Creates and returns a dictionary with each element containing the 
    dictionary form of each reaction from the reaction definition text 
    file 'filename'
    """
    reactions = {}
    f = open(filename, 'r')
    textlines = f.readlines()
    i = 0
    for line in textlines:
        rxn = read_reaction(line)
        if rxn is not None:
            reactions[str(i+1)] = rxn
            i += 1
    return reactions
        

def get_species_list(reactions):
    """
    From the reaction dictionary 'reactions', get the list of species that 
    will have to be kept track of.
    """
    spc = []
    for key in reactions:
        for r in reactions[key]['R']:
            if r not in spc:
                spc.append(r)
        for p in reactions[key]['P']:
            if p not in spc:
                spc.append(p)
    return spc
    

def get_rate_coeff(rxnk, T):
    """
    Arguments: rxnk: 'k' section of a reaction's dictionary, T: temperature
    Return the temperature-adjusted rate coeff 
    """
    if rxnk[-1] == 1:
        return rxnk[0]*np.exp(T/rxnk[1])
    elif rxnk[-1] == -1:
        return rxnk[0]*np.exp(rxnk[1]/T)
    elif rxnk[-1] == 2:
        return rxnk[0]*(T/rxnk[1])
    elif rxnk[-1] == -2:
        return rxnk[0]*(rxnk[1]/T)
        

if __name__=='__main__':
    reax = ['C8F17CH2CH2OH + OH -> C8F17CH2CHO     k = 3.2e-11exp(-1000/T)',#-1
    'C8F17CH2C(O)OO + NO2 -> C8F17CH2C(O)OONO2    k = 1.1e-11(298/T)',#-2
    'C8F17CH2C(O)OONO2 -> C8F17CH2C(O)OO + NO2    k = 2.8e16exp(-13580/T)',#-1
    'C8F17CH2C(O)OO + HO2 -> C8F17CHO + CO2    k = 1.72e-13exp(1040/T)',#-1
    'C8F17OH -> C7F15COOH           k = 2.3e-6exp(0/T)',#-1
    'GAS <--> PART               ']
    
    for r in reax:
        print( read_reaction(r))
        
