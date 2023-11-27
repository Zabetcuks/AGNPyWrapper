# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:10:59 2023

A GNPy wrapper for the network simulation with GNPy in Python

@author: Maksims Zabetcuks
"""

import sys
import json
from numpy import linspace
from pathlib import Path
import gnpy.core.ansi_escapes as ansi_escapes
from gnpy.core.elements import Fiber
from gnpy.core.equipment import trx_mode_params
from gnpy.core.parameters import FiberParams
import gnpy.core.exceptions as exceptions
from gnpy.core.network import build_network
from gnpy.core.utils import db2lin, lin2db, automatic_nch, per_label_average, pretty_summary_print
from gnpy.topology.request import (PathRequest, compute_constrained_path, propagate)
from gnpy.tools.json_io import (Span, Roadm, SI, Amp, _automatic_spacing, load_equipment, load_network)
from gnpy.tools.cli_examples import load_common_data

def load_data(path_to_eqpt, path_to_topology):
    """uploads both the equipment library and the network topology from .json files
    
        Parameters:
            path_to_eqpt:               a string with the path to the equipment library file [.json]
            path_to_topology:           a string with the path to the network topology file [.json]
            
        Outputs:
            equipment:                  a dictionary containing the equipment configuration
            network:                    a dictionary containing the network topology
    """
    (equipment, network) = load_common_data(Path(path_to_eqpt), Path(path_to_topology) , None, None)   
    return equipment, network

def load_eqpt(path_to_equipment):
    """uploads the equipment library from .json file
    
        Parameters:
            path_to_eqpt:               a string with the path to the equipment library file [.json]
            
        Output:
            equipment:                  a dictionary containing the equipment configuration
    """
    equipment = load_equipment(Path(path_to_equipment))
    return equipment

def load_net(path_to_network, equipment):
    """uploads the network topology from .json file
        
        Parameters:
            path_to_network:            a string with the path to the network topology file [.json]
            equipment:                  a Python dictionary containing the equipment library (can be created with the load_eqpt function)
            
        Output:
            network:                    a dictionary containing the network topology
    """
    network = load_network(Path(path_to_network), equipment)
    return network
    

def create_path_request(equipment, path_req_params):
    """creates a path request for GNPy
    
        Parameters:
            equipment:                      a dictionary containing the equipment configuration (uploaded with the load_data/load_eqpt functions)
            path_req_params:                a dictionary containing path request and transceiver parameters
                                                
                                                        Full dictionary example:
                                                            
                                                            path_req_params = {
                                                                'request_id': 0,
                                                                'trx_type': '',                             # 'custom' to use custom parameters/'trx_type' from the equipment(requires also trx_mode specification)
                                                                'trx_mode': '',                             # needed and matters if 'trx_type' is not 'custom'
                                                                'source': 'trx_source',                     # uid of the source node in the network topology
                                                                'destination': 'trx_destination',           # uid of the destination node in the network topology
                                                                'bidir': False,                             # bidirectional
                                                                'nodes_list': ['trx_destination'],          # nodes list must contain at least the destination uid from the network topology (see compute_constrained_path)
                                                                'loose_list': ['strict'],                   # ‘strict’ value means that the list of nodes should be strictly followed, while any other value means that the constraint may be relaxed if the node is not reachable
                                                                'format': '',                               # name of CUSTOM constellation, no impact if 'trx_type' & 'trx_mode specified'
                                                                'path_bandwidth': 0,
                                                                'effective_freq_slot': None,
                                                                'baud_rate': 34.17e9, 
                                                                'OSNR': 10.5, 
                                                                'bit_rate': 100e9, 
                                                                'roll_off': 0.15, 
                                                                'tx_osnr': 36, 
                                                                'min_spacing': 50e9, 
                                                                'cost': 1, 
                                                                'penalties': {}, 
                                                                'f_min': 191.3e12, 
                                                                'f_max': 195.1e12, 
                                                                'power': 0.001, 
                                                                'spacing': None                             # set to None if you want automatic spacing or 50GHz = 50000000000.0 = 50e9
                                                                                                            # 'initial_spectrum' can also be assigned with SI from gnpy.tools.json_io (see also example-data initial_spectrum1/2 in gnpy)
                                                                }
                                                            
                                                        Minimal dictionary example:
                                                                
                                                            path_req_params = {
                                                                'request_id': 0,
                                                                'trx_type': '',                             
                                                                'trx_mode': '',                             
                                                                'source': 'trx_source',                     
                                                                'destination': 'trx_destination',           
                                                                'bidir': False,
                                                                'nodes_list': ['trx_destination'],          
                                                                'loose_list': ['strict'],
                                                                'format': '',                               
                                                                'path_bandwidth': 0,
                                                                'effective_freq_slot': None,
                                                                }
                                                            
                                                        If there is no 'trx_type' or 'trx_mode' specified (or they don't match those in the equipment configuration),
                                                        GNPy will use SPECTRAL INFORMATION (default from the equipment configuration/specified, see si_params) as a
                                                        basis for the transceiver parameters.
        Output:
            req:                        a dictionary containing the simulation options
            
    """
    #HERE: if structure for custom/existing transponder
    if path_req_params['trx_type'] != 'custom':
            trx_params = trx_mode_params(equipment, path_req_params['trx_type'], path_req_params['trx_mode'])
            path_req_params.update(trx_params)      
     
    #HERE: if structure for auto/custom spacing
    if path_req_params['spacing'] == None:                                                  #careful: min_spacing must be > baud rate
         path_req_params['spacing'] = _automatic_spacing(path_req_params['baud_rate'])      #automatic spacing from baud rate
        
    #HERE: number of channels computation
    nb_channels = automatic_nch(path_req_params['f_min'], path_req_params['f_max'], path_req_params['spacing'])
    path_req_params['nb_channel'] = nb_channels 
    
    #HERE: creating path request
    req = PathRequest(**path_req_params)
    
    return req

def customize_si(equipment, si_params):
    """customizes the spectral information in the equipment library
    
    Parameters:
        equipment:                          a dictionary containing the equipment configuration (uploaded with the load_data/load_eqpt functions)
        si_params:                          a dictionary containing spectral information parameters. Will be used if no/wrong 'trx_type'/'trx_mode' in path_req_params
                                            
                                            Dictionary example:
                                            
                                                si_parameters = {                           
                                                    "f_min": 191.3e12,
                                                    "f_max": 196.1e12,
                                                    "baud_rate": 32e9,
                                                    "spacing": 50e9,
                                                    "power_dbm": 0,
                                                    "power_range_db": [0, 0, 0.5],
                                                    "roll_off": 0.15,
                                                    "tx_osnr": 40,
                                                    "sys_margins": 0                        # in dB. Added margin on min required transceiver OSNR.
                                                   }
                                                
                                            Also possible implementation:
                                                
                                                si_parameters = {
                                                  "f_min": 191.4e12,
                                                  "f_max":193.1e12,
                                                  "baud_rate": 32e9,
                                                  "slot_width": 50e9,
                                                  "roll_off": 0.15,
                                                  "tx_osnr": 40
                                                },
                                                {
                                                  "f_min": 193.1625e12,
                                                  "f_max":195e12,
                                                  "baud_rate": 64e9,
                                                  "delta_pdb": 3,                           # optional
                                                  "slot_width": 75e9,
                                                  "roll_off": 0.15,
                                                  "tx_osnr": 40
                                                }
    """
    equipment['SI']['default'] = SI(**si_params)
    
def customize_amp(equipment, amp_params):
    """customizes the amplifier to be used in the simulation (amplifiers that will be set automatically after each fiber span)
            
            Parameters:
                equipment:                          a dictionary containing the equipment configuration (uploaded with the load_data/load_eqpt functions)
                amp_params:                         a dictionary containing amplifier parameters
                                                    
                                                    Possible noise figure models:
                                                        
                                                        from collections import namedtuple
                                                        
                                                        Model_vg = namedtuple('Model_vg', 'nf1 nf2 delta_p orig_nf_min orig_nf_max')
                                                        Model_fg = namedtuple('Model_fg', 'nf0')
                                                        Model_openroadm_ila = namedtuple('Model_openroadm_ila', 'nf_coef')
                                                        Model_hybrid = namedtuple('Model_hybrid', 'nf_ram gain_ram edfa_variety')
                                                        Model_dual_stage = namedtuple('Model_dual_stage', 'preamp_variety booster_variety')
                                                    
                                                    Dictionary example ('openroadm_ila_low_noise' from OpenROADM Ver.5 equipment configuration: eqpt_config_openroadm_ver5.json in GNPy example data):
                                                    
                                                        amplifier_parameters = {
                                                            'f_min': 191.35e12,
                                                            'f_max': 196.1e12,
                                                            'type_variety': '',                     # name of the amplifier
                                                            'type_def': 'openroadm',                # possible values: fixed_gain, advanced_model, variable_gain, openroadm, openroadm_preamp, openroadm_booster, dual_stage
                                                            'gain_flatmax': 27,
                                                            'gain_min': 0,
                                                            'p_max': 22,                            # limit for resulting total power during propagation
                                                            'nf_model': Model_openroadm_ila([-8.104e-4, -6.221e-2, -5.889e-1, 37.62]),    #must match with 'type_def'
                                                            'dual_stage_model': None,               # dual-stage amplifier combines two distinct amplifiers. Vendors which provide an accurate description of their preamp and booster stages separately can use the dual-stage model for an aggregate description of the whole amplifier.
                                                            'nf_fit_coeff': None,                   # only for polynomial fit based nf calculation
                                                            'nf_ripple': None,                      # is uploaded from 'default_edfa_config.json'
                                                            'dgt': [0],                             # is uploaded from 'default_edfa_config.json', dynamic gain tilt,
                                                            'gain_ripple': None,                    # is uploaded from 'default_edfa_config.json'
                                                            'out_voa_auto': False,                  # auto-design feature to optimize the amplifier output VOA (variable optical attenuator). If true, output VOA is present and will be used to push amplifier gain to its maximum, within EOL power margins
                                                            'allowed_for_design': True,             # always set to True (or it will be not used)
                                                            'raman': False,
                                                            'pmd': 0,
                                                            'pdl': 0
                                                        }
    
    """
    #structure to upload dgt, nf_ripple, gain_ripple from 'default_edfa_config.json'
    with open('default_edfa_config.json') as default_edfa_config:
        amp_dgt = json.load(default_edfa_config)
        amp_params['dgt'] = list(amp_dgt['dgt'])
        amp_params['nf_ripple'] = list(amp_dgt['nf_ripple']) 
        amp_params['gain_ripple'] = list(amp_dgt['gain_ripple']) 
    #making sure only custom edfa will be used
    for n in equipment['Edfa'].values():
        n.allowed_for_design = False
    #adding custom edfa to dictionary
    equipment['Edfa'].update({'custom': Amp(**amp_params)}) 
    
    
def customize_one_fiber(network, fiber_uid, fiber_params):
    """customizes fiber parameters for a single fiber based on its uid
            
            Parameters:
                network:                            a dictionary containing the network topology (uploaded with the load_data/load_net functions)
                fiber_uid:                          a string containing an uid of the fiber to be customized from the network topology
                fiber_params:                       a dictionary containing fiber parameters
                                                    
                                                    Dictionary example ('SSMF' from standard GNPy equipment configuration):
                                                        
                                                        fiber_parameters = {            
                                                            'length': 450,
                                                            'length_units': 'km',               # possible units: km, m
                                                            'dispersion': 1.67e-05,
                                                            'effective_area': 83e-12,
                                                            'pmd_coef': 1.265e-15,
                                                            'loss_coef': 0.2
                                                            # 'ref_wavelength' OR 'ref_frequency' can also be specified 
                                                        }
                                                        
                                                    'ref_wavelength' = 1550e-9      # conventional central C band wavelength [m]
                                                    'ref_frequency' is calculated with 'ref_frequency' = c/'ref_wavelength'     # with c = speed of light
    """
    for n in network.nodes(): #setting all fibers in network to custom 
        if isinstance(n, Fiber) and n.uid == fiber_uid:            
            n.params = FiberParams(**fiber_params)
            
            
def customize_all_fiber(network, fiber_params):
    """customizes fiber parameters of all fibers in the network (i.e. all fibers will have these values)
    
            Parameters:
                network:                            a dictionary containing the network topology (uploaded with the load_data/load_net functions)
                fiber_params:                       a dictionary containing fiber parameters
                                                    
                                                    Dictionary example ('SSMF' from standard GNPy equipment configuration):
                                                        
                                                        fiber_parameters = {            
                                                            'length': 450,
                                                            'length_units': 'km',               # possible units: km, m
                                                            'dispersion': 1.67e-05,
                                                            'effective_area': 83e-12,
                                                            'pmd_coef': 1.265e-15,
                                                            'loss_coef': 0.2
                                                            # 'ref_wavelength' OR 'ref_frequency' can also be specified 
                                                        }
                                                        
                                                    'ref_wavelength' = 1550e-9      # conventional central C band wavelength [m]
                                                    'ref_frequency' is calculated with 'ref_frequency' = c/'ref_wavelength'     # with c = speed of light
    """
    #HERE: assignment of fiber parameters
    for n in network.nodes(): #setting all fibers in network to custom 
        if isinstance(n, Fiber):            
            n.params = FiberParams(**fiber_params)
            
            
def customize_span(equipment, span_params):
    """customizes span parameters
    
            Parameters:
                equipment:                          a dictionary containing the equipment configuration (uploaded with the load_data/load_eqpt functions)
                span_params:                        a dictionary containing span parameters
                                                    
                                                    Dictionary example: 
                                                                                        
                                                        span_parameters = {
                                                             'power_mode': True,
                                                             'delta_power_range_db': [-2,3,0.5],        #  [min, max, step] power excursion/span
                                                             'max_fiber_lineic_loss_for_raman': 0.25,
                                                             'target_extended_gain': 2.5,
                                                             'max_length': 90,                          # possibe span length from 50 to 90 km
                                                             'length_units': 'km',                      # possible units: km, m
                                                             'max_loss': 28,                            # not used in the current code implementation (from GNPy documentation)
                                                             'padding': 10,                             # in dB. Min span loss before putting an attenuator before fiber
                                                             'EOL': 0,                                  # All fiber span loss ageing. The value is added to the con_out (fiber output connector). So the design and the path feasibility are performed with span_loss + EOL
                                                             'con_in': 0,                               # fiber input connector
                                                             'con_out': 0                               # fiber output connector
                                                         }
    """
    equipment['Span']['default'] = Span(**span_params) 
    
def customize_roadm(equipment, roadm_params):
    """customizes roadm parameters
    
        Parameters:
            equipment:                          a dictionary containing the equipment configuration (uploaded with the load_data/load_eqpt functions)
            roadm_params:                       a dictionary containing roadm parameters
                                                
                                                Dictionary example:
                                                    
                                                    roadm_parameters = {
                                                        'target_pch_out_db': -20,           # other allowed equalisations: 'target_psd_out_mWperGHz', 'target_out_mWperSlotWidth'
                                                        'add_drop_osnr': 33,                # OSNR contribution from the add/drop port
                                                        'pmd': 3e-12,                       # polarization mode dispersion
                                                        'pdl': 1.5,                         # polarization dependent loss
                                                        'restrictions': {
                                                            'preamp_variety_list': [],      # uid of the preferred preamplifier from the equipment library
                                                            'booster_variety_list': []      # uid of the preferred booster from the equipment library
                                                        }
                                                    }

    """
    equipment['Roadm']['default'] = Roadm(**roadm_params)

def simulate(equipment, network, req):
    """builds and simulates the network using the equipment configuration according to the path request
    
            Parameters:
                equipment:                      a dictionary containing the equipment configuration (uploaded with the load_data/load_eqpt functions)
                network:                        a dictionary containing the network configuration (uploaded with the load_data/load_net functions)
                req:                            a path request (created with create_path_request)
                
            Outputs:
                path:                           a Python list containing each crossed network element on the path together with the gains/losses they introduce
                infos:                          an object of the GNPy class gnpy.core.info.SpectralInformation. It contains in a form of Python lists the parameters of each WDM channel:
                                                frequency, baud_rate, slot_width, signal, nli, ase, roll_off, chromatic_dispersion, pmd, pdl, delta_pdb_per_channel, tx_osnr, ref_power, label.
    """
    #HERE: power mode status
    power_mode = equipment['Span']['default'].power_mode
    
    # Keep the reference channel for design: the one from SI, with full load same channels
    pref_ch_db = lin2db(req.power * 1e3)  # reference channel power / span (SL=20dB)
    pref_total_db = pref_ch_db + lin2db(req.nb_channel)  # reference total power / span (SL=20dB)
      
    #HERE: set all the parameters before building the network below
    try:
          build_network(network, equipment, pref_ch_db, pref_total_db, None)
    except exceptions.NetworkTopologyError as e:
          print(f'{ansi_escapes.red}Invalid network definition:{ansi_escapes.reset} {e}')           
          sys.exit(1)
    except exceptions.ConfigurationError as e:
            print(f'{ansi_escapes.red}Configuration error:{ansi_escapes.reset} {e}')
            sys.exit(1)
    
    
    #HERE: actual path computation
    path = compute_constrained_path(network, req)
     
    power_range = [0]
    if power_mode:
            # power cannot be changed in gain mode
     try:
                p_start, p_stop, p_step = equipment['SI']['default'].power_range_db
                p_num = abs(int(round((p_stop - p_start) / p_step))) + 1 if p_step != 0 else 1
                power_range = list(linspace(p_start, p_stop, p_num))
     except TypeError:
                print('invalid power range definition in eqpt_config, should be power_range_db: [lower, upper, step]')
     for dp_db in power_range:
            req.power = db2lin(pref_ch_db + dp_db) * 1e-3
            # if initial spectrum did not contain any power, now we need to use this one.
            # note the initial power defines a differential wrt req.power so that if req.power is set to 2mW (3dBm)
            # and initial spectrum was set to 0, this sets a initial per channel delta power to -3dB, so that
            # whatever the equalization, -3 dB is applied on all channels (ie initial power in initial spectrum pre-empts
            # "--power" option)
            
            infos = propagate(path, req, equipment) 
            """
            print('Infos: ')
            print(infos)
            """
    else: 
        infos = propagate(path, req, equipment)
            
    return path, infos

def final_gsnr(path):
    """computes per label average GSNR over propogated labels at the destination
        
        Parameters:
            path:           a GNPy path created with the 'simulate' function
            
        Output:
            final_gsnr:     average GSNR of the entire WDM spectrum at the destination
    """
    x = per_label_average(path[-1].snr_01nm, path[-1].propagated_labels)
    final_gsnr = float(pretty_summary_print(x))
    return final_gsnr

def print_path(path):
    """prints the path data (from GNPy) in the console
    
        Parameters:
            path            a GNPy path created with the 'simulate' function    
    
    """
    for elem in path:
        print(elem)
        
def sim_reach(path_to_eqpt, path_to_topology, path_req_params, fiber_parameters, step_size):
    """calculates the maximal feasible reach for a certain transceiver (use case example of GNPyWrapper)
    
        Parameters:
            path_to_eqpt:       a string with the path to the equipment library file [.json]
            path_to_topology:   a string with the path to the network topology file [.json]
            path_req_params:    a dictionary containing path request and transceiver parameters (see create_path_request function)
            fiber_parameters:   a dictionary containing fiber parameters (see the customize_one_fiber/customize_all_fiber functions)
            step_size:          an integer defining steps size in which the fiber length will be incremented after each iteration
            
        Outputs:
            fin_gsnr:           final GSNR at maximal optical reach
            fin_fib_length:     maximal optical reach of a transceiver in km
    """
    feasible = True
    equipment = load_eqpt(path_to_eqpt)
    path_request = create_path_request(equipment, path_req_params)
    while feasible:
        network = load_net(path_to_topology, equipment)
        customize_all_fiber(network, fiber_parameters)
        path, infos = simulate(equipment, network, path_request)
        if final_gsnr(path) > path_request.OSNR:
            fiber_parameters['length'] += step_size
        else:
            feasible = False
    fin_gsnr = final_gsnr(path)
    fin_fib_length = fiber_parameters['length']
    return fin_gsnr, fin_fib_length
 
    

