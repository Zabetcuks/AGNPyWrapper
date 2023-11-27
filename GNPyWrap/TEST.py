# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:49:07 2023

@author: Maksims Zabetcuks
"""

from GNPyWrap.GNPyWrapper import (load_data, load_net, load_eqpt, create_path_request, customize_si, simulate, final_gsnr, customize_amp, 
                                customize_all_fiber, customize_roadm, customize_span, customize_one_fiber, print_path, sim_reach)
from collections import namedtuple

#Params:
eqpt_path = 'eqpt_config_IKR_wrap.json'
topology_path = 'meshTopologyExampleV2.json'
lin_topology_path = 'LinearTopology.json'
#Use 1:
# equipment, network = load_data(eqpt_path, topology_path)
#Use 2:
equipment = load_eqpt(eqpt_path)
network = load_net(topology_path, equipment)


"""TESTING SI"""
# #Params:
# si_parameters = { 
#         "f_min": 191.3e12,
#         "f_max": 196.1e12,
#         "baud_rate": 32e9,
#         "spacing": 50e9,
#         "power_dbm": 0,
#         "power_range_db": [0, 0, 0.5],
#         "roll_off": 0.15,
#         "tx_osnr": 20,
#         "sys_margins": 0
#         }
# #Use:
# customize_si(equipment, si_parameters) 
"""WORKS"""


"""TESTING FIBER"""
#Params:
fiber_parameters = { #SSMF
    'length': 90,
    'length_units': 'km',
    'dispersion': 1.67e-05,
    'effective_area': 83e-12,
    'pmd_coef': 1.265e-15,
    'loss_coef': 0.2
    #'ref_wavelength' OR 'ref_frequency' can also be specified
}
# fiber_uid = 'fiber (Quimper → Lorient_KMA)-'
# #Use:
# customize_all_fiber(network, fiber_parameters)
# # customize_one_fiber(network, fiber_uid, fiber_parameters)

"""WORKS"""

"""TESTING SPAN"""
# # Params:
# span_parameters = {
#       'power_mode': True,
#       'delta_power_range_db': [-2,3,0.5],
#       'max_fiber_lineic_loss_for_raman': 0.25,
#       'target_extended_gain': 2.5,
#       'max_length': 50,  
#       'length_units': 'km',
#       'max_loss': 28,
#       'padding': 10,
#       'EOL': 0,
#       'con_in': 0,
#       'con_out': 0
#   }
# #Use:
# customize_span(equipment, span_parameters)
"""WORKS"""

"""TESTING ROADM"""
# # Params:
# roadm_parameters = {
#     'target_pch_out_db': -20,            
#     'add_drop_osnr': 33,
#     'pmd': 3e-12,
#     'pdl': 1.5,
#     'restrictions': {
#         'preamp_variety_list': [],
#         'booster_variety_list': []
#     }
# }
# #Use:
# customize_roadm(equipment, roadm_parameters)
"""WORKS"""

"""TESTING AMP"""
# #Params:
# Model_vg = namedtuple('Model_vg', 'nf1 nf2 delta_p orig_nf_min orig_nf_max')
# Model_fg = namedtuple('Model_fg', 'nf0')
# Model_openroadm_ila = namedtuple('Model_openroadm_ila', 'nf_coef')
# Model_hybrid = namedtuple('Model_hybrid', 'nf_ram gain_ram edfa_variety')
# Model_dual_stage = namedtuple('Model_dual_stage', 'preamp_variety booster_variety')

# amplifier_parameters = {
#     'f_min': 191.35e12,
#     'f_max': 196.1e12,
#     'type_variety': '', 
#     'type_def': 'openroadm', 
#     'gain_flatmax': 27,
#     'gain_min': 0,
#     'p_max': 22,
#     'nf_model': Model_openroadm_ila([-0.0008104, -0.06221, -0.5889, 37.62]), 
#     'dual_stage_model': None,
#     'nf_fit_coeff': None,
#     'nf_ripple': None, 
#     'dgt': [0], 
#     'gain_ripple': None, 
#     'out_voa_auto': False,
#     'allowed_for_design': True, 
#     'raman': False,
#     'pmd': 0,
#     'pdl': 0
# }
# #Use:
# customize_amp(equipment, amplifier_parameters)
"""WORKS"""

"""ALL TOGETHER WORK AS WELL"""
         

# Params:
req_params_lin = {
    'request_id': 0,
    'trx_type': 'OpenZR+', 
    'trx_mode': '100ZR+, DP-QPSK',
    'source': 'trx_source', 
    'destination': 'trx_destination', 
    'bidir': False,
    'nodes_list': ['trx_destination'], 
    'loose_list': ['strict'],
    'format': '100 Gbit/s, DP-QPSK',
    'path_bandwidth': 0,
    'effective_freq_slot': None,
    'baud_rate': 34170000000.0, 
    'OSNR': 10.5, 
    'bit_rate': 100000000000.0, 
    'roll_off': 0.15, 
    'tx_osnr': 36, 
    'min_spacing': 50000000000.0, 
    'cost': 1, 
    'penalties': {}, 
    'f_min': 191.3e12, 
    'f_max': 195.1e12, 
    'power': 0.001, 
    'spacing': None
    } 

req_params_mesh = {
    'request_id': 0,
    'trx_type': 'OpenZR+', 
    'trx_mode': '300ZR+, DP-8QAM',
    'source': 'trx Brest_KLA', 
    'destination': 'trx Vannes_KBE', 
    'bidir': False,
    'nodes_list': ['trx Vannes_KBE'], 
    'loose_list': ['strict'],
    'format': '100 Gbit/s, DP-QPSK',
    'path_bandwidth': 0,
    'effective_freq_slot': None,
    'baud_rate': 34170000000.0, 
    'OSNR': 10.5, 
    'bit_rate': 100000000000.0, 
    'roll_off': 0.15, 
    'tx_osnr': 36, 
    'min_spacing': 50000000000.0, 
    'cost': 1, 
    'penalties': {}, 
    'f_min': 191.3e12, 
    'f_max': 195.1e12, 
    'power': 0.001, 
    'spacing': None
    }     
#Use:   
path_request = create_path_request(equipment, req_params_mesh)

#Use:
path, infos = simulate(equipment, network, path_request)

print_path(path)

x = final_gsnr(path)
print('Final GSNR in dB')
print(x)

step_size = 90
fin_gsnr, fib_length = sim_reach(eqpt_path, lin_topology_path, req_params_lin, fiber_parameters, step_size)
print('Final GSNR ')
print(fin_gsnr)
print('Final fiber length in km: ')
print(fib_length)

