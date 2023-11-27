How to install GNPyWrapper:

1. Install GNPy (https://gnpy.readthedocs.io/en/master/install.html)
2. Copy the GNPyWrapper folder to your local repository
3. Import the functions you need:   from GNPyWrap.GNPyWrapper import ...(functions you need)
    where 'GNPyWrap.GNPyWrapper' is the path to the GNPyWrapper.py file in your local repository
    
The GNPyWrapper package (in a form it can be uploaded to pip: https://pypi.org/project/pip/) can be found in GNPyWrap/package
    
How to use the GNPyWrapper:

1. Upload the equipment library and the network topology
    Functions needed: 
        load_data OR load_eqpt AND load_net
    Example data: 
        Network Topology:   a) oopt-gnpy/gnpy/example-data/edfa_example_network.json
                            b) GNPyWrap/LinearTopology.json
                            c) GNPyWrap/meshTopologyExampleV2.json
        Equipment Configuration:    a) oopt-gnpy/gnpy/example-data/eqpt_config.json
                                    b) GNPyWrap/eqpt_config_IKR_wrap.json
2. Customize the equipment and the network if needed
    Functions:
        customize_si
        customize_amp
        customize_one_fiber
        customize_all_fiber
        customize_roadm 
        customize_span
3. Create a path request
    Functions needed:
        create_path_request
4. Simulate
    Functions needed:
        simulate
5. Use the results for your purpose
    Functions:
        final_gsnr
        print_path
        
        
The whole simulation example can be seen on GNPyWrap/TEST.py.