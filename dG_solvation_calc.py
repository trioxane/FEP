#!/usr/bin/env python3

import sys
import os
import time
import pprint
import glob
import shutil

from pathlib import Path
import helper_functions as hf
from Bio.PDB import PDBParser, PDBIO

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_md.gromacs.pdb2gmx import pdb2gmx
from biobb_md.gromacs.editconf import editconf
from biobb_md.gromacs.solvate import solvate
from biobb_md.gromacs.grompp import grompp
from biobb_md.gromacs.genion import genion
from biobb_md.gromacs.mdrun import mdrun
from biobb_analysis.gromacs.gmx_energy import gmx_energy
from biobb_pmx.pmx import pmxanalyse


# Receiving the input configuration file (YAML)
conf = settings.ConfReader('dG_solvation_config.yaml')  # sys.argv[1]
global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
parameters = conf.get_prop_dic(global_log=global_log)
global_paths = conf.get_paths_dic()

wd = Path(conf.get_working_dir_path())
hf.make_dir(wd)

reader = PDBParser()
structure = reader.get_structure(id='input_complex',
                                 file=conf.properties['complex_structure'])

for protein_chain in conf.properties['selected_chains']:

    chain_dir = wd / protein_chain
    hf.make_dir(chain_dir)
    mol_path = str(chain_dir / f'protein_{protein_chain}.pdb')

    io = PDBIO()
    io.set_structure(structure)
    io.save(mol_path, hf.Selector(protein_chain))

    chain_parameters = conf.get_prop_dic(protein_chain, global_log=global_log)
    chain_paths = conf.get_paths_dic(protein_chain)
    chain_paths['step_pdb2gmx']['input_pdb_path'] = mol_path

    global_log.info("step_pdb2gmx: Generate the topology")
    pdb2gmx(**chain_paths["step_pdb2gmx"], properties=chain_parameters["step_pdb2gmx"])

    global_log.info("step_editconf: Create the solvent box")
    editconf(**chain_paths["step_editconf"], properties=chain_parameters["step_editconf"])

    global_log.info("step_solvate: Fill the solvent box with water molecules")
    solvate(**chain_paths["step_solvate"], properties=chain_parameters["step_solvate"])

    global_log.info("step_grompp_genion: Preprocess ion generation")
    grompp(**chain_paths["step_grompp_genion"], properties=chain_parameters["step_grompp_genion"])

    global_log.info("step_genion: Ion generation")
    genion(**chain_paths["step_genion"], properties=chain_parameters["step_genion"])

    lambdas = list(zip(conf.properties['free_energy_parameters']['properties']['mdp']['vdw-lambdas'].split(),
                       conf.properties['free_energy_parameters']['properties']['mdp']['coul-lambdas'].split()))

    for i, lbds in enumerate(lambdas):

        window = f'lambdas-vdw_{lbds[0]}-coul_{lbds[1]}'
        current_lambda_dir = chain_dir / f'window_{window}'
        window_parameters = conf.get_prop_dic(os.path.join(protein_chain, f'window_{window}'), global_log=global_log)
        window_paths = conf.get_paths_dic(os.path.join(protein_chain, f'window_{window}'))

        window_paths['step_grompp_min']['input_gro_path'] = chain_paths['step_genion']['output_gro_path']
        window_paths['step_grompp_min']['input_top_zip_path'] = chain_paths['step_genion']['output_top_zip_path']
        window_paths['step_grompp_nvt']['input_top_zip_path'] = chain_paths['step_genion']['output_top_zip_path']
        window_paths['step_grompp_npt']['input_top_zip_path'] = chain_paths['step_genion']['output_top_zip_path']
        window_paths['step_grompp_md']['input_top_zip_path'] = chain_paths['step_genion']['output_top_zip_path']

        window_parameters = hf.add_free_parameters(window_parameters,
                                                   ['step_grompp_min', 'step_grompp_nvt',
                                                    'step_grompp_npt', 'step_grompp_md'],
                                                   additional_dict={'init-lambda-state': i,
                                                                    'couple-moltype': f'Protein_chain_{protein_chain}',
                                                                    'nstenergy': '250'
                                                                    # '': ''
                                                                    })
        # pprint.pprint(window_paths["step_energy_min"], indent=2)

        global_log.info(f"step_grompp_min: Starting minimization for window {i}")
        grompp(**window_paths["step_grompp_min"], properties=window_parameters['step_grompp_min'])

        global_log.info("step_mdrun_min: Execute energy minimization")
        mdrun(**window_paths["step_mdrun_min"], properties=window_parameters['step_mdrun_min'])

        hf.make_dir(current_lambda_dir / 'step_energy_min')
        global_log.info("step_energy_min: Check energy minimization result")
        gmx_energy(**window_paths["step_energy_min"], properties={'terms': ["Potential"]})
        hf.plot_step_energy_min(window_paths['step_energy_min']['output_xvg_path'], save_dir=current_lambda_dir)


        global_log.info(f"step_grompp_nvt: Preprocess system temperature equilibration for window {i}")
        grompp(**window_paths["step_grompp_nvt"], properties=window_parameters['step_grompp_nvt'])

        global_log.info("step_mdrun_nvt: Execute system temperature equilibration")
        mdrun(**window_paths["step_mdrun_nvt"], properties=window_parameters['step_mdrun_nvt'])

        hf.make_dir(current_lambda_dir / 'step_energy_nvt')
        global_log.info("step_energy_nvt: Check temperature during NVT equilibration")
        gmx_energy(**window_paths["step_energy_nvt"], properties={'terms': ["Temperature"]})
        hf.plot_step_energy_nvt(window_paths["step_energy_nvt"]['output_xvg_path'], save_dir=current_lambda_dir)


        global_log.info(f"step_grompp_npt: Preprocess system pressure equilibration for window {i}")
        grompp(**window_paths["step_grompp_npt"], properties=window_parameters['step_grompp_npt'])

        global_log.info("step_mdrun_npt: Execute system pressure equilibration")
        mdrun(**window_paths["step_mdrun_npt"], properties=window_parameters['step_mdrun_npt'])

        hf.make_dir(current_lambda_dir / 'step_energy_npt')
        global_log.info("step_energy_npt: Check Density & Pressure during NPT equilibration")
        gmx_energy(**window_paths["step_energy_npt"], properties={'terms': ["Pressure", "Density"]})
        hf.plot_step_energy_npt(window_paths["step_energy_npt"]['output_xvg_path'], save_dir=current_lambda_dir)


        global_log.info(f"step_grompp_md: Preprocess free dynamics for window {i}")
        grompp(**window_paths["step_grompp_md"], properties=window_parameters['step_grompp_md'])

        global_log.info("step_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**window_paths["step_mdrun_md"], properties=window_parameters['step_mdrun_md'])

        hf.make_dir(current_lambda_dir / 'step_energy_md')
        global_log.info("step_energy_md: Check Density & Pressure during MD run")
        gmx_energy(**window_paths["step_energy_md"], properties={'terms': ["Temperature", "Pressure", "Density"]})
        hf.plot_step_energy_md(window_paths["step_energy_md"]['output_xvg_path'], save_dir=current_lambda_dir)

    # pprint.pprint(window_parameters['step_mdrun_min'], indent=2)
    # print('#' * 100)

    global_log.info("Analyse MD run data")
    for n, xvg_file in enumerate(glob.glob(f'{wd}/{protein_chain}/window_*/step_mdrun_md/dhdl.xvg')):
        shutil.copy2(xvg_file, f'{wd}/{protein_chain}/dhdl_{n}.xvg')
    os.system(f'gmx bar -f {wd}/{protein_chain}/dhdl_*.xvg -o -oi -oh')
