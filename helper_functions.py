import os
import shutil

import pandas as pd
import matplotlib.pyplot as plt


class Selector:
    def __init__(self, chain_name):
        self.selected_chain = chain_name
    def accept_model(self, model):
        return True
    def accept_chain(self, chain):
        return True
    def accept_residue(self, residue):
        return True       
    def accept_atom(self, atom):
        if atom.get_full_id()[2] == self.selected_chain:
            return True
        else:
            return False


def make_dir(dir_name):
    try:
        os.mkdir(str(dir_name))
        print(f'Folder {dir_name.name} is created')
    except:
        print(f'Folder {dir_name.name} already exists')
        # shutil.rmtree(str(dir_name))
        # os.mkdir(str(dir_name))


def add_free_parameters(dict, params_to_update, additional_dict={}):
    if isinstance(params_to_update, str):
        params_to_update = [params_to_update]

    for step_param in params_to_update:
        dict[step_param]['mdp'] = {**dict[step_param]['mdp'],
                                   **dict['free_energy_parameters']['mdp'],
                                   **additional_dict}
    return dict


def plot_step_energy_min(data_xvg, save_dir):
    plt.figure(figsize=(8, 5))
    df = pd.read_csv(data_xvg, sep=r'\s+', header=None, names=['time', 'energy'])
    df.plot('time')
    plt.title("Energy change")
    plt.xlabel("Energy Minimization Step")
    plt.ylabel("Potential Energy, kJ/mol$^{-1}$")
    plt.savefig(str(save_dir / 'step_energy_min.png'), dpi=200)


def plot_step_energy_nvt(data_xvg, save_dir):
    plt.figure(figsize=(8, 4))
    df = pd.read_csv(data_xvg, sep=r'\s+', header=None, names=['time', 'Temperature'])
    df.plot('time')
    plt.title("Temperature during NVT equilibration")
    plt.xlabel("Time, ps")
    plt.ylabel("Temperature, K")
    plt.savefig(str(save_dir / 'step_energy_nvt.png'), dpi=200)


def plot_step_energy_npt(data_xvg, save_dir):
    plt.figure()
    df = pd.read_csv(data_xvg, sep=r'\s+', header=None, names=['time', 'Pressure', 'Density'])
    axes = df.plot('time', subplots=True, layout=(2, 1), figsize=(8, 5), sharex=True)
    plt.suptitle("Pressure and density during NPT equilibration")
    axes[1, 0].set_xlabel("Time, ps")
    axes[0, 0].set_ylabel("Pressure, bar")
    axes[1, 0].set_ylabel("Density, g/L")
    plt.savefig(str(save_dir / 'step_energy_npt.png'), dpi=200)


def plot_step_energy_md(data_xvg, save_dir, ma=5):
    plt.figure()
    df = pd.read_csv(data_xvg, sep=r'\s+', header=None, names=['time', 'Temperature', 'Pressure', 'Density'])
    # axes = df.plot('time', subplots=True, layout=(3, 1), figsize=(10, 8), sharex=True)
    # plt.suptitle("Temperature, pressure and density during MD run")
    # axes[1, 0].set_xlabel("Time, ps")
    # axes[0, 0].set_ylabel("Temperature, K")
    # axes[1, 0].set_ylabel("Pressure, bar")
    # axes[2, 0].set_ylabel("Density, g/L")
    # plt.savefig(str(save_dir / 'step_energy_md.png'), dpi=200)

    df[f'T_MA({ma})'] = df['Temperature'].rolling(ma).mean()
    df[f'p_MA({ma})'] = df['Pressure'].rolling(ma).mean()
    df[f'd_MA({ma})'] = df['Density'].rolling(ma).mean()

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    df.plot(x='time', y='Temperature', ax=ax1)
    df.plot(x='time', y=f'T_MA({ma})', ax=ax1, color='k')
    df.plot(x='time', y='Pressure', ax=ax2)
    df.plot(x='time', y=f'p_MA({ma})', ax=ax2, color='k')
    df.plot(x='time', y='Density', ax=ax3)
    df.plot(x='time', y=f'd_MA({ma})', ax=ax3, color='k')

    fig.suptitle("Temperature, pressure and density during MD run")
    ax1.set_ylabel("Temperature, K")
    ax2.set_ylabel("Pressure, bar")
    ax3.set_ylabel("Density, g/L")
    ax3.set_xlabel("Time, ps")
    fig.tight_layout()
    fig.savefig(str(save_dir / 'step_energy_md.png'), dpi=200)