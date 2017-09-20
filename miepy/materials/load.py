"""
Various functions to load materials and display metadata from the material database
"""

import numpy as np
import matplotlib.pyplot as plt
import yaml
import os
import itertools
from matplotlib.markers import MarkerStyle
import miepy
from miepy.material_functions import data_material

def get_filepath(name, author):
    """get the absolute filepath to the material data of name and author"""
    root = miepy.__path__[0]
    filepath = f"{root}/materials/database/main/{name}/{author}.yml"
    return filepath

def load_material(name, author):
    """Load a material from the database
              name         material name
              author       author of the experimental data
    """

    # return load_material(miepy.__path__[0] + "/materials/ag.npy")
    filepath = get_filepath(name, author)
    with open(filepath, 'r') as f:
        docs = yaml.load(f)
        str_data = docs['DATA'][0]['data']
        list_data = str_data.splitlines()
        for i,item in enumerate(list_data):
            vals = item.split(' ')
            vals = list(map(float, vals))
            list_data[i] = vals

        data = np.array(list_data)
        wav = data[:,0]*1e-6
        n = data[:,1]
        k = data[:,2]
        eps = (n + 1j*k)**2

    return data_material(wav, eps)

def get_authors(material_name):
    """Get a list of possible authors for a given material"""
    filepath = get_filepath(material_name, "")
    directory = filepath[:filepath.rfind('/')]
    files = os.listdir(directory)
    authors = map(lambda f: os.path.splitext(f)[0], files)
    return authors

def wavelength_filter(df, min_wav, max_wav):
    """return a material dataframe for wavelengths between min and max"""
    mask = (df['wavelength']>min_wav)&(df['wavelength']<max_wav)
    return df[mask]

def plot_material_by_author(material_name, wavelength_min=0, wavelength_max=np.inf):
    """For a given material, plot all permittivity data for all authors between a given wavelength range"""
    fig1,(ax1,ax2) = plt.subplots(ncols=2, figsize=plt.figaspect(1/2))

    markers = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p', 'x'))
    markers = itertools.cycle(MarkerStyle.filled_markers)
    colors = itertools.cycle([f'C{i}' for i in range(9)])
    authors = get_authors(material_name)

    for author in authors:
        mat = load_material(material_name, author)
        data = wavelength_filter(mat.data, wavelength_min, wavelength_max)

        color = next(colors)
        marker = next(markers)
        
        ax1.plot(data['wavelength']*1e9, data['eps'].real, 
            label=author, color=color, linewidth=1, marker=marker)
        ax2.plot(data['wavelength']*1e9, data['eps'].imag, 
            label=author, color=color, linewidth=1, marker=marker)
    
    for ax in (ax1,ax2):
        ax.set(xlabel="wavelength (nm)", ylabel="permitivitty")
        ax.legend()
    
    ax1.set_title("Real part")
    ax2.set_title("Imaginary part")

    plt.show()

def material_info_by_author(material_name, wavelength_min=0, wavelength_max=np.inf):
    """Print info on all authors of material between a given wavelength range"""

    authors = get_authors(material_name)

    print(f"Information of {material_name} by author between {wavelength_min*1e9} nm and {wavelength_max*1e9} nm:")
    for author in authors:
        mat = load_material(material_name, author)
        data = wavelength_filter(mat.data, wavelength_min, wavelength_max)
        if not data.empty:
            print(f"\t{author}:  {len(data)} datapoints")

            filepath = get_filepath(name, author)
            with open(filepath, 'r') as f:
                docs = yaml.load(f)
                reference = docs['REFERENCES']
                try:
                    comments = docs['COMMENTS']
                except KeyError:
                    comments = ""
            print('\t\t comments: ', comments) 
            print('\t\t reference: ', reference, '\n') 


if __name__ == "__main__":
    mat = load_material('Ag', 'Johnson')
    # plot_material_by_author('Ag', 300e-9, 1000e-9)
    material_info_by_author('Ag')