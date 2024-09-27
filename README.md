# excitation_inhibition_topology
Analyzing the topology of inhibitory to excitatory connectivity
This is the custom analysis code for our preprint:
_Specific inhibition and disinhibition in the higher-order structure of a cortical connectome_.

## Requirements
Analyses are using python (3.8 or above) in jupyter notebooks. Additionally, the following python packages are required. All can be installed from pypi.
- numpy
- scipy
- h5py
- pandas
- Connectome-Utilities

Additionally, for advanced topological analyses, the [Connectome Analysis](https://github.com/danielaegassan/connectome_analysis/) package is required.

## Getting started
This project conducts topological analyses of specified connectomes. The basic mode of operation is as follows:
1. Define source connectomes (data and any number of controls)
2. Define subnetworks of the connectomes that are analyzed separately
3. Run analyses for all subnetworks of all connectomes individually. For a specified subnetwork the analysis will generate a set of plots and extract a number of parameters that will be stored in an hdf5 file.
4. Run an analysis that uses the exported data in the hdf5 files to compare source connectomes and characterize variability over subnetworks.

Steps (1) and (2) are conducted using .json-formatted configuration files. They specifcy connectome, subnetwork, where to store results and a number of additional analysis parameters. We provide the config files used in our manuscript as part of this repository in ./configs.

Step (3) is conducted by running the notebook [Topology of neurons...](./Topology%20of%20neurons%20relative%20to%20simplices.ipynb) and entering the path to a chosen configuration file where prompted.

After running (3) for several subvolumes of several connectomes, step (4) is conducted by running the notebook [Plot summary statistics](./Plot%20summary%20statistics.ipynb).

Some additional analyses are conducted in the notebooks [Compare simplices...](./Compare%20simplices%20and%20edge%20participation.ipynb) and [Complexity of...](./Complexity%20of%20neighborhoods%20in%20different%20simplex%20positions.ipynb).

For the exploration of a rewiring rule in Figure 6 of the manuscript, see [below](#rewiring-the-model-connectome).

## Data sources for connectomes
As mentioned, using the configuration files provided will repeat the analyses we conducted for the manuscript. However, you will have to update the reference to the source connectome with the location of where it can be found on your system.

Currently, data can be loaded from two different formats: [SONATA](https://doi.org/10.1371/journal.pcbi.1007696) and [Conntility](https://github.com/BlueBrain/ConnectomeUtilities).

#### Conntility
This is the recommended way to run the analyses. Because the data in this format is already optimized for the types of analyses we conduct here. It is configured as follows:
```
"loading": {
      "conntility": "/path/to/connectome.h5",
      "args": [
        "name_of_connectome"
      ]
    },
```
"name_of_connectome" is the name of a group inside the referenced hdf5 file that contains the connectome. For repeating the analyses we conducted in the manuscript, keep this entry as in the provided config files. Otherwise, see [conntility documentation](https://github.com/BlueBrain/ConnectomeUtilities/blob/main/README.md).

The connectome hdf5 files themselves are available on Zenodo. In the manuscript we compare the Microns connectome to the connectivity in the model of Isbister et al.
- [Microns connectome on Zenodo](https://zenodo.org/doi/10.5281/zenodo.8364069). 
- [Isbister et al. model connectome on Zenodo](https://doi.org/10.5281/zenodo.10079406).

#### Sonata
Loading from SONATA is configured as follows:
```
"loading": {
        "snap": {
        "circuit": "/path/to/circuit_config.json",
        "connectome": "name_of_connectome_to_load",
        "flatmap": "/path/to/flatmap.nrrd"
        }
    }
```
- The circuit_config.json is a file representing an entire model in the Sonata format (see the [Sonata documentation](https://libsonata.readthedocs.io/en/stable/)).
- The "connectome" is a set of synaptic connections that are conceptually bundled together in Sonata. Currently, this code only allows analysis of one connectome at a time.
- The "flatmap" is a vometric file that assigns each voxel a position in cortical coordinates, i.e. horizontal "x" and "y" and a "depth" coordinate.

__Note__: All required files, including the flatmap, are part of our release of [a morphologically detailed model of brain tissue](https://doi.org/10.5281/zenodo.8026353) on Zenodo. To repeat the analyses of our manuscript, keep the "connectome" entry as it is in the provided config files.

## Rewiring the model connectome
One analysis we conducted in our manuscript was making the connectivity of the model more similar to the Microns connectivity by employing a plasticity-inspired rewiring rule. This is run as follows:
1. Parameterize which connectome to rewire and rewiring-specific parameters.
2. Run rewiring and save the rewired connectivity as a new connectome
3. Analyze the new connectome as described [above](#getting-started).

Step (1) is once again conducted in a configuration file of the same format as before. The one we used in the manuscript can be found under ./configs/rewiring
