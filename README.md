# Nucleosome resolution Chromatin Model (NuChroM)
Python scripts for conducting nucleosome resolution MD simulations of chromatin structures.

This script uses [OpenMiChroM](https://open-michrom.readthedocs.io/en/latest/) and [openMM](https://openmm.org/) packages. Please refer to these packages for installation guide and dependency requirements.

This script is adapted from the tutorial of OpenMiChroM, [Single Chromosome Simulation](https://open-michrom.readthedocs.io/en/latest/Tutorials/Tutorial_Single_Chromosome.html).

## Usage

```
python3 NuChroM.py [--iclist] [--type2type] [--inputseq] [--sysname] [--outpath] [--outfilename]

```
### Required arguments
* `--iclist`: Parameters for setting ideal chromosome(IC) potentials. Trained IC parameters for NuChroM: `Trained_parameters/IC_gammas`.
* `--type2type`: Parameters for setting type-to-type interactions. For example, when only A and B compartments are used in the simulation, interactions between A-A, A-B and B-B are defined in this file. Trained type-to-type parameters for NuChroM: `Trained_parameters/Type_to_type_alphas`.
* `--inputseq`: The file that contains the sequence of types for the simulated chromatin segment. An example file for chromatin 7 39.5-42.5 Mb region can be found: `examples/chr7_39.5-42.5_AB_beads.txt`.

### Optional arguments
* `--sysname`: Name the chromatin segments that is simulated. Default: `chrom_unknown`.
* `--outpath`: Path for saving the output files. Default: `output`.
* `--outfilename`: Prefix of the output files. Default: `chrom_unknown`.


### Example:
```
python3 NuChroM.py --iclist Trained_parameters/IC_gammas --type2type Trained_parameters/Type_to_type_alphas  --sysname test_chr7 --outpath test_out --outfilename test_chr7_file --inputseq ./examples/chr7_39.5-42.5_AB_beads.txt
```
