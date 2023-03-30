#!/usr/bin/env python coding: utf-8

## This script enables MD simulations of chromatin segments at 200 bp resolution
## This script is adjusted from the tutorial of OpenMiChroM simulation of single chromosome.

## Dependencies of using this script:
## OpenMM
## OpenMiChroM
## (To install OpenMM and OpenMiChroM, please follow their instalation guide)

from OpenMiChroM.ChromDynamics import MiChroM  # OpenMiChroM simulation module
from OpenMiChroM.CndbTools import cndbTools  # analysis tools module

import numpy as np
import argparse


class NuChroM(MiChroM):

    """
    Nucleosome resolution Chromatin Model (NuChroM) inherits the MiChroM Class
    from OpenMiChroM and adjusts the polymer definition of chromatin according
    to the real size and connections of nucleosomes on the chromatin.
    """

    def _initFENEBond(self, kfb=30):
        """
        Internal function that inits FENE bond force. Here, we adjusted the FENE
        bond distances such that the distances between neighboring nucleosomes
        is consistent with the length of linker DNA.
        """

        if "FENEBond" not in list(self.forceDict.keys()):
            force = (
                "- 0.5 * kfb * r0 * r0 * log(1-((r-cutfene)/r0)*((r-cutfene)/r0)) * "
                "step(r - cutfene) + (4 * e * ((s/r)^12 - (s/r)^6) + e) * step(cut - r)"
            )

            bondforceGr = self.mm.CustomBondForce(force)
            bondforceGr.addGlobalParameter("kfb", kfb)
            bondforceGr.addGlobalParameter("r0", 1.5)
            bondforceGr.addGlobalParameter("e", 1.0)
            bondforceGr.addGlobalParameter("s", 1.0)
            bondforceGr.addGlobalParameter("cut", 2.0 ** (1.0 / 6.0))
            bondforceGr.addGlobalParameter("cutfene", 1.5)

            self.forceDict["FENEBond"] = bondforceGr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--iclist",
        type=str,
        required=True,
        help="Path to the file contains ideal chromosome (IC) parameters.",
    )
    parser.add_argument(
        "--type2type",
        type=str,
        required=True,
        help="Path to the file contains parameters for type to type interactions.",
    )
    parser.add_argument(
        "--sysname",
        type=str,
        default="chrom_unknown",
        help="The name of the simulated structure.",
    )
    parser.add_argument(
        "--outpath",
        type=str,
        default="output",
        help="Path to save the output files.",
    )
    parser.add_argument(
        "--outfilename",
        type=str,
        default="chrom_unknown",
        help="Prefix of the output trajectory files.",
    )
    parser.add_argument(
        "--inputseq", type=str, required=True, help="Path to the input compartment sequence file."
    )

    # Parse argument for simulation
    args = parser.parse_args()

    # Set up simulation with system name
    sim = NuChroM(name=args.sysname, temperature=1.0, time_step=0.01)

    # Choose opencl acceleration for the simulation. Other options: 'cpu',
    # 'cuda'
    sim.setup(platform="opencl")

    sim.saveFolder(args.outpath)

    # Initialize the chromatin structure with a spiral spring shape
    mychro = sim.create_springSpiral(ChromSeq=args.inputseq)

    # Load the initial structure into "sim" object
    sim.loadStructure(mychro, center=True)

    # Include the force field for the simulated chromatin **Homopolymer
    # Potentials**
    sim.addFENEBonds(kfb=30.0)
    sim.addAngles(ka=2.0)
    sim.addRepulsiveSoftCore(Ecut=4.0)

    # **Chromosome Potentials**
    #
    # Trained parameters are passed to the simulation here.
    sim.addCustomTypes(mu=1.79, rc=3.43, TypesTable=args.type2type)
    sim.addCustomIC(mu=1.79, rc=3.43, dinit=10, dend=5500, IClist=args.iclist)

    # Note: these valeus for mu and rc were calculated for human GM12878
    # lymphoblastoid cells at nucleosome resolution and can be changed for other
    # species.

    # The last potential to be added is the spherical restraint in order to
    # collapse the initial structure
    conf = (len(sim.type_list_letter) * 200 / (np.pi * 4 / 3 * 0.002)) ** (1 / 3) / 10

    sim.addFlatBottomHarmonic(kr=5 * 10**-3, n_rad=conf)

    # Run a short simulation in order to get a collapsed structure.

    # There are two variables that control the chromosomes simulation steps:
    #
    # **block:** The number of steps performed in each cycle (n_Blocks)<br>
    # **n_blocks:** The number of blocks that will be perfomed. <br>
    block = 5 * 10**2
    n_blocks = 2 * 10**4

    # We can save the radius of gyration of each block to observe the
    # convergence into the collapsed state (the time required here depends on
    # the size of your chromosome)
    rg = []

    for _ in range(n_blocks):
        sim.runSimBlock(block, increment=False)
        rg.append(sim.chromRG())

    # save a collapsed structure in pdb format for inspection
    sim.saveStructure(mode="gro")

    # Perform the production run of the simulation. **run_n_block** defines the
    # number of frames will be saved into the trajectory.
    block = 1000
    run_n_blocks = 125000

    sim.initStorage(filename=args.outfilename)
    for _ in range(run_n_blocks):
        sim.runSimBlock(block, increment=True)  # perform 1 block of the simulation
        sim.saveStructure()

    sim.saveStructure(filename="%s_block_%d" % (args.outfilename, run_n_blocks), mode="gro")

    sim.storage[0].close()


if __name__ == "__main__":
    main()
