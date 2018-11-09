# RIEP
Rotational interaction energy profiling -- useful information

The user needs to write their own batch script (which will depend on their
system) to loop over rotation angles.

Protein and membrane coordinates must be provided in separate .gro files.  Each
file must have the same box size. This box size must be sufficiently large to
encompass a system comprising both the protein and the membrane.

Protein and membrane .gro file names should exclude the '.gro' suffix.

GROMACS location must be defined through environment variable 'EBROOTGROMACS'.
Suffix to executables can be defined through '-suffix' option.
Prefix ('g_' or 'gmx') must be specified (depends on GROMACS version).

Protein-membrane distance is calculated using GROMACS mindist.  The two groups
between which the distance is calculated must be passed by number.  For a
protein/membrane system using the index file created by the script these are 1
('Protein') and 12 ('Other').  These values are hard-coded into the RIEP script.
This may require changing if GROMACS changes.  Extraction of the computed
distance value is  hard-coded to work for the current mindist file format.

Initial placement of protein and membrane puts protein into box centre,
calculated from box dimensions.

The minimum distance is sometimes through the periodic boundary conditions
whereas the translation command (and the analysis of the interaction energies)
assumes it is within the box.  This can result in numerous iterations of the
mindist+translation commands.

Protein-membrane interaction energy is extracted from .edr file from 25-step MD.
Energy groups must be passed by number.  The following values are hard-coded: 47
Coul-SR:Protein-Other 48  LJ-SR:Protein-Other 49  LJ-LR:Protein-Other

Commands for calculating per-residue interaction energies are supplied as
comments.  Future versions will have the choice of rotational profiling vs
per-residue energies for a given rotation as a user argument.
