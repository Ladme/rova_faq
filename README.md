# Frequently asked questions

## Gromacs

### How to run a gromacs simulation

1) Prepare a `tpr` file.

`gmx grompp -p {TOP_FILE}.top -c {GRO_FILE}.gro -t {CPT_FILE}.cpt -f {MDP_FILE}.mdp -n {NDX_FILE}.ndx -r {REFERENCE_FILE}.gro -o simulation.tpr`

`-t {CPT_FILE}.cpt` is optional. `-r {REFERENCE_FILE}.gro` is only needed when `itp` restraints are used.

2) Run `mdrun`.

`mpirun -np {NUMBER_OF_MPI_RANKS} gmx_mpi mdrun -v -deffnm simulation -ntomp 1`

### How to center a protein in a gromacs trajectory

1) Using gromacs trjconv:
`yes 1 0 | gmx trjconv -f {XTC_FILE}.xtc -s {TPR_FILE}.tpr -n {NDX_FILE}.ndx -pbc mol -center -o output.xtc`

2) Using [center](https://github.com/Ladme/center):
`center -f {XTC_FILE}.xtc -c {GRO_FILE}.gro -n {NDX_FILE}.ndx -o output.xtc`

### How to fit protein to the reference frame using RMSD

It is better to use a trajectory that is already centered.

`yes 1 0 | gmx trjconv -f {XTC_FILE}.xtc -s {TPR_FILE}.tpr -n {NDX_FILE}.ndx -fit rotxy+transxy -o output.xtc`

or

`yes 1 0 | gmx trjconv -f {XTC_FILE}.xtc -s {TPR_FILE}.tpr -n {NDX_FILE}.ndx -fit rot+trans -o output.xtc`

### How to calculate free energy from umbrella sampling windows

`gmx wham -it tpr_files.dat -ix pullx_files.dat -o G -hist hist -unit kJ`

`tpr_files.dat` contains a list of paths to the `tpr` files from the individual windows.

`pullx_files.dat` contains a list of paths to the `pullx` files from the individual windows.

### How to calculate peptide tilt in a membrane

`gmx gangle -s {TPR_FILE}.tpr -f {XTC_FILE}.xtc -g1 vector -g2 z -group1 'atomnr {FIRST_ATOM} {LAST_ATOM}' -oh angle.xvg`

`{FIRST_ATOM}` and `{LAST_ATOM}` is the first and last CA atom (or backbone bead in Martini) of the peptide.

This assumes that the membrane is built in the xy-plane (with its normal oriented along the z-axis) and that the peptide is a single membrane-spanning alpha-helix.

### How to ionize a membrane system

`gmx genion -s {TPR_FILE}.tpr -p {TOP_FILE}.top -pname {POSITIVE_ION_NAME} -nname {NEGATIVE_ION_NAME} -np {NUMBER_OF_POSITIVE_IONS} -nn {NUMBER_OF_NEGATIVE_IONS} -neutral -o output.gro`

Do **not** use the flag `-conc` because `gmx genion` calculates the number of ions to add to the system based on the volume of the simulation box. However, as the membrane occupies a substantial part of the system, the number of molecules of water is much lower than `gmx genion` assumes and therefore also the number of ions must be lower than calculated by `gmx genion`. When preparing membrane systems, you should always calculate the salt concentration manually (or by using `insane` in case of Martini simulations).

### How to calculate number of ions

_N_ = _c_ (_W_ / _c_<sub>W</sub>)

Where _N_ is the number of molecules of salt, _c_ is the molar concentration of salt, _W_ is the number of molecules of water in the system, _c_<sub>W</sub> is the number of moles of water in one dm<sup>-3</sup> (that is 55.56).

Note that one Martini bead corresponds to 4 water molecules.

### What is physiological salt concentration

0.154 mol dm<sup>-3</sup>

### How to get free energies from AWH simulation

`gmx awh -f {EDR_FILE}.edr -s {TPR_FILE}.tpr -o awh.xvg`

In case your simulation was divided into multiple cycles and/or was run with multiple walkers, `{EDR_FILE}` is `.edr` file from any walker and the last cycle, `{TPR_FILE}` is `.tpr` file from any walker and any cycle.

***

## Martini

### Where to get martinize2

From https://github.com/marrink-lab/vermouth-martinize

### How to martinize an alpha-helical peptide

`martinize2 -f {PDB_FILE}.pdb -o topol.top -x peptide_martinized.pdb -ss {H*PEPTIDE_LENGTH} -p backbone -ff {FORCE_FIELD}`

You can also use `gro` file for the `-f` option.

`{H*PEPTIDE_LENGTH}` is a sequence of `H` characters with a length corresponding to the number of amino acids of the peptide.

`{FORCE_FIELD}` is `martini3001` for Martini 3 or `martini22` for Martini 2.2.

### How to martinize a protein with elastic network 

`martinize2 -f {PDB_FILE}.pdb -x protein_martinized.pdb -o topol.top -elastic -ef {FORCE_CONSTANT} -p backbone -ff {FORCE_FIELD}`

Default `{FORCE_CONSTANT}` is 500.

`{FORCE_FIELD}` is `martini3001` for Martini 3 or `martini22` for Martini 2.2.

### How to set neutral termini for a martini peptide

Add option `-nt` to `martinize2`.

### Where to get insane

From https://github.com/Tsjerk/Insane

### How to prepare a protein-membrane system using insane

Symmetric membrane composed of one lipid type:

`insane -f {PROTEIN}.pdb -m {MARTINI_PHOSPHOLIPIDS}.itp -l {LIPID} -d {DISTANCE_BETWEEN_PERIODIC_IMAGES} -z {Z_DIMENSION} -sol W -o system.gro -p system.top -pbc square -center -salt 0.154`

Symmetric membrane composed of multiple lipid types:

`insane -f {PROTEIN}.pdb -m {MARTINI_PHOSPHOLIPIDS}.itp -l {LIPID1:RELATIVE_ABUNDANCE1} -l {LIPIDS2:RELATIVE_ABUNDANCE2} (...) -d {DISTANCE_BETWEEN_PERIODIC_IMAGES} -z {Z_DIMENSION} -sol W -o system.gro -p system.top -pbc square -center -salt 0.154`

Asymmetric membrane:

`insane -f {PROTEIN}.pdb -m {MARTINI_PHOSPHOLIPIDS}.itp -l {LOWER_LIPID} -u {UPPER_LIPID} -d {DISTANCE_BETWEEN_PERIODIC_IMAGES} -z {Z_DIMENSION} -sol W -o system.gro -p system.top -pbc square -center -salt 0.154`

### How to visualize a Martini peptide using Licorice

(How to add bond information into a pdb file.)

1) `gmx trjconv -f {GRO_FILE}.gro -s {TPR_FILE}.tpr -pbc whole -conect -o output.pdb`

2) In `output.pdb`, remove the line containing `ENDMDL`.

3) Visualize with VMD.

***

## Infinity

### How to activate a module

`module add {MODULE}`

### How to submit a simulation to robox zeros

1) CPU-only simulation

`psubmit cpu {RUN_SCRIPT} ncpus=X,walltime=Y -y`

2) GPU simulation

`psubmit gpu {RUN_SCRIPT} ncpus=X,ngpus=Z,walltime=Y -y`

### How to submit a simulation to sokar cluster

1) CPU-only simulation

`psubmit default {RUN_SCRIPT} props=kraken,ncpus=X,walltime=Y -y`

2) GPU simulation

`psubmit default {RUN_SCRIPT} props=gr_compnode,ncpus=X,ngpus=Z,walltime=Y -y`

### How to submit a simulation to metacentrum

`psubmit default@M {RUN_SCRIPT} ncpus=X,walltime=Y -y`

### How to get a list of your jobs submitted to a cluster

`pjobs`

### How to get a list of all jobs submitted to a cluster

`pqstat`

### How to change a site

`site activate {SITE}`

### How to get a list of available sites

`site`

### How to get the list of nodes of a cluster

`pnodes`

### How to kill a job

`cd {JOB_DIRECTORY} && pkill -y`

or 

`qdel {JOB_NUMBER}`

You can also kill the job softly (the files will be synchronized with the working directory).

`cd {JOB_DIRECTORY} && pkill -s -y`

### How to get to working directory of a job

`pgo`

If the job crashed:

`pgo --force`

### How to get information about a job

`pinfo`


