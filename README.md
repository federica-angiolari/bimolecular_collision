Hello!

I'm here to help you :) Here are some useful information about this code.

INTRODUCTION

This program is created to simulate and analyse bimolecular collision. It is necessary to have an idea in mind for your project: you have to know the distance, the energy collision, if there is allowed vibration and rotation before collision, if you need information about fragmentation or energy transfer.
The script was made by a student, so don't be so cruel about the code's syntax, since it was her first time.

For all information, questions or problems please contact: 

federica.angiolari@gmail.com


PYTHON

Version: 3.6; numpy, os, sys, math, system


PACKAGE'S STRUCTURE

DIRECTORY:

- script: there are all the scripts. They have to be copied in your working dir

  - start.py

  - system.py

  - sm.py

  - coll.py

- input: there are two inputs. It must be chosen considering the external program that you'll use for simulations. It must be copied in working dir


  - start_input_dftbplus

  - start_input_gamess

- examples:Some useful examples to view output and input file.

  - H2+H

  - OH+H

  - Prop+OH

  - ms


EXTERNAL PROGRAM

For the simulation part, the code (structured in 4 scripts, located in /script dir) calls external programs: GAMESS or DFTB+. 
GAMESS: important: you have to specify the scr dir in code's input file. The scr dir must NOT be in the working dir: the script "system.py" will take all the files you need for analysis in scr dir and move it into the working dir when the dynamic end. The basis set is 6-311++G(d,p) and the ab-initio method is B3LYP, if you want to change it, you have to edit the script "start_tot.py", in "def input_gamess". 
DFTB+: important: the skf are 3ob-3-1. It should be parametrized for all the atoms in these sk files.
If you want change parameters such as atom's masses or parameters for sk files, look at the first lines of start.py script.



HOW IT WORKS

The code is really simple. There are two different input:

a) input_dftb+
b) input_gamess

It depends on the external program you will use. When you'll write on it, you must consider ONE SPACE after ":". So for example, in the first row of input_*:

path_target: /dir/target.xyz 

You'll have to write the path of your molecules (in xyz format, already OPTIMIZED), the distance [Angstrom], the limit of impact factor and the impact factor [Angstrom], number angles of target's rotation, vibro-rotation temperature [Kelvin] for both molecule, collision energy[kcal/mol for GAMESS and eV for DFTB+] in CM frame, dt, steps to add (the script calculate min step for collision), path of external program, and other parameters for the external program.

After you have prepared the input you can launch the code (look next section, HOW TO USE & MOD). At this point (if everything is correct), you'll see some printed sentences on screen such as "Number Dynamic: X": this means that the external program is doing the dynamic. Of course, you can check always the job, since the output of the external program is not printed on the screen, it is saved in an external file ("output_program"). After the dynamic, it will be printed on screen "Traject Analyse: X", so this means that the code it's running the analysis. At the end the code will exit automatically. You'll have in your working dir different files of output, that we will see in next section.

IMPORTANT: projectile and target must be choosen wisely! The molecule that MUST be rotated is the TARGET. This means this cannot bemade up of only one atom.
For now, the test was made for these systems (considering AB = Target, and CD or C = Projectile)
-	AB + CD -> A + CDB 
-	AB + CD -> ABC + D
-	AB + C  -> A + CB
-	AB + C  -> ABC


HOW TO USE & MOD

How to launch the code: > python start.py <program> <option>
  
 All possible examples:
a)python start.py DFTB+ coll
b)python start.py DFTB+ ms
c)python start.py GAMESS coll

a) and c): the external program is DFTB+ and GAMESS. 
OUTPUT:
- before dynamic:

  - geometry_i.xyz:  0<i<num_tot_traj. Initial geometry. There will be n geometry for each impact factor, that will have the target in different orientations. This format xyz is faster to check the geometry, since it's easy to read.
  
  - geometry_i.gen:   0<i<num_tot_traj. Format read by DFTB+, initial geometry.
  
  - VELOC.DAT: file for initial velocities for DFTB+, in m/s.
  
  - summ_traj:  there are written all the number of initial geometries, with their impact factor(b) and target's angles orientation.
  
  - dftb_in.hsd: input file for dynamic DFTB+.
  
  - geometry_i.inp: format for input for GAMESS. It has everything!(Velocities [bohr/ps], coordinates, and input for dynamic part).
  
- after dynamic:

  - NVE_i.xyz: trajectories file from DFTB+.
  
  - md_i.out: trajectories's information file from DFTB+.
  
  - geometry_i.trj: trajectories file from GAMESS.
  
  - en_pot_i: potential energy graphic during dynamic. The unit choosen are the default for external program.
  
  - out: file with a column, made by "N","Y". It depends on the type of collision. It's in order of your number traj.
  
  - bi_rotovib_trasf_*, bi_etrasf_trasl_*, bi_etrasf_tot_*: graphic of bi(impact factor) vs energy transfer (%). The energy could be transaltional, roto-vib, or roto-vib + translational.
  
  - map_3d: file with three columns, represented by vx, vy and vz: final velocities of CM of product.
  
  - angle_vs_bi: projectile's deviation angle after collision.


There will be written n_angle geometries for ALL the impact factor chosen. Some points to note: first of all the number of angle chosen for target's rotation are RANDOM. If you need some restriction and want to put angles manually, you have to edit "start.py" in "def angle". The b = 0 will not be considered, since the probability for cross section reaction would be NULL for this impact factor.

b): the external program is DFTB+:
- before dynamic:

  See a) and c).
  The only difference from the other mod ("collision") is made for initial geometries: there are NO different orientation for target, and the impact factor range is from (-max impact factor) (+max impact factor) and the bi = 0 is considered. 

- after dynamic:

  - spectrum_cation: mass spectrum for cation species.
  
  - spectrum_anion: mass spectrum for anion species.
  
  - spectrum_netrual: mass spectrum for neutrasl species.
  
  - spectrum: mass spectrum for all species.
  
  - Fragments: there are written all the numbers of traj with associated fragments and their charges.
  
  - table: table with all fragments generated during all trajectories, with assocciated their charges (as average), occurrence and mass.
