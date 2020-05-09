# turbulence-isotropic-flow
Code to generate velocities and passive scalar turbulent homogeneous flows in a cube w periodic boundary conditions.
Used in IBMs Sp2 with the message passing interface library, coded in Fortran-90
May output realizations of the passive scalar embedded in the flow, velocities.

T
nadmontu
bjobs -q dedicated -u all -lp
jm_status -j
For status:
bhosts

(1) bsub -Is /bin/tcsh -l
module load ...(modules)

module load mpt_default MIPSpro_default


lsload
bjobs

bsub -n N -Is -q dedicatedi tcsh -lbmgroup

bjobs -u all
bjobs -p

(2) f90 -O3 spectral.f -lmpi -lcomplib.sgimath

mpirun -p n03 16, n03 8, ... a.out

bsub -n 8 mpirun a.out
(3) bsub -q dedicated -n 32 -W 720 -o output -m n09 mpijob -sgi mpirun
a.out

or bsub -q dedicated -n 32 -W 720 -o output -m n13 mpijob -sgi  mpirun
a.out

bsub -n 256 -R "span[ptile=128]" -L /bin/tcsh -q dedicated -W 200 -o
output -m "n13 n09" mpijob -sgi mpirun a.out

(if not dedicated: type:
bsub -q shared -n 16 -W 60 -o output -m n02 mpijob -sgi mpirun a.out
                                       (or n03)

bkill jobid
bhosts
(4) bjobs
(5) lsload (look availability of machine)

