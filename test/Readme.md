# the file format
* .gjf, this is the Gauss's input format.
* .hif, this is Heisenberg's input format, note that Heisenberg can use Gauss's input format, but latter may cannot
use the former's.
* .log/.out, this is the Gauss's output format.
* .hof, this is the Heisenberg's output famat.

# the job type
now, the program only support the single point energy calculation job, other job may be ignored and just do sp
calculation. the following job will be added in future:
* scan
* opt
* freq
* sp+freq
* QM/MM

# the basis type
now only support GTO, such as 3-21g, 6-31g, 6-311g etc. But cannot support adding polarization/diffusion.
the sto-3g, here should revise the notation "3-3g" for recent test, or the program would not recognize it.
the following basis would add in the future:
* polarization/diffussion (this can be easy!)
* cc-PVXZ
* STO

# the Integration
the Integration now can only treat s-s type integrals. 
* the s-p and p-p type integration function should add later. (the total frame is completed)
