# Heisenberg
Simple chemistry quantum computation (SCF), version-0.1  
> (when you see it's 0.0 version, means it still not correct. However, its 0.1 version, meaning that it lacks
some functions)  
> Author, XShinHe (He Xin) <1500011805@pku.edu.cn>  
> this work is parallel with the team work "2018QC-Project-Ab-initio-wavefunction-program", a private github
project of Utenaq; this program is first plan to as the programming work for the Quantum Chemistry lecture, 
of Prof. W.J. Liu, CCME of Peking University.
> After the hand in the programmimg work for the QC class, the former basic Hatree Fock method would be put 
in the "tyro" directory with additional annotations, which can be referenced by beginers of Quantum Chemistry.  

# compile
in this directory, by command `make`, or in this child directory type `make`, your computer should be LINUX platform.

# usage
* \[ Hsbg -h \], help
* \[ Hsbg -d \], default
* \[Hsbg -f *.hif/*.gjf\], read from a file, here the file with the Gauss09's format (or with .hif suffixes).

# the files structure
## Hsbg_Const.h, basic constant definations.  
* basic math/physics constant
* unit conversion constant
* elements label
* elements mass
if necessary, you can revise here to define your type of "atoms"
  
## Hsbg_Tools.h, some useful functions
(mainly string processing function) without in original c++.
* \[ getArrayLen \], template getting length of an array.
* \[ trim \], trim a string.
* \[ replace_recursive \], a replace method of string, not recommend for it is recursive
* \[ replace_distinct \], a replace method of string, not recursive. 
  
## Hsbg_Point.h
define the namespace __Hsbg__, where defines __Point__, __GPoint__ basic classes.
* Point class:  
with members __x__, __y__, __z__, __aname__(the atom name);  
with methods __dist__, __dist\_A__, __dist\_A2__, __dist\_AB__, __dist\_AB2__(though most of them are not used);  
with overload operators __\+__, __\-__, __\*__, __>>__.  
* GPoint class: is a child class of Point (inherit from Point).  
with addtional data __znum__, __mass__, __perd__, __fmly__, __indx__, __frag__.  
with addtional methods __get\_Point__(return Point object, polymorphism), __link\_info__(complete the atom information).
  
## Hsbg_Geom.h
* HGeom class: inherit from GPoint. It is a list of GPoint.  
with menbers __Gname__, __Natom__, __Nfrag__, __iatom__, __ifrag__, __geom__(pointer as list of GPoint).  
with method __set\_Geom__, __read\_Geom__.  
with overload operator __<<__.  
  
## Hsbg_Orbital.h
* HOrbital class:  
basis orbital class, can be use as GTO/STO.  
with members __x__, __y__, __z__, __alpha__, __coeff__, __L__, __M__, __N__.(note it not inherit from Point)  
with methhods __init\_Orbit__, __set\_CA__, __get\_CA__, __set\_LMN__, __get\_LMN__, __set\_XYZ__, __get\_XYZ__,
 __get\_Point__(return a Point object), __normGTO__.
* HOrbital_cgto class:  
contract GTO, inherit from HOrbital class.  
with additional members __cn__, __alphas__(list), __coeffs__(list), __orbits__(abandoned!).  
with additional methods __set\_cgto__, __get\_CA__(overload), __set\_Pars__, __get\_BohrL__, __conv\_AUnit__,
 __normGTO__(overload).    
  
## Hsbg_Basis.h
* BPoint class:  
inherit from GPoint  
with additional menbers __ncgto__, __split__, __numcs__, __idxcs__.
with additional methods __setfrom\_GPoint__, __set\_Space__, __set\_BPoint__, __read\_BPoint__, __conv\_AUint__.  
with overload operator __<<__.
* HBasis class:
inherit from BPoint, use as a "list" of BPoint.  
with additional menbers __Bname__, __Natom__, __Allzum__, __basis__(pointer of BPoint), __idmap__(pointer of int).  
with additional methods __set\_Basis__, __set\_Map__, __count\_Allzum__.  
with overload operator __<<__.  
  
## Hsbg_Parser.h
HTask class:  
with members __Hiffile__, __Logfile__, __Maxmem__, __Job__, __Method__, __Basis__, __Title__, __Charge__, __Smulti__,
 __TaskGeom__(HGeom object), __TaskBasis__(Hbasis object), __Natom__, __Nelec__, __N\_set__.  
with methods __set\_IO__, __set\_Job__, __read\_Predo__, __read\_Task__, __read\_Postdo__, __Taskparser__.  
with overload operator __<<__.  
  
## Hsbg_InteG.h/(Hsbg_InteS.h)
for the integral of Gauss function or Slater function (the former only supports ss type, and the latter is 
not completed).  
  
## Hsbg_SCF.h
HScf class:  
the solver of self-consistent field.  
with members (pointer)__tasklink__, __HBlink__, __Scf\_name__, __nsz__, __nocc__, __threshold__, __eigen\_S__, 
__eigen\_H__, __eigen\_ERI__, __eigen\_J__, __eigen\_K__, __eigen\_G__, __eigen\_F__, __eigen\_C__, 
__eigen\_E__, __eigen\_Fp__, __eigen\_Cp__, __eigen\_X__, __eigen\_Y__, __eigen\_P__, __eigen\_P\_old__, 
__E__, __E\_old__, __do\_loop__, __list__.    
with methods __set\_Threshold__, __set\_Space__, __clac\_SHEERI__, __guess\_P__, __calc\_XY__, __calc\_Fock__, 
__calc\_Cprim__, __calc\_PE__, __check\_loop__, __report\_SCF__, __calc_NE__, __m\_Diff__.  
  
  
# the process of the procedure
1. __main__ read the arguments, pass to Parser(or Tasker).  
2. Parser preread the file, determine how many atoms it contains. And pass to Geom, creat object.  
3. Parser read file, initialize task and geom parameters.  
4. Parser postread and deal with the Geom, Creat Basis object, and init each basis with coefficients 
and exponents.  
5. Pass Task to SCF, SCF creat the array space.  
6. Do integral of S, H, ERI integration.
7. Geuss a density matrix P (just use H to from a P).
8. calculate Fock matrix and Fock'.
9. calculate Coefficient Matrix C'.
10. transfer C' to C.
11. calculate density matrix P and Energy E.
12. Judge whether converge of E or P, otherwise back to calculation Fock matrix again using new P.
  
# Acknowledgement 
to Prof. W.J. Liu,  and my teammates in the Quantum Chemistry class.  


