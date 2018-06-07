# Heisenberg
Simple chemistry quantum computation (SCF), named after Heisenberg, version-0.1  
  
> (when you see it's 0.0 version, means it still not correct. However, now it's 0.1 version, that means it lacks
some functions)  
  
> Author, XShinHe (He Xin) <1500011805@pku.edu.cn>  
  
> this work is parallel with the team work "2018QC-Project-Ab-initio-wavefunction-program", a private github
project of Utenaq; this program is first plan to as the programming work for the Quantum Chemistry lecture, 
of Prof. W.J. Liu, CCME of Peking University.  

> After the handing in the programmimg work for the QC class, the basic Hatree Fock method (closed-RHF) would be put 
in the "tyro" directory with additional annotations, which can be referenced by beginers of Quantum Chemistry.  
  
# compile
in ./ directory, by command `make`, or in this child directory ./src type `make`. The excutable file will be 
generated in the ./src directory. Note your computer should be LINUX platform because the default compiler is g++.  
  
# usage
* \[ Hsbg -h \], help information
* \[ Hsbg -d \], to calculate HeH+ example by default.
* \[Hsbg -f \*.hif/\*.gjf\], read from a file with the Gauss09's format (or with .hif[^format] format, it's similiar 
with the .gjf format).  
[^format]: this file-format is called Heisenberg Input file .hif, which is similiar to Gauss's <>.gjf file, but 
'%' just stand for an annotation. Similiar, its output file is called Heisenberg output file type (.hof).  
  
# tests and problems
the closed shell of s type Hatree-Fock calculation is fine. But it can not treat with open shell RHF equation 
now or UHF equation. Later will add them.  
|test	|result						| remarks											|  
|------ |-----: 					|:------:											|  
| H2    | consistent with Gauss		| closed-shell RHF works well						|  
| HeH+	| consistent with Gauss		| closed-shell RHF works well						|  
| H		| incorrect					| open-shell RHF / UHF needed[^a] 					|  
| CH4	| (no result)				| integration fucntion of s-p/p-p type is needed	|  
| H4	| incorrect					| maybe convergence/degenerate problem[^b] 			|  
  
[^a]: the test of H, is not correct.(though we not we just diagonize the H matrix to obtain the correct result, but
 the RHF equation of closed shell is not fit to treat it. It's need open-shell RHF SCF method.)  
[^b]: the test of H4, has convergence problem with oscillation ! It seems that it quite need DIIS to help it to
 obtain a correct result.(or to say, if there are degenerate orbitals in the molecule? So the occupation of orbitals should be treat more carefully! P.S., according to Jahn-Teller theorem, non-linear molecule of some symmetry must be degenerate.)  


  
# files structure
## new structure (now updating)
new structure with few files, mainly contains six files, and move them a copy to tyro directory now. From them,
 it's easy to understand how does this procedure work!  
#### Hsbg_Const.h
just as the old file structure.  
* basic math and physics constants, mainly a copy from Pysi4 with little modification.  
* \[static const char*\] elements_label  
* \[static const char*\] elements_mass  

#### Hsbg_Tools.h
just as the old file structure, mainly process the strings.  
* __tempalate \<class T\> int getArrayLen(T& array)__  
* __string& trim(string &s)__ 
* __string& replace_recursive(string& str, const string& old_value, const string& new_value)__  
* __string& replace_distinct(string& str, const string& old_value, const string& new_value)__  

#### Hsbg_Global.h
all modeling class of Heisenberg are put into this file:  
here defines:
* class Point  
> _with members_  
__string name__  
__double x__  
__double y__  
__double z__  
> _with methods_  
__double norm()__  
__double norm2()__  
__Point ref_Point(Point& A)__  
__double ref_norm(Point A)__  
__double re_norm2(Point A)__  
__friend double dist_AB(Point A, Point B)__  
__friend double dist_AB2(Point A, Point B)__
> _with operators_  
__+__  
__-__  
__*__  
__*__  
__<<__  
  
* class Orbital: inherit Point  
> _extra members_  
__int L__  
__int M__  
__int N__  
__double alpha__  
> _extra methods_  
__int set_Alpha(double aa)__  
__int set_LMN(int ll,int mm, int nn)__  
__int set_XYZ(double xx, double yy, double zz)__  
__int set_XYZ(Point P)__  
__Point& get_Point()__  
__int get_LMN(int& ll,int& mm, int& nn)__  
__int get_XYZ(double& xx, double& yy, double& zz)__  
__int conv_AUnit()__  
__double normGTO()__  
__double normGTO( double a)__  
> _with operators_  

* class Orbital_cgto: inherit from Orbital  
> _extra members_  
__int cn__  
the number of (Gauss) Orbital in this contracted GTO.  
__int* coeffs__  
list of coefficients.  
__int* alphas__  
list of exponents.(to be abandoned)  
__Orbital* gtos__  
list of Orbital objects.  
> _extra/overload methods_  
__int set_Cgto(int num)__  
__int get_CA(int idx, double& cc, double& aa)__  
__int set_LMN(int ll,int mm, int nn)__  
__int set_XYZ(double xx, double yy, double zz)__  
__int set_XYZ(Point P)__  
__int conv_AUnit()__  
__int get_BohrL()__  
__double normGTO(int k)__  
> _with operators_  

* class Atom: inherit from Point  
> _extra members_  
__int znum__   
z-number of the atom.  
__double mass__  
the mass of atom (1mol C12 = 12g).  
__int perd__  
__int fmly__  
__int indx__  
__int frag__  
__int ncgto__  
__Orbital_cgto* cgto__  
> _with methods_  
__Point& get_Point()__
__int get_Point(Point& P)__  
__int read_Atom(string line)__  
__int link_Info(string my_name)__  
__int set_Basis(string my_bname, int mysplit, int mynumcs[], int myidxcs[])__  
__int set_Basisspace(int split, int numcs[], int idxcs[])__  
__int read_Basis( fstream& input)__  
__int conv_AUnit()__  
> _with operators_  
__<<__  

* class Molecule: inherit from Atom  
> _extra members_  
__int Natom__  
__int iatom__  
__Atom* atoms__  
> _extra methods_  
(none)  
> _with operators_  
(none)  
* class System: inherit from Molecule  
> _extra members_  
__int Nmol__  
__int imol__  
__Molecule moles__  
__string Bname__  
__int Nbasis__  
__int split__  
__int* numcs__  
__int* idxcs__  
__int* idmap__  
> _extra/overload methods_  
__int set_Natom(int num)__  
__int read_Atom(string line)__  
__int set_Basis(string my_bname)__  
__int set_Map()__  
__Orbital_cgto& OrbC(int &idx)__  
__int count_Znum()__  
__int conv_AUnit()__  
__int solve_Top(string line)__  
> _with operators_  
__<<__
__\[\]__
  
	    
		#################################################################################  
		#                                                                               #  
		#       [ note: ----- inherit, ===== list, =-=-= inherit & list ]               #  
		#                                                                               #  
		#       Point ----- Orbital =-=-= Orbital_cgto                                  #  
		#		   '                     //                                             #  
		#               '               //                                              #  
		#                    '         //                                               #  
		#                         ' Atom  =-=-=-=-= Molecule =-=-=-= System             #  
		#                                                                               #  
		#################################################################################  
	    
  
#### Hsbg_Parser.h
a parser of reading input file.  
* HTask class:  
> _with members_  
__string Hiffile__  
location and name of input file.  
__string Logfile__  
location and name of output file.  
__int Maxmem__  
max memory to be used (not used for now).  
__string Job__  
type of job (now only support sp; scan, opt, freq, qm/mm are not supported now).  
__string Method__  
type of method (now only support closed-shell RHF SCF; 
open-shell RHF, UHF, DFT, MP2, CISD, CCSD are not supported now).  
__string Basis__  
type of basis, such as 6-31g (should be lowercase! later will add a convector).  
__string Title__ 
title, optional.    
__int Charge__  
charge of system.  
__int Smulti__  
spin-multiplicity of system.  
__System Sys__  
object of system.  
__int Natom__  
total atoms of system.  
__int Nelec__
total electrons of system.  
__int Nbasis__
total (contracted) orbital of system.  
  
> _with methods_  
__int set_IO(string Hiffile, string logfile)__  
setup input file and ouput file.  
__int set_Job(string job, string method, string basis)__  
setup job, methed and basis for the task.  
__int read_Predo()__  
pre-read the input file and determine the atom number of the system.  
__int read_Task()__  
read each line of input file, and do analyzation.  
__int Taskparser(string term)__  
a function analyize a line-string, is the job parser!  

> _with overload operator_  
__<<__. 

#### Hsbg_InteG.h  
the integration function. new type of s-p/p-p integration method add in the annotation block (the V & ERI of
them are still not compleled. next only list the completed functions:)  
> _with function_  
__double integral_S(Orbital_cgto& cgto1,  Orbital_cgto& cgto2)__  
__double integral_ERI(  Orbital_cgto& cgto1, Orbital_cgto& cgto2, Orbital_cgto& cgto3, Orbital_cgto& cgto4)__  
__double integral_T(  Orbital_cgto& cgto1,  Orbital_cgto& cgto2)__  
__double integral_V( 	Orbital_cgto& cgto1,  Orbital_cgto& cgto2, Atom& P )__  
__double integral_S_sstype(Orbital_cgto& cgto1,  Orbital_cgto& cgto2)__  
__double integral_ERI_sstype(  Orbital_cgto& cgto1, Orbital_cgto& cgto2, Orbital_cgto& cgto3, Orbital_cgto& cgto4)__  
__double integral_T_sstype(  Orbital_cgto& cgto1,  Orbital_cgto& cgto2)__  
__double integral_V_sstype( 	Orbital_cgto& cgto1,  Orbital_cgto& cgto2, Atom& P )__  

#### Hsbg_SCF.h
the self-consistent field iteration procedure.  
* class HScf  
> _with members_
__ofstream report__   
__HTask* tasklink__  
__System* SYSlink__  
__string Scf_name__  
__int nsz__  
__int nocc__  
__double threshold__  
__MatrixXd 	eigen_S__  
__MatrixXd 	eigen_H__  
__MatrixXd	eigen_G__  
__MatrixXd	eigen_J__  
__MatrixXd	eigen_K__  
__Tensor4D	eigen_ERI__  
__MatrixXd	eigen_F__  
__MatrixXd 	eigen_C__  
__MatrixXd 	eigen_P__  
__MatrixXd	eigen_X__  
__MatrixXd	eigen_Y__  
__MatrixXd	eigen_Fp__  
__MatrixXd	eigen_Cp__  
__bool		do_loop__  
__VectorXd	eigen_E__  
__double		E__  
__double		E_old__  
__MatrixXd	eigen_P_old__  
__int* 		list__  
> _with methods_  
__int set_Threshold(double myeps)__  
__int set_Space(HTask& HT)__  
allocate the matrix size.  
__int calc_SHERI(System &SYS, MatrixXd &S, MatrixXd &H, Tensor4D &G)__  
calculate the H, S, ERI integrations.  
__int guess_P()__  
guess P (density matrix) by a pre-SCF using H as F (Fock Matrix) directly.  
__int calc_XY()__  
calculate the transform matrix, X=S^(-1/2), Y=S^(1/2).  
__int calc_Fock()__  
calculate Fock matrix from H, ERI and P (density) matrix. By the way, calculate F'=X`*F*X .  
__int calc_Cprim()__  
calculate the solution of F'C'= e C' .   
__int tr_Cprim2C()__  
from C' calculate C, by C=X*C' .  
__int calc_PE()__  
from C, calculate P, then from P and eigenvalues of F' calculate E (and add nuclues exclusion energy!).  
__int check_Loop(int cnt)__  
check if it is consistent with itself. (SCF)  
__int report_SCF()__  
print the SCF information.  
__double get_NE()__  
calculate nuclues exclusion energy.  
__double m_Diff(MatrixXd &M, MatrixXd &N)__  
give a kind of norm of matrix. (the elements with max absolution).  
__int loop_SCF()__  
control total SCF cycles.  
 
### Old files structure ( to be adandoned )!
###### Hsbg_Const.h, basic constant definations.  
* basic math/physics constant
* unit conversion constant
* elements label
* elements mass
if necessary, you can revise here to define your type of "atoms"
  
###### Hsbg_Tools.h, some useful functions
(mainly string processing function) without in original c++.
* \[ getArrayLen \], template getting length of an array.
* \[ trim \], trim a string.
* \[ replace_recursive \], a replace method of string, not recommend for it is recursive
* \[ replace_distinct \], a replace method of string, not recursive. 
  
###### Hsbg_Point.h
define the namespace __Hsbg__, where defines __Point__, __GPoint__ basic classes.
* Point class:  
with members __x__, __y__, __z__, __aname__(the atom name);  
with methods __dist__, __dist\_A__, __dist\_A2__, __dist\_AB__, __dist\_AB2__(though most of them are not used);  
with overload operators __\+__, __\-__, __\*__, __>>__.  
* GPoint class: is a child class of Point (inherit from Point).  
with addtional data __znum__, __mass__, __perd__, __fmly__, __indx__, __frag__.  
with addtional methods __get\_Point__(return Point object, polymorphism), __link\_info__(complete the atom information).
  
###### Hsbg_Geom.h
* HGeom class: inherit from GPoint. It is a list of GPoint.  
with menbers __Gname__, __Natom__, __Nfrag__, __iatom__, __ifrag__, __geom__(pointer as list of GPoint).  
with method __set\_Geom__, __read\_Geom__.  
with overload operator __<<__.  
  
###### Hsbg_Orbital.h
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
  
###### Hsbg_Basis.h
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
  
###### Hsbg_Parser.h
HTask class:  
with members __Hiffile__, __Logfile__, __Maxmem__, __Job__, __Method__, __Basis__, __Title__, __Charge__, __Smulti__,
 __TaskGeom__(HGeom object), __TaskBasis__(Hbasis object), __Natom__, __Nelec__, __N\_set__.  
with methods __set\_IO__, __set\_Job__, __read\_Predo__, __read\_Task__, __read\_Postdo__, __Taskparser__.  
with overload operator __<<__.  
  
###### Hsbg_InteG.h/(Hsbg_InteS.h)
for the integral of Gauss function or Slater function (the former only supports ss type, and the latter is 
not completed).  
  
###### Hsbg_SCF.h
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
2. Parser preread the file, determine how many atoms it contains. And pass to creat object.  
3. Parser read file, initialize task and system parameters.  
4. Parser postread and deal with the modeling, Creat System/Atom object, and init basis with coefficients 
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


