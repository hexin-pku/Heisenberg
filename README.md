# Heisenberg
Simple chemistry quantum computation (SCF), version-0.0
(when you see it's 0.0 version, means it still not completed.)
Author, XShinHe <1500011805@pku.edu.cn>

# usage
* Hsbg -d, default
* Hsbg -f *gjf, read from a file, here the file with the Gauss09's format

# the files structure
* Hsbg_Const.h, basic constant definations.
* Hsbg_Tools.h, basic and useful function (mainly string processing function) without in original c++.
* Hsbg_Global.h, creat the namespace __Hsbg__. And defines __Point__, __GPoint__ classes.
* Hsbg_Geom.h, defines __HGeom__ class.
* Hsbg_Orbital.h, defines __HOrbital__ class and __HOrbital\_cgto__ class.
* Hsbg_Basis.h, defines __BPoint__ and __HBasis__ class.
* Hsbg_Parser.h, defines __HTask__ class.
* Hsbg_InteG.h/Hsbg_InteS.h, for the integral of Gauss function or Slater function.
* Hsbg_Slover.h, defines __HSolver__ calss.

# the process of the procedure
1. __main__ read the file name, pass to Parser(or Tasker).
2. Parser preread the file, determine how many atoms it contains. and pass to Geom, creat object.
3. Parser read file, init task and geom parameters.
4. Parser postread deal with the Geom, Creat Basis object, and init each basis with coefficients
and exponents.
5. Pass Task to Solver, Solver creat the array space.
6. Do integral of S, H, ERI integration.
7. Geuss a density matrix P.
8. calculate Fock matrix.
9. transfer Fock to Fock'
10. calculate Coefficient Matrix C'.
11. transfer C' to C.
12. calculate density matrix P and Energy E.
13. Judge whether converge of E or P, otherwise back to calculation Fock matrix again using new P.

# Acknowledgement 


