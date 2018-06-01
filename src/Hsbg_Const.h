#ifndef _Hsbg_Const_H_
#define _Hsbg_Const_H_

#define c_pi    3.14159265358979323846264338327950288
#define c_twopi 6.2831853071795862320E0
#define c_h 6.62606896E-34                   /*- The Planck constant (Js) -*/
#define c_c 2.99792458E8                     /*- Speed of light (ms$^{-1}$) -*/
#define c_kb 1.3806504E-23                   /*- The Boltzmann constant (JK$^{-1}$) -*/
#define c_R 8.314472                         /*- Universal gas constant (JK$^{-1}$mol$^{-1}$) -*/
#define c_bohr2angstroms 0.52917720859       /*- Bohr to Angstroms conversion factor -*/
#define c_bohr2m 0.52917720859E-10           /*- Bohr to meters conversion factor -*/
#define c_bohr2cm 0.52917720859E-8           /*- Bohr to centimeters conversion factor -*/
#define c_amu2g 1.660538782E-24              /*- Atomic mass units to grams conversion factor -*/
#define c_amu2kg 1.660538782E-27             /*- Atomic mass units to kg conversion factor -*/
#define c_au2amu 5.485799097E-4              /*- Atomic units (m$@@e$) to atomic mass units conversion factor -*/
#define c_hartree2J 4.359744E-18             /*- Hartree to joule conversion factor -*/
#define c_hartree2aJ 4.359744                /*- Hartree to attojoule (10$^{-18}$J) conversion factor -*/
#define c_cal2J 4.184                        /*- Calorie to joule conversion factor -*/
#define c_dipmom_au2si    8.47835281E-30     /*- Atomic units to SI units (Cm) conversion factor for dipoles -*/
#define c_dipmom_au2debye 2.54174623         /*- Atomic units to Debye conversion factor for dipoles -*/
#define c_dipmom_debye2si 3.335640952E-30    /*- Debye to SI units (Cm) conversion factor for dipoles -*/
#define c_c_au 137.035999679                 /*- Speed of light in atomic units -*/
#define c_hartree2ev 27.21138                /*- Hartree to eV conversion factor -*/
#define c_hartree2wavenumbers 219474.6       /*- Hartree to cm$^{-1}$ conversion factor -*/
#define c_hartree2kcalmol 627.5095           /*- Hartree to kcal mol$^{-1}$ conversion factor -*/
#define c_hartree2kJmol 2625.500             /*- Hartree to kilojoule mol$^{-1}$ conversion factor -*/
#define c_hartree2MHz 6.579684E9             /*- Hartree to MHz conversion factor -*/
#define c_kcalmol2wavenumbers 349.7551       /*- kcal mol$^{-1}$ to cm$^{-1}$ conversion factor -*/
#define c_e0 8.854187817E-12                 /*- Vacuum permittivity (Fm$^{-1}$)-*/
#define c_na 6.02214179E23                   /*- Avagadro's number -*/
#define c_me 9.10938215E-31                  /*- Electron rest mass (in kg)-*/

#define MAX_ELEMENT 118

static const char *elements_label[]=
{"X","H","HE","LI","BE","B","C","N","O","F","NE","NA","MG","AL","SI",
"P","S","CL","AR","K","CA","SC","TI","V","CR","MN","FE","CO","NI",
"CU","ZN","GA","GE","AS","SE","BR","KR","RB","SR","Y","ZR","NB","MO",
"TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I","XE","CS","BA",
"LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
"LU","HF","TA","W","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO",
"AT","RN","FR","RA","AC","TH","PA","U","NP","PU","AM","CM","BK","CF",
"ES","FM","MD","NO","LR","RF","DB","SG","BH","HS","MT","DS","RG",
  "UUB","UUT","UUQ","UUP","UUH","UUS","UUO"};



#endif
