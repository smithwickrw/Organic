# Semiempirical MO calculation for organic molecules
  See also Supplemental Information of “Semiempirical molecular-orbital calculations of the atomization energies of organic molecules”, Smithwick and Roy (2026), Computational Chemistry, Scientific Research Publishing.
  
  To calculate the 117 molecules (plus nicotine and caffeine and CO and C2), one needs only three files:  1) qcp5.f, which has built-in parameters, 2) input2, which has the molecular inputs followed by ‘FINI’ in capital letters to end the program, and 3) FORT8, which is the pre-formed output file that is overwritten.  After the program is compiled and program starts it takes only ~1 second total time to calculate all the molecules.  Hydrogen has one orbital per atom; C, N, O, and F have four orbitals per atom.  After execution of the program, FORT8 should match file FORT8R.
  
  To generate Table S1 of the Supplemental Information, we used QCP6.F with FORT22, INPUT2 and FORT9.  The parameters are read from FORT22.  FORT9 is the pre-formed output file that is overwritten.  After execution of the program, FORT9 should match FORT9R.
  
  To modify the program and to determine new parameters, QCP7.F, FORT22A, FORT22, and FORT8 are enabled.  This is a slower program and one has to set how many times (ILM =1,9 for example) one wants to cycle through the parameters for each of the data sets.   Some sets were omitted in the past; therefore, sets ILL=1, 12 are cycled with sets #5 and #9 skipped to end up with 10 sets.  “ILK=1,12” initializes several parameters prior to the calculation of each molecule.  After the cycles are run, Fort22 should be copied into Fort22A to upgrade the initial parameters to newer parameters.  (Fort22 and Fort22A are different files because an error is generated if the same file is used for both reading and writing.)
  
  The GNU Fortran (GFortran) is an open-source Fortran compiler.  The program that we used is an inexpensive commercial version of GNU Fortran Compiler called Simply FORTRAN.  We found this platform to be convenient for editing and running our legacy code.  The files should all be located in the same computer folder.  Under project options, we used fortran-library.a.exe and Executable mode and 64 bit (x64).
  
  These files may be copied, but please reference our contribution.
