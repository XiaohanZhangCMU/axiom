
About this script:
 Require disl_nuc_hetero.tcl in the same folder. All configurations are dumped to MD++ dirname folder.
 Remove atoms between two neighboring glide set planes. Specify glide-set plane by selecting two points (x,z) in the reference configuration. One on top surface one in the body. Both on (111) glide-set such that their middle points lie on glide set.
 status = 0: Create free surface, surface reconstructed. applied strain.
 status = 1: Read trenched configuration, create ellipse loop of different sizes and relax with a little strain help

Some notations:
 cfg: atom configuration for all atoms without any trench
 atomIdx: [0,1,2...NP], NP = # of atoms of cfg
 nbrlist: nbrlist for all atoms on the slice (totIdx)
 totIdx is the total atoms of two neighboring glide set
 nucleus: indices for atoms that are removed from the slice
            frank partial loop (nucleus) is a subset of totIdx
 dataset: a dictionary { bitarr string : potential energy }

Function list:
  select_atoms_within_ecllipse(): return nucleus of an eclipse shape
  select_totalatoms_on_slice()  : return all atom id on a slice
  write_rawdata()               : write position to .dat file for a group of atoms
  writecfg_fromdata()           : write cn and cfg file for a group of atoms
  removeNucleusAtoms()          : remove a group of atoms from perfect lattice
  find_bdy_atoms()              : find atoms on the boundary of the nucleus ( to be improved)
  build_nbrlist()               : 2d array, each line contains neighboring atoms in slice of a nucleus atom
  set_bits_for_nucleus()        : get bit representation of a nucleus

# exec(open("./scripts/work/frankMEP/global_search.py").read())
