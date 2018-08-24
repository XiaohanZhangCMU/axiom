# -*-shell-script-*-
# compute energy barrier of homogeneous dislocation nucleation in BCC W
# relax all other stress components (except yz) to zero
#
# make fs build=R SYS=mc2
# sw_mc2_mpich scripts/disl_nuc_hetero.tcl 0 0 001 1     -generates 0.01 0.02 0.03 strain perfect states A
# sw_mc2_mpich scripts/disl_nuc_hetero.tcl 1 0.030 001 1 -generates 0.03 state B.
# mpirun -np $ncpu sw_mc2_mpich scripts/disl_nuc_hetero.tcl 4 0.030 001 1 ----string_relax_parallel
source "$::env(MDPLUS_DIR)/scripts/Examples/Tcl/startup.tcl"

#*******************************************
# Definition of procedures
#*******************************************
proc initmd { n } {
#MD++ setnolog
MD++ setoverwrite
#max number neigbors
set myname [ MD++_Get "myname" ]
readmeam-according-to-myname $myname
#max number neigbors
MD++ NNM = 200
}

proc readpot { } { MD++ {
#--------------------------------------------
#Read in potential file
#
potfile = $::env(MDPLUS_DIR)/potentials/w_pot readpot
} }

#------------------------------------------------------------
proc readmeam-lammps { } { 
#Read in MEAM potential (Baskes format)
MD++ meamfile = "~/Planet/Libs/MD++UMB.svn3/potentials/MEAMDATA/meamf"
#MD++ meafile = "~/Planet/LIbs/MD++UMB.svn3/potentials/MEAMDATA/AuSi2nn.meam" 
MD++ nspecies = 2  element0 = "Siz" element1 = "Ge" 
MD++ {
rcut = 4.5  readMEAM
NNM = 300
} }

# make sure the coordinate is right hand sided.

proc make_perfect_crystal { nx ny nz } {
    MD++ crystalstructure = diamond-cubic latticeconst =  5.4309529817532409 #(A) for Si
    MD++ latticesize = \[  1 1 0  $nx  -1 1 0  $ny  0 0 1  $nz \]
    MD++ makecrystal #finalcnfile = perf.cn writecn #eval
}

#0: Si. 1 : Ge.
proc set_all_atoms_species { id } {
    MD++ fixedallatoms
    MD++ input = $id setfixedatomsspecies
    MD++ freeallatoms
}

proc readmeam-according-to-myname { myname  } {
 if { [ string match "*baskes*" $myname ] } {
  puts "readmeam-baskes"
  readmeam-baskes
 } elseif { [ string match "*lammps*" $myname ] } {
  puts "readmeam-lammps"
  readmeam-lammps
 } elseif { [ string match "*meam*" $myname ] } {
  puts "readmeam"
  readmeam
 } elseif { [ string match "*eam*" $myname ] } {
  puts "readeam"
  readeam
 } else {
  puts "not an eam potential, not reading any files"
 }
}


#--------------------------------------------
proc relax_fixbox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 2e-6 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 1
relax
} }
#end of proc relax_fixbox

#--------------------------------------------
proc relax_freebox { } { MD++ {
# Conjugate-Gradient relaxation
conj_ftol = 1e-4 conj_itmax = 1000 conj_fevalmax = 1000
conj_fixbox = 0
conj_fixboxvec = [ 0 1 1
                   1 1 1
                   1 1 0 ]
relax
} }
#end of proc relax_fixbox

proc setup_window { } { MD++ {
#------------------------------------------------------------
#colors for Central symmetry view
color00 = "red" color01 = "blue" color02 = "green"
color03 = "magenta" color04 = "cyan" color05 = "purple"
color06 = "gray80" color07 = "white" color08 = "orange"
#--------------------------------------------
# Plot Configuration
#
atomradius = [1.0 0.78] bondradius = 0.3 bondlength = 2.8285 #for Si
win_width=600 win_height=600
#atomradius = 0.9 bondradius = 0.3 bondlength = 0 #2.8285 #for Si
atomcolor = orange highlightcolor = purple  bondcolor = red
fixatomcolor = red backgroundcolor = gray70
#atomcolor = lightgrey highlightcolor = purple  bondcolor = darkgrey
plot_color_windows = [ 0
                 -5   -4   0
                 -5.5 -5   8   #color00 = red
                 -6  -5.5  1
               -6.55 -6    2
               -6.72 -6.55 3
                 -6.5 -6   6   #color06 = gray80
                 -7.5 -6.5 4
                ]

#xiaohan, if you want to plot the dislocations, uncomment the following

plot_color_axis = 2  NCS = 4
plot_color_windows = [ 0
                       0.6 9.4   1  
                       9.4  10   5
                       10 20   6
                       20 50   8
                       0  0.6  4
                     ]

#plot_limits = [ 1 -10 10 -10 10.0 0.05 10 ]


#
#xiaohan

#plot_color_windows = 5
#plot_color_windows = 0 
plot_atom_info = 1 # reduced coordinates of atoms
plot_limits = [ 1 -10 10  -0.2 0.2 0 10 ]
#plot_atom_info = 2 # real coordinates of atoms
#plot_atom_info = 3 # energy of atoms
#plot_highlight = [ 0 0 1 2 3 4 5 6 7 8 9 ]
plotfreq = 10
#

rotateangles = [ -0 90 0 1.2 ]

#rotateangles = [ 0 0 0 1.7 ]
#rotateangles = [ 0 -90 0 1.7 ]
#openwin alloccolors rotate saverot plot
#plot_color_axis = 0 input = [ -8 -3 10] GnuPlotHistogram
#plot_color_axis = 2 input = [ 0.6 50 50 ] GnuPlotHistogram
} }

proc openwindow { } { 
setup_window
MD++ openwin alloccolors rotate saverot eval plot
}

#--------------------------------------------
proc exitmd { } { MD++ quit }
#end of proc exitmd
#--------------------------------------------

#--------------------------------------------
proc setup_md { } { MD++ {     
T_OBJ = 300 #Kelvin #add by xiaohan

equilsteps = 0  totalsteps = 5000 timestep = 0.0001 # (ps)
atommass = 28.0855 # (g/mol)
DOUBLE_T = 1
saveprop = 1 savepropfreq = 100 openpropfile #run
savecn = 1 savecnfreq = 10000 openintercnfile
plotfreq = 100 printfreq = 100
#ensemble_type = "NPH" integrator_type = "Gear6" implementation_type = 0
ensemble_type = "NVE" integrator_type = "VVerlet" implementation_type = 0
vt2 = 1e28  #1e28 2e28 5e28
wallmass = 2e3     # atommass * NP = 14380
boxdamp = 1e-3     # optimal damping for 216 atoms and wallmass 1e-3
saveH # Use current H as reference (H0), needed for specifying stress
fixboxvec = [ 0 0 1
              1 0 1
              0 0 0 ]
output_fmt = "curstep EPOT KATOM Tinst HELM HELMP TSTRESS_xx TSTRESS_yy TSTRESS_zz H_11 H_22 H_33" 
} }
#end of proc setup_md

proc myRand {min max} { return [expr int(rand()*($max-$min+1)) + $min] }

#*******************************************
# Main program starts here
#*******************************************
# status 0:
#        1:
#        2:
#
# read in status from command line argument
if { $argc <= 0 } {
 set status 0
} elseif { $argc > 0 } {
 set status [lindex $argv 0]
}
puts "status = $status"

if { $argc <= 1 } {
 set n 0
} elseif { $argc > 1 } {
 set n [lindex $argv 1]
}
if { $argc <= 2 } {
 set flag 0
} elseif { $argc > 2 } {
 set flag [lindex $argv 2]
}
if { $argc <= 3 } {
 set opt1 0
} elseif { $argc > 3 } {
 set opt1 [lindex $argv 3]
}
if { $argc <= 4 } {
 set opt2 0
} elseif { $argc > 4 } {
 set opt2 [lindex $argv 4]
}
if { $argc <= 5 } {
 set opt3 0
} elseif { $argc > 5 } {
 set opt3 [lindex $argv 5]
}
if { $argc <= 6 } {
 set opt4 0
} elseif { $argc > 6 } {
 set opt4 [lindex $argv 6]
}
if { $argc <= 7 } {
 set opt5 0
} elseif { $argc > 7 } {
 set opt5 [lindex $argv 7]
}
if { $argc <= 8 } {
 set opt6 0
} elseif { $argc > 8 } {
 set opt6 [lindex $argv 8]
}
if { $argc <= 9 } {
 set opt7 0
} elseif { $argc > 9 } {
 set opt7 [lindex $argv 9]
}

if { $argc <= 10 } {
 set opt8 0
} elseif { $argc > 10 } {
 set opt8 [lindex $argv 10 ]
}

if { $status == 0 } {
  MD++ setnolog
  initmd $n
  if { $argc > 4 } {
    MD++ dirname = $opt2
  } else { 
    puts "Need to specify dirname from python"
    MD++ quit
  }
 
  readpot
  make_perfect_crystal 12 12 16
  MD++ eval saveH conj_fixboxvec =  \[ 0 1 0  0 0 0   1 1 0 \]
  MD++ conj_fixbox = 1 conj_ftol = 2e-2 conj_itmax = 1000 conj_fevalmax = 2000
  MD++ relax finalcnfile = "0K_perfect.cn" writecn
  MD++ eval

  set dighole 0
  set make_free_surface 1  
  set surface_reconstruct 1
  
  if { $dighole } {
     MD++ input = \[ 1 -10 -0.25  -10 10   0.25 10 \] fixatoms_by_position
     MD++ removefixedatoms
     MD++ freeallatoms      
  }

  if {$make_free_surface} {
    MD++ vacuumratio = [expr 1.0-1.0/(1.0 + 0.2)]
    MD++ conj_fixbox = 0 
    if { $flag == 001 } {
       # create 001 surface
       MD++ input = \[ 3 3 0.2 \] changeH_keepR relax finalcnfile = "0K_surf${flag}.cn" writecn
       set strain_data_file "strain_surf_001.dat"
       set stress_data_file "stress_surf_001.dat"
       MD++ { conj_fixboxvec = [ 0 1 1
                                 1 1 1
                                 1 1 1 ] }
    }
    MD++ conj_fixbox = 1
  }

  #  move top and bot layer atoms to form dimers
  #  Top layer
  #  0.2440 -- 0.2381
  if { $surface_reconstruct } {
   set ny  [ MD++_Get latticesize(7) ]
   for { set i 0 } { $i < $ny } { incr i 1 } {
      set ymin [ expr -0.5006+1.0/$ny*$i ]
      set ymax [ expr -0.5006+1.0/$ny*($i+0.5) ]
      MD++ input = \[ 1 -10 10 $ymin $ymax 0.403   10 \]  fixatoms_by_position
      MD++ input = \[ 1 -10 10 $ymin $ymax  -10  -0.416 \]  fixatoms_by_position
   }
   MD++ input = 1  setfixedatomsgroup  freeallatoms
  
   for { set i 0 } { $i < $ny } { incr i 1 } {
      set ymin [ expr -0.5006+1.0/$ny*($i+0.5) ]
      set ymax [ expr -0.5006+1.0/$ny*($i+1) ]
      MD++ input = \[ 1 -10 10 $ymin $ymax 0.403   10 \]  fixatoms_by_position
      MD++ input = \[ 1 -10 10 $ymin $ymax -10 -0.416 \]  fixatoms_by_position
   }
   MD++ input = 2  setfixedatomsgroup  freeallatoms
  
   MD++ input = \[ 1  0  0.8 0  1 \] movegroup  
   MD++ input = \[ 1  0 -0.8 0  2 \] movegroup  
  }
  #end of surf reconstruct

  relax_fixbox
  MD++ finalcnfile = "0K_0.0_relaxed_surf${flag}.cn" writecn
  MD++ finalcnfile = "0K_0.0_relaxed_surf${flag}.cfg" writeatomeyecfg

  set epsilon $opt3
  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]
  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR
  MD++ eval
  set fp [ open "EPOT_1.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
  MD++ finalcnfile = "0K_${epsilon}_strained.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_strained.cfg" writeatomeyecfg

  set NP [MD++_Get NP]
  MD++ refreshnnlist
  MD++ input = $NP fprintnnlist

  exitmd


  #This is same as status ==0 but creates a surface pit structure
} elseif { $status == 500 } {
  MD++ setnolog
  initmd $n
  if { $argc > 4 } {
    MD++ dirname = $opt2
  } else { 
    puts "Need to specify dirname from python"
    MD++ quit
  }
 
  readpot
  make_perfect_crystal 12 12 21
  MD++ eval saveH conj_fixboxvec =  \[ 0 1 0  0 0 0   1 1 0 \]
  MD++ conj_fixbox = 1 conj_ftol = 2e-2 conj_itmax = 1000 conj_fevalmax = 2000
  MD++ relax finalcnfile = "0K_perfect.cn" writecn

  MD++ eval
# setup_window
# openwindow


  set dighole 0
  set make_free_surface 1  
  set surface_reconstruct 1
  
  if { $dighole } {
     MD++ input = \[ 1 -0.25 0.25  -0.25 0.25   0.4 10 \] fixatoms_by_position
     MD++ removefixedatoms
     MD++ freeallatoms      
  }
   #end of dig hole

  if {$make_free_surface} {
    MD++ vacuumratio = [expr 1.0-1.0/(1.0 + 0.2)]
    MD++ conj_fixbox = 0 
    if { $flag == 001 } {
       # create 001 surface
       MD++ input = \[ 3 3 0.6 \] changeH_keepR relax relax finalcnfile = "0K_surf${flag}.cn" writecn
       set strain_data_file "strain_surf_001.dat"
       set stress_data_file "stress_surf_001.dat"
       MD++ { conj_fixboxvec = [ 0 1 1
                                 1 1 1
                                 1 1 1 ] }
    }
    MD++ conj_fixbox = 1
  }
   # end of make free surface

  set NP [ MD++_Get "NP" ]
  set xholemin [ expr -0.5 + 0.2083 ]
  set yholemin [ expr -0.5 + 0.2083 ]
  set xholemax [ expr 0.5 - 0.2083 ]
  set yholemax [ expr 0.5 - 0.2083 ]
  set zholemax [ expr 0.5*5.0/6.0 - 0.14 ]
  
  for { set ip 0 } { $ip < $NP } { incr ip 1 } {
       set sxi [ MD++_GetVector SR $ip x ] 
       set syi [ MD++_GetVector SR $ip y ]
       set szi [ MD++_GetVector SR $ip z ]
       if { $sxi < $xholemax && $sxi > $xholemin && $syi < $yholemax && $syi > $yholemin && $szi > $zholemax } {
          MD++ fixed($ip) = 1
       }
  }
  
  MD++ removefixedatoms
  MD++ freeallatoms
  
  #setup_window
  #openwindow
  #MD++ plot
  #MD++ sleep 

  #  move top and bot layer atoms to form dimers
  #  Top layer

  #  0.2440 -- 0.2381


  set zholemin_tol [expr $zholemax-0.02]
  set zholemax_tol [expr $zholemax+0.02]
  if { $surface_reconstruct } {
   set ny  [ MD++_Get latticesize(7) ]
   for { set i 0 } { $i < $ny } { incr i 1 } {
      set ymin [ expr -0.5006+1.0/$ny*$i ]
      set ymax [ expr -0.5006+1.0/$ny*($i+0.5) ]
      MD++ input= \[ 1 -10 10 $ymin $ymax 0.3025   10 \]  fixatoms_by_position
      MD++ input= \[ 1 -10 10 $ymin $ymax  -10  -0.310 \]  fixatoms_by_position
      MD++ input = \[ 1 -10 10 $ymin $ymax  $zholemin_tol  $zholemax_tol \]  fixatoms_by_position
   }
   MD++ input = 1  setfixedatomsgroup  freeallatoms
  
   for { set i 0 } { $i < $ny } { incr i 1 } {
      set ymin [ expr -0.5006+1.0/$ny*($i+0.5) ]
      set ymax [ expr -0.5006+1.0/$ny*($i+1) ]
      MD++ input = \[ 1 -10 10 $ymin $ymax 0.3025   10 \]  fixatoms_by_position
      MD++ input = \[ 1 -10 10 $ymin $ymax -10 -0.310 \]  fixatoms_by_position
      MD++ input = \[ 1 -10 10 $ymin $ymax $zholemin_tol $zholemax_tol \]  fixatoms_by_position
   }
   MD++ input = 2  setfixedatomsgroup  freeallatoms
  
   MD++ input = \[ 1  0  0.8 0  1 \] movegroup  
   MD++ input = \[ 1  0 -0.8 0  2 \] movegroup  
  }
  #end of surf reconstruct

  relax_fixbox
  MD++ finalcnfile = "0K_0.0_relaxed_surf${flag}.cn" writecn
  MD++ finalcnfile = "0K_0.0_relaxed_surf${flag}.cfg" writeatomeyecfg

  set epsilon $opt3
  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]
  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR
  MD++ eval
  set fp [ open "EPOT_1.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
  MD++ finalcnfile = "0K_${epsilon}_strained.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_strained.cfg" writeatomeyecfg

  set NP [MD++_Get NP]
  MD++ refreshnnlist
  MD++ input = $NP fprintnnlist

  exitmd


# Read 0K_0.0_relaxed_surf${flag}.cn, Read nucleus-0.dat. 
# Move neighboring atoms of nucleus and apply strain to close trench.
# Relax. Evaluate energy of the frank partial. Write to EPOT_2.dat

} elseif { $status == 1 } {
  puts "n = $n"
  puts "flag = $flag"
  puts "opt1 = $opt1"
  puts "opt2 = $opt2"
  puts "opt3 = $opt3"
  puts "opt4 = $opt4"
  puts "opt5 = $opt5"
  puts "opt6 = $opt6"
  puts "opt7 = $opt7"
  puts "opt8 = $opt8"

  #MD++ setnolog
  initmd $status
  if { $argc > 9 } {
    MD++ dirname = $opt7
    exec cp ../../../scripts/work/frankMEP/disl_nuc_hetero.tcl .
    exec cp ../../../scripts/work/frankMEP/global_search.py .
  } else { 
    puts "Need to specify dirname from python"
    MD++ quit
  }
  set sliceid 1

  set A $flag
  set B $opt1
  set C $opt2
  set D $opt3

  set x0 $opt4
  set y0 $opt5
  set z0 $opt6

  # Read in structure with a slice of atoms removed
  MD++ incnfile = "0K_0.0_relaxed_surf001.cn" readcn relax
  set NP [ MD++_Get NP ] 

  # If atom is in tol of a nucleus atom, move this atom
  set tol [expr 0.05]

  # Read in nucleus to be removed
  # nbratoms.dat is atoms near the slice to speed up search
  set fpnbr [ open "nbratoms.dat" r ]
  set fp    [ open "nucleus-$n.dat" r ]
  set file_data [ read $fp ]
  set nbr_data  [ read $fpnbr ]
  close $fp
  close $fpnbr
  set data [split $file_data "\n"];
  set nbrdata [split $nbr_data "\n"];

  # Move atoms on each side of the removed plane
  set heterogeneous 1
  if { $heterogeneous == 1 } { 
    set Nnuc [ llength $data ]
    # the last line is empty
    set Nnuc [ expr $Nnuc -1 ]
    set Nnbr [ llength $nbrdata ]
    # the last line is empty
    set Nnbr [ expr $Nnbr -1 ] 
    puts "Nnuc = $Nnuc"
    puts "Nnbr = $Nnbr"

    set cnt1 0
    foreach line_nbr $nbrdata { 
      set nbrlinedata [ split $line_nbr " " ]
      lassign $nbrlinedata i xi yi zi
      if { $cnt1 < $Nnbr } {
        set cnt2 0
        foreach line $data {
          if { $cnt2 < $Nnuc } { 
            set linedata [ split $line " " ]
            lassign $linedata j x y z
            set dx [ expr $xi - $x ]
            set dy [ expr $yi - $y ]
            set dz [ expr $zi - $z ]
            if { [expr abs($dx) ] < $tol && [expr abs($dy) ] < $tol && [expr abs($dz)]<$tol } {
               set tmp [expr $A * $xi + $B * $yi + $C * $zi - $D ]
               if { $tmp >0 } {
	         set kk [ expr $i ]
		 MD++ fixed($kk) = 1
               } 
            }
          }
   	  set cnt2 [expr $cnt2 + 1]
        }
      }
      set cnt1 [expr $cnt1 + 1 ]
    }
    MD++ input = 1 setfixedatomsgroup freeallatoms

    set cnt1 0
    foreach line_nbr $nbrdata { 
      set nbrlinedata [ split $line_nbr " " ]
      lassign $nbrlinedata i xi yi zi
      if { $cnt1 < $Nnbr } {
        set cnt2 0
        foreach line $data {
          if { $cnt2 < $Nnuc } { 
            set linedata [ split $line " " ]
            lassign $linedata j x y z
            set dx [ expr $xi - $x ]
            set dy [ expr $yi - $y ]
            set dz [ expr $zi - $z ]
            if { [expr abs($dx) ] < $tol && [expr abs($dy) ] < $tol && [expr abs($dz)]<$tol } {
               set tmp [expr $A * $xi + $B * $yi + $C * $zi - $D ]
               if { $tmp <0 } {
	         set kk [ expr $i ]
		 MD++ fixed($kk) = 1
               }
            }
          }
	  set cnt2 [expr $cnt2 + 1]
        }
      }
      set cnt1 [expr $cnt1 + 1 ]
    }
    MD++ input = 2 setfixedatomsgroup freeallatoms

if { 0 } {  
   MD++ input = \[ 1 1 \] fixatoms_by_group
   MD++ input = \[ 1 2 \] fixatoms_by_group
   set cnt 0
   for { set i 0 } { $i < $NP } { incr i 1 } {
     set tmp [ MD++_Get fixed($i) ]
     if { $tmp == 0 } { 
       MD++ fixed($i) = 1 
     } else { 
       MD++ fixed($i) =  0
       set cnt [ expr $cnt + 1 ]
     }
   }
  MD++ removefixedatoms
  MD++ finalcnfile = "tmp.cn" writecn
  MD++ finalcnfile = "tmp.cfg" writeatomeyecfg
  puts "cnt = $cnt"
  setup_window
  openwindow
  MD++ plot
  MD++ sleepseconds = 10 sleep 
}
    set mag  0.9
    set magx [ expr $A*$mag ]
    set magy [ expr $B*$mag ]
    set magz [ expr $C*$mag ]

    puts "Move neighboring atoms closer"
    MD++ input = \[ 1 -$magx -$magy -$magz 1 \] movegroup
    MD++ input = \[ 1  $magx  $magy  $magz 2 \] movegroup

    puts "Remove atoms belonging to nucleus"
    set cnt2 0
    foreach line $data {
      if { $cnt2 < $Nnuc } { 
        set linedata [ split $line " " ]
        lassign $linedata j x y z
	MD++ fixed($j) = 1
      }
      set cnt2 [ expr $cnt2 + 1 ]
    }
    MD++ removefixedatoms freeallatoms
  }

  puts "Now apply strain and close the trench"
  set epsilon $opt8
  MD++ finalcnfile = "0K_${epsilon}_trench_${n}.cfg" writeatomeyecfg

  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]

  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR

  MD++ finalcnfile = "0K_${epsilon}_relaxed_${n}.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_relaxed_${n}.cfg" writeatomeyecfg
  MD++ eval 
  set fp [ open "EPOT_2.dat" w ]; puts $fp [ MD++_Get EPOT ]; close $fp

  exitmd

# This status is the same as status =1, but use topological information to apply heterogeneous strain
# That is, move the whole planes near by. 
} elseif { $status == 11 } {
  puts "n = $n"
  puts "flag = $flag"
  puts "opt1 = $opt1"
  puts "opt2 = $opt2"
  puts "opt3 = $opt3"
  puts "opt4 = $opt4"
  puts "opt5 = $opt5"
  puts "opt6 = $opt6"
  puts "opt7 = $opt7"
  puts "opt8 = $opt8"

#  MD++ setnolog
  initmd $status
  if { $argc > 9 } {
    MD++ dirname = $opt7
    exec cp ../../../scripts/work/frankMEP/disl_nuc_hetero.tcl .
    exec cp ../../../scripts/work/frankMEP/global_search.py .
  } else { 
    puts "Need to specify dirname from python"
    MD++ quit
  }
  set sliceid 1

  set A $flag
  set B $opt1
  set C $opt2
  set D $opt3

  set x0 $opt4
  set y0 $opt5
  set z0 $opt6

  # Read in structure with a slice of atoms removed
  MD++ incnfile = "0K_0.0_relaxed_surf001.cn" readcn relax
  set NP [ MD++_Get NP ] 

  # If atom is in tol of a nucleus atom, move this atom
  set tol [expr 0.05]

  # Read in nucleus to be removed
  # nbratoms.dat is atoms near the slice to speed up search
  set fpnbr_u [ open "nucleus_nbrlist_u-$n.dat" r ]
  set fpnbr_d [ open "nucleus_nbrlist_d-$n.dat" r ]
  set fp    [ open "nucleus-$n.dat" r ]
  set file_data [ read $fp ]
  set nbr_data_u  [ read $fpnbr_u ]
  set nbr_data_d  [ read $fpnbr_d ]
  close $fp
  close $fpnbr_u
  close $fpnbr_d
  set data [split $file_data "\n"];
  set nbrdata_u [split $nbr_data_u "\n"];
  set nbrdata_d [split $nbr_data_d "\n"];
  set Nnuc [ llength $data ]
# the last line is empty
  set Nnuc [ expr $Nnuc -1 ] 

  # Move atoms on each side of the removed plane
  set heterogeneous 1
  if { $heterogeneous == 1 } { 
    set Nnbr_u [ llength $nbrdata_u ]
    set Nnbr_d [ llength $nbrdata_d ]
    # the last line is empty
    set Nnbr_u [ expr $Nnbr_u -1 ] 
    set Nnbr_d [ expr $Nnbr_d -1 ] 
    puts "Nnbr_u = $Nnbr_u"
    puts "Nnbr_d = $Nnbr_d"

    set cnt1 0
    foreach line_nbr $nbrdata_u { 
      if { $cnt1 < $Nnbr_u } {
        set nbrlinedata [ split $line_nbr " " ]
        lassign $nbrlinedata i xi yi zi
        MD++ fixed($i) = 1
        set cnt1 [expr $cnt1 + 1 ]
      }
    }
    MD++ input = 1 setfixedatomsgroup freeallatoms

    set cnt1 0
    foreach line_nbr $nbrdata_d { 
      if { $cnt1 < $Nnbr_d } {
        set nbrlinedata [ split $line_nbr " " ]
        lassign $nbrlinedata i xi yi zi
        MD++ fixed($i) = 1
        set cnt1 [expr $cnt1 + 1 ]
      }
    }
    MD++ input = 2 setfixedatomsgroup freeallatoms

if { 0 } {  
   MD++ input = \[ 1 1 \] fixatoms_by_group
   MD++ input = \[ 1 2 \] fixatoms_by_group
   set cnt 0
   for { set i 0 } { $i < $NP } { incr i 1 } {
     set tmp [ MD++_Get fixed($i) ]
     if { $tmp == 0 } { 
       MD++ fixed($i) = 1 
     } else { 
       MD++ fixed($i) =  0
       set cnt [ expr $cnt + 1 ]
     }
   }
  MD++ removefixedatoms
  MD++ finalcnfile = "tmp.cn" writecn
  MD++ finalcnfile = "tmp.cfg" writeatomeyecfg
  puts "cnt = $cnt"
  setup_window
  openwindow
  MD++ plot
  MD++ sleepseconds = 10 sleep 
}
    set mag  0.9
    set magx [ expr $A*$mag ]
    set magy [ expr $B*$mag ]
    set magz [ expr $C*$mag ]

    puts "Move neighboring atoms closer"
    MD++ input = \[ 1 -$magx -$magy -$magz 1 \] movegroup
    MD++ input = \[ 1  $magx  $magy  $magz 2 \] movegroup

    puts "Remove atoms belonging to nucleus"
    set cnt2 0
    foreach line $data {
      if { $cnt2 < $Nnuc } { 
        set linedata [ split $line " " ]
        lassign $linedata j x y z
	MD++ fixed($j) = 1
      }
      set cnt2 [ expr $cnt2 + 1 ]
    }
    MD++ removefixedatoms freeallatoms
  }

  puts "Now apply strain and close the trench"
  set epsilon $opt8
  MD++ finalcnfile = "0K_${epsilon}_trench_${n}.cfg" writeatomeyecfg

  set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]

  MD++ saveH
  set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
  MD++ H_11 = $H11_fix
  MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
  MD++ conj_fixbox = 1  relax
  MD++ SHtoR

  MD++ finalcnfile = "0K_${epsilon}_relaxed_${n}.cn" writecn
  MD++ finalcnfile = "0K_${epsilon}_relaxed_${n}.cfg" writeatomeyecfg
  MD++ eval 
  set fp [ open "EPOT_2.dat" w ]; puts $fp [ MD++_Get EPOT ]; close $fp

  exitmd

} elseif { $status == 2 } { 
  initmd $status
  MD++ setnolog
  MD++ incnfile = "runs/frankMEP/db_0.04/0K_0.0_relaxed_surf001.cn" readcn
  #make_perfect_crystal 12 12 16
  MD++ relax eval
  set EPOT_0 [MD++_Get EPOT]
  puts "EPOT_0 = $EPOT_0"

  set Get_strain_dependent_cohsv 1 
  if { $Get_strain_dependent_cohsv == 1 } { 
    set EPOT_i [MD++_Get EPOT_IND(1000)]
    puts "Unstrained EPOT_i = $EPOT_i"

    set epsilon 0.04
    set H11_0 [ MD++_Get H_11 ]; set H22_0 [ MD++_Get H_22 ]; set H33_0 [ MD++_Get H_33 ]

    MD++ saveH
    set H11_fix [ expr $H11_0*(1.0-$epsilon) ]
    MD++ H_11 = $H11_fix
    MD++ conj_ftol = 1e-4 conj_itmax = 3800 conj_fevalmax = 6000
    MD++ conj_fixbox = 1  relax
    MD++ SHtoR
    set EPOT_i [MD++_Get EPOT_IND(10)]
    puts "Strained EPOT_i = $EPOT_i"
    MD++ sleep
  }


  set NP [MD++_Get NP]
  set n [myRand 0 $NP]
  MD++ input=\[ 1 $n \] fixatoms_by_ID
  MD++ removefixedatoms
  MD++ freeallatoms
  MD++ relax
  MD++ eval 
  set EPOT_1 [MD++_Get EPOT]
  puts "EPOT_1 = $EPOT_1"
} elseif { $status == 100 } {
  MD++ incnfile = "${n}.cn" readcn
  MD++ finalcnfile = "${n}.cfg" writeatomeyecfg
} elseif { $status == 1000 } {
  #for { set i 0 } { $i < 20 } { incr i } { 
    MD++ incnfile = "runs/frankMEP/temp/0K_0.05_strained.cn" readcn
    MD++ eval
    set fp [ open "EPOT.dat" a+ ]; puts $fp [ MD++_Get EPOT ]; close $fp
  #}
  #MD++ incnfile = "runs/frankMEP/0K_0.04_relaxed_0.cn" readcn
  MD++ eval
} else {
 puts "unknown status = $status"
 exitmd 
} 
