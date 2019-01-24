source "/scratch/users/hfernandes/BioSIM-repository/energySplit/energySplit/externalLib/linalg.tcl"
package require math::linearalgebra

global colori colorf
set colori "\033\[1;33m"
set colorf "\033\[0m"

proc bond {inputFile} {
    global colori colorf
    # Read the inputFile. The following variables are set: numberAtoms pi types charges parameters connectivity xyz fragList
    source "$inputFile"

    #Get HrmStr1
    set bondParam {}
    foreach line $parameters {
        set var [regexp -inline {HrmStr1\s+"(\S+)"\s+"(\S+)"\s+(\S+)\s+(\S+)} $line]
        if {$var == ""} {
            continue
        } else {
            lappend bondParam [lrange $var 1 4]
        }
    }

    set bondEnergyTable {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        # Get all bonds involving the atom $index
        #set bonds [lsort -increasing -unique -real [join [lsearch -all -inline -regexp $connectivity "(^| )${index}($| )"]]]
        set bonds [lindex $connectivity $index]

        # Get atomtype for the atom $index
        set atomType [lindex $types $index]

        foreach bond $bonds {
            if {$bond > $index} {
                set parameter ""

                # Get atomtype for the atom $bond
                set atomType1 [lindex $types $bond]

                # Getting atom coordinates
                set xi [lindex [lindex $xyz $index] 0]
                set yi [lindex [lindex $xyz $index] 1]
                set zi [lindex [lindex $xyz $index] 2]
                set xj [lindex [lindex $xyz $bond] 0]
                set yj [lindex [lindex $xyz $bond] 1]
                set zj [lindex [lindex $xyz $bond] 2]

                # Calculating the distance between the atoms $index and $bond
                set l [expr {( ($xi - $xj)**2 + ($yi - $yj)**2 + ($zi - $zj)**2 )**(0.5)}]

                # Get bond parameters
                set parameter [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $bondParam $atomType] $atomType1]
                if {$parameter == ""} {
                    set parameter [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $bondParam $atomType1] $atomType]
                }

                # Calculate the energy associated to a certain bond
                if {$parameter != ""} {
                    set l0 [lindex [lindex $parameter 0] 3]
                    set kb [lindex [lindex $parameter 0] 2]

                    set energy [expr {((($l - $l0)**2) * $kb) / 627.50947415151515}]

                    # Check where the atoms $index and $bond belong
                    set logic0 [lsearch -all -regexp $fragList "(^| )${index}($| )"]
                    set logic1 [lsearch -all -regexp $fragList "(^| )${bond}($| )"]

                    if {$logic0 == $logic1} {
                        # Atoms belonging to the same fragment
                        if {$logic0 == ""} {
                            # Both atoms belonging to the remaining system
                            if {[lsearch -index 0 $bondEnergyTable "none"] == -1} {
                                lappend bondEnergyTable [list "none" 0]
                            }
                            set pos [lsearch -index 0 $bondEnergyTable "none"]
                            set prevEnergy [lindex [lindex $bondEnergyTable $pos] 1]
                            set finalEnergy [expr $prevEnergy + $energy]
                            set bondEnergyTable [lreplace $bondEnergyTable $pos $pos [list "none" $finalEnergy]]
                        } else {
                            # Atoms belonging to a certain fragment
                            if {[lsearch -index 0 $bondEnergyTable "$logic0-$logic1"] == -1} {
                                 lappend bondEnergyTable [list "$logic0-$logic1" 0]
                            }
                            set pos [lsearch -index 0 $bondEnergyTable "$logic0-$logic1"]
                            set prevEnergy [lindex [lindex $bondEnergyTable $pos] 1]
                            set finalEnergy [expr $prevEnergy + $energy]
                            set bondEnergyTable [lreplace $bondEnergyTable $pos $pos [list "$logic0-$logic1" $finalEnergy]]
                        }
                    } else {
                        # Atoms belonging to different fragments
                        if {[lsearch -index 0 -exact $bondEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $bondEnergyTable "$logic1-$logic0"] == -1} {
                            lappend bondEnergyTable [list "$logic0-$logic1" 0]
                        }
                        set pos [lsearch -index 0 -exact $bondEnergyTable "$logic0-$logic1"]
                        if {$pos == -1} {
                            set pos [lsearch -index 0 -exact $bondEnergyTable "$logic1-$logic0"]
                        }
                        set prevEnergy [lindex [lindex $bondEnergyTable $pos] 1]
                        set finalEnergy [expr $prevEnergy + $energy]
                        set bondEnergyTable [lreplace $bondEnergyTable $pos $pos [list "$logic0-$logic1" $finalEnergy]]
                    }

                   
                } else {
                    puts "Missing parameter for bond $atomType-$atomType1 (indexes $index $bond)"
                }

            }
        }
    }

    # Calculating the final energy for bonds
    puts "####################\n## $colori Bonds $colorf\n####################"
    puts "\tEnergy Decomposition:"
    set fe 0
    set bondEnergyTable [lsort -index 0 $bondEnergyTable]
    foreach line $bondEnergyTable {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Remaining system\t"
        } else {
            set a [split [lindex $line 0] "-"]
            if {[lindex $a 0] == [lindex $a 1]} {
                set description "Fragment [lindex $a 0]\t\t"
            } elseif {[lindex $a 0] == ""} {
                set description "Fragment [lindex $a 1] + Remaining system"
            } elseif {[lindex $a 1] == ""} {
                set description "Fragment [lindex $a 0] + Remaining system"
            } else {
                set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t"
            }
        }
        puts "\t\t$description\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"

}

proc angleDihedral {inputFile} {
    global colori colorf
    # Read the inputFile. The following variables are set: numberAtoms pi types charges parameters connectivity xyz fragList
    source "$inputFile"

    set uniqueAngles {}
    set uniqueDihedral {}
    set angleEnergyTable {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set atoms1 [lindex $connectivity $index]

        foreach atom $atoms1 {
            set atoms2 [lindex $connectivity $atom]

            foreach atom1 $atoms2 {
                if {$atom1 != $index} {
                    set atoms3 [lindex $connectivity $atom1]

                    # Dihedral
                    foreach atom2 $atoms3 {
                        if {$atom2 != $index && $atom2 != $atom} {
                            set newLine [list $index $atom $atom1 $atom2]
                            set newLine1 [list $atom2 $atom1 $atom $index]

                            if {[lsearch $uniqueDihedral $newLine] == -1 && [lsearch $uniqueDihedral $newLine1] == -1} {
                                lappend uniqueDihedral $newLine
                            }
                        }
                    }

                    # Angles
                    set newLine [list $index $atom $atom1]
                    set newLine1 [list $atom1 $atom $index]

                    if {[lsearch $uniqueAngles $newLine] == -1 && [lsearch $uniqueAngles $newLine1] == -1} {
                        lappend uniqueAngles $newLine
                    }
                }
            }
        }
    }

    ###### Angles

    ## Get HrmBnd1
    set angleParam {}
    foreach line $parameters {
        set var [regexp -inline {HrmBnd1\s+"(\S+)"\s+"(\S+)"\s+"(\S+)"\s+(\S+)\s+(\S+)} $line]
        if {$var == ""} {
            continue
        } else {
            lappend angleParam [lrange $var 1 5]
        }
    }

    # Calculate energy
    set angleEnergy 0
    foreach angle $uniqueAngles {
        set parameter ""
        set index0 [lindex $angle 0]
        set index1 [lindex $angle 1]
        set index2 [lindex $angle 2]

        set type0 [lindex $types $index0]
        set type1 [lindex $types $index1]
        set type2 [lindex $types $index2]

        set parameter [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $angleParam $type0] $type1] $type2]
        if {$parameter == ""} {
            set parameter [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $angleParam $type2] $type1] $type0]
        }

        if {$parameter != ""} {
            # Getting atom coordinates
            set xa [lindex [lindex $xyz $index0] 0]
            set ya [lindex [lindex $xyz $index0] 1]
            set za [lindex [lindex $xyz $index0] 2]
            set xb [lindex [lindex $xyz $index1] 0]
            set yb [lindex [lindex $xyz $index1] 1]
            set zb [lindex [lindex $xyz $index1] 2]
            set xc [lindex [lindex $xyz $index2] 0]
            set yc [lindex [lindex $xyz $index2] 1]
            set zc [lindex [lindex $xyz $index2] 2]

            # Calculate vectors
            set xyza [list $xa $ya $za]
            set xyzb [list $xb $yb $zb]
            set xyzc [list $xc $yc $zc]
            set vec0 [::math::linearalgebra::sub $xyza $xyzb]
            set vec1 [::math::linearalgebra::sub $xyzc $xyzb]

            # Calculate Angle
            set teta [expr [::math::linearalgebra::angle $vec0 $vec1] * (180/$pi)]

            set kangle [lindex [lindex $parameter 0] 3]
            set teta0 [lindex [lindex $parameter 0] 4]

            # Calculate energy
            set energy [expr {(((($teta - $teta0)*($pi/180))**2) * $kangle ) / 627.50947415151515}]

            # Check where the atoms $index0 $index1 and $index2 belong
            set logic0 [lsearch -all -regexp $fragList "(^| )${index0}($| )"]
            set logic1 [lsearch -all -regexp $fragList "(^| )${index1}($| )"]
            set logic2 [lsearch -all -regexp $fragList "(^| )${index2}($| )"]

            if {$logic0 == $logic1 && $logic1 == $logic2} {
                # Atoms belonging to the same fragment
                if {$logic0 == ""} {
                    # All the atoms belonging to the remaining system
                    if {[lsearch -index 0 $angleEnergyTable "none"] == -1} {
                        lappend angleEnergyTable [list "none" 0]
                    }
                    set pos [lsearch -index 0 $angleEnergyTable "none"]
                    set prevEnergy [lindex [lindex $angleEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set angleEnergyTable [lreplace $angleEnergyTable $pos $pos [list "none" $finalEnergy]]
                } else {
                    # All the atoms belong to a certain fragment
                    if {[lsearch -index 0 -exact $angleEnergyTable "$logic0"] == -1} {
                        lappend angleEnergyTable [list "$logic0" 0]
                    }
                    set pos [lsearch -index 0 $angleEnergyTable "$logic0"]
                    set prevEnergy [lindex [lindex $angleEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set angleEnergyTable [lreplace $angleEnergyTable $pos $pos [list "$logic0" $finalEnergy]]
                }
            } else {
                # Atoms belong to different fragments
                if {$logic0 == $logic1} {
                    if {[lsearch -index 0 -exact $angleEnergyTable "$logic0-$logic2"] == -1 && [lsearch -index 0 -exact $angleEnergyTable "$logic2-$logic0"] == -1 } {
                        lappend angleEnergyTable [list "$logic0-$logic2" 0]
                    }
                    set pos [lsearch -index 0 -exact $angleEnergyTable "$logic0-$logic2"]  
                    if {$pos == -1} {
                        set pos [lsearch -index 0 -exact $angleEnergyTable "$logic2-$logic0"]
                    }
                    set prevEnergy [lindex [lindex $angleEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set angleEnergyTable [lreplace $angleEnergyTable $pos $pos [list "$logic0-$logic2" $finalEnergy]]
                } elseif {$logic0 == $logic2} {
                    if {[lsearch -index 0 -exact $angleEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $angleEnergyTable "$logic1-$logic0"] == -1 } {
                        lappend angleEnergyTable [list "$logic0-$logic1" 0]
                    }
                    set pos [lsearch -index 0 -exact $angleEnergyTable "$logic0-$logic1"]  
                    if {$pos == -1} {
                        set pos [lsearch -index 0 -exact $angleEnergyTable "$logic1-$logic0"]
                    }
                    set prevEnergy [lindex [lindex $angleEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set angleEnergyTable [lreplace $angleEnergyTable $pos $pos [list "$logic0-$logic1" $finalEnergy]]
                } elseif {$logic1 == $logic2} {
                    if {[lsearch -index 0 -exact $angleEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $angleEnergyTable "$logic1-$logic0"] == -1 } {
                        lappend angleEnergyTable [list "$logic0-$logic1" 0]
                    }
                    set pos [lsearch -index 0 -exact $angleEnergyTable "$logic0-$logic1"]  
                    if {$pos == -1} {
                        set pos [lsearch -index 0 -exact $angleEnergyTable "$logic1-$logic0"]
                    }
                    set prevEnergy [lindex [lindex $angleEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set angleEnergyTable [lreplace $angleEnergyTable $pos $pos [list "$logic0-$logic1" $finalEnergy]]
                } else {
                    set logicList [lsort [list $logic0 $logic1 $logic2]]
                    if {[lsearch -index 0 -exact $angleEnergyTable "[lindex $logicList 0]-[lindex $logicList 1]-[lindex $logicList 2]"] == -1} {
                        lappend angleEnergyTable [list "[lindex $logicList 0]-[lindex $logicList 1]-[lindex $logicList 2]" 0]
                    }
                    set pos [lsearch -index 0 -exact $angleEnergyTable "[lindex $logicList 0]-[lindex $logicList 1]-[lindex $logicList 2]"]
                    set prevEnergy [lindex [lindex $angleEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set angleEnergyTable [lreplace $angleEnergyTable $pos $pos [list "[lindex $logicList 0]-[lindex $logicList 1]-[lindex $logicList 2]" $finalEnergy]]
                }
            }
   
        } else {
            puts "Missing parameter for angle $type0-$type1-$type2 (indexes $index0 $index1 $index2)"
        }

    }

    # Calculating the final energy for angles
    puts "####################\n##$colori Angles $colorf\n####################"
    puts "\tEnergy Decomposition:"
    set fe 0
    set angleEnergyTable [lsort -index 0 $angleEnergyTable]
    foreach line $angleEnergyTable {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Remaining system\t\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            set aLength [llength $a]

            if {$aLength == 1} {
                set description "Fragment [lindex $a 0]\t\t\t\t"
            } elseif {$aLength == 2} {
                if {[lindex $a 0] != "" && [lindex $a 1] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t\t\t"
                } else {
                    set b [lindex $a 0]
                    if {$b == ""} {
                        set b [lindex $a 1]
                    }
                    set description "Fragment $b + Remaining system\t\t"
                }
            } elseif {$aLength == 3} {
                if {[lindex $a 0] != "" && [lindex $a 1] != "" && [lindex $a 2] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1] + Fragment [lindex $a 2]\t\t"
                } else {
                    set b [lindex $a 0]
                    if {$b == ""} {
                        set b [lindex $a 1]
                        set c [lindex $a 2]
                    } else {
                        set c [lindex $a 1]
                        if {$c == ""} {
                            set c [lindex $a 2]
                        }
                    }
                    set description "Fragment $b + Fragment $c + Remaining system"
                }
            } else {
                puts "Something went wrong!"
            }

        }
        puts "\t\t$description\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"


    ##### Dihedral

    ## Get AmbTrs
    set dihedralParam {}
    foreach line $parameters {
        set var [regexp -inline {AmbTrs\s+"(\S+)"\s+"(\S+)"\s+"(\S+)"\s+"(\S+)"\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)} $line]
        if {$var == ""} {
            continue
        } else {
            lappend dihedralParam [lrange $var 1 13]
        }
    }

    # Calculate energy
    set dihedralEnergyFrag 0
    set dihedralEnergyInteraction 0
    set dihedralEnergyNotFrag 0
    foreach dihedral $uniqueDihedral {

        set parameter ""
        set index0 [lindex $dihedral 0]
        set index1 [lindex $dihedral 1]
        set index2 [lindex $dihedral 2]
        set index3 [lindex $dihedral 3]

        set type0 [lindex $types $index0]
        set type1 [lindex $types $index1]
        set type2 [lindex $types $index2]
        set type3 [lindex $types $index3]

        set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $dihedralParam $type0] $type1] $type2] $type3]
        if {$parameter == ""} {
            set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $dihedralParam $type3] $type2] $type1] $type0]
        }

        if {$parameter != ""} {
            set teta [measure dihed [list $index0 $index1 $index2 $index3]]
            set po1 [lindex [lindex $parameter 0] 4]
            set po2 [lindex [lindex $parameter 0] 5]
            set po3 [lindex [lindex $parameter 0] 6]
            set po4 [lindex [lindex $parameter 0] 7]
            set mag1 [lindex [lindex $parameter 0] 8]
            set mag2 [lindex [lindex $parameter 0] 9]
            set mag3 [lindex [lindex $parameter 0] 10]
            set mag4 [lindex [lindex $parameter 0] 11]
            set npaths [lindex [lindex $parameter 0] 12]

            set energyDihed [expr {(($mag1*(1+cos(((1*$teta-$po1)*($pi/180)))))/$npaths) + (($mag2*(1+cos(((2*$teta-$po2)*($pi/180)))))/$npaths) + (($mag3*(1+cos(((3*$teta-$po3)*($pi/180)))))/$npaths) + (($mag4*(1+cos(((4*$teta-$po4)*($pi/180)))))/$npaths)}]

            set logic0 [lsearch $fragList $index0]
            set logic1 [lsearch $fragList $index1]
            set logic2 [lsearch $fragList $index2]
            set logic3 [lsearch $fragList $index3]

            if {[llength [lsearch -all -inline [list $logic0 $logic1 $logic2 ] "-1"]] == 4} {
                set dihedralEnergyNotFrag [expr {$dihedralEnergyNotFrag + $energyDihed}]
            } elseif {[llength [lsearch -all -inline [list $logic0 $logic1 $logic2 $logic3] "-1"]] == 0} {
                set dihedralEnergyFrag [expr {$dihedralEnergyFrag + $energyDihed}]
            } else {
                set dihedralEnergyInteraction [expr {$dihedralEnergyInteraction + $energyDihed}]
            }
            
        } else {
            puts "Missing parameter for dihedral $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
        }

    }

}

### Start procedure

puts "################################################################################\n### Energy Split\n################################################################################\n### Output\n################################################################################\n"
puts "####################\n### $colori Fragments $colorf\n####################"
source $argv
set i 0
foreach line $fragList {
    puts "$colori\tFragment $i $colorf: $line"
    incr i
}
puts "\n"

bond $argv
angleDihedral $argv

puts "################################################################################\n### The calculation finished successfully on [clock format [clock seconds] -format %Y_%b_%d] at [clock format [clock seconds] -format %H:%M:%S]\n################################################################################"