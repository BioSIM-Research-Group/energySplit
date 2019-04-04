source "[file dirname [file normalize [info script]]]/externalLib/linalg.tcl"
package require math::linearalgebra
package require Thread

global colori colorf bondTotalList angleTotalList dihedralTotalList impTotalList vdwTotalList coloumbTotalList
set colori ""
set colorf ""

proc bond {inputFile outputFile} {
    global colori colorf bondTotalList
    # Read the inputFile. The following variables are set: numberAtoms pi types charges parameters connectivity xyz fragList
    source "$inputFile"

    puts $outputFile "####################\n## Bonds \n####################"
    puts "####################\n## $colori Bonds $colorf\n####################"
    
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
                            # Both atoms belonging to the Fragment X
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
    puts $outputFile "\tEnergy Decomposition:"
    puts "\tEnergy Decomposition:"
    set fe 0
    set bondEnergyTable [lsort -index 0 $bondEnergyTable]
    foreach line $bondEnergyTable {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            if {[lindex $a 0] == [lindex $a 1]} {
                set description "Fragment [lindex $a 0]\t\t\t"
            } elseif {[lindex $a 0] == ""} {
                set description "Fragment [lindex $a 1] + Fragment X\t\t"
            } elseif {[lindex $a 1] == ""} {
                set description "Fragment [lindex $a 0] + Fragment X\t\t"
            } else {
                set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t\t"
            }
        }
        puts "\t\t$description\t\t\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t\t\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

    # Export variable
    set bondTotalList $bondEnergyTable
}

proc angleDihedral {inputFile outputFile} {
    global colori colorf angleTotalList dihedralTotalList

    puts "####################\n##$colori Angles $colorf\n####################"
    puts $outputFile "####################\n## Angles \n####################"
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
                    # All the atoms belonging to the Fragment X
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
            puts $outputFile "Missing parameter for angle $type0-$type1-$type2 (indexes $index0 $index1 $index2)"
        }

    }

    # Calculating the final energy for angles
    puts "\tEnergy Decomposition:"
    puts $outputFile "\tEnergy Decomposition:"
    set fe 0
    set angleEnergyTable [lsort -index 0 $angleEnergyTable]
    foreach line $angleEnergyTable {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t\t\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            set aLength [llength $a]

            if {$aLength == 1} {
                set description "Fragment [lindex $a 0]\t\t\t\t\t"
            } elseif {$aLength == 2} {
                if {[lindex $a 0] != "" && [lindex $a 1] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t\t\t\t"
                } else {
                    set b [lindex $a 0]
                    if {$b == ""} {
                        set b [lindex $a 1]
                    }
                    set description "Fragment $b + Fragment X\t\t\t\t"
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
                    set description "Fragment $b + Fragment $c + Fragment X\t\t"
                }
            } else {
                puts "Something went wrong!"
                puts $outputFile "Something went wrong!"
            }

        }
        puts "\t\t$description\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

    # Export variable
    set angleTotalList $angleEnergyTable


    ##### Dihedral
    puts "####################\n##$colori Torsions $colorf\n####################"
    puts $outputFile "####################\n## Torsions \n####################"
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
    set dihedEnergyTable {}
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
            set xd [lindex [lindex $xyz $index3] 0]
            set yd [lindex [lindex $xyz $index3] 1]
            set zd [lindex [lindex $xyz $index3] 2]

            # Calculate vectors
            set xyza [list $xa $ya $za]
            set xyzb [list $xb $yb $zb]
            set xyzc [list $xc $yc $zc]
            set xyzd [list $xd $yd $zd]
            set vec0 [::math::linearalgebra::sub $xyza $xyzb]
            set vec1 [::math::linearalgebra::sub $xyzb $xyzc]
            set vec2 [::math::linearalgebra::sub $xyzc $xyzd]
            set vec3 [::math::linearalgebra::sub $xyzc $xyzb]

            set norm0 [::math::linearalgebra::crossproduct $vec0 $vec1]
            set norm1 [::math::linearalgebra::crossproduct $vec2 $vec3]

            set teta [expr [::math::linearalgebra::angle $norm0 $norm1] * (180/$pi)]
            set po1 [lindex [lindex $parameter 0] 4]
            set po2 [lindex [lindex $parameter 0] 5]
            set po3 [lindex [lindex $parameter 0] 6]
            set po4 [lindex [lindex $parameter 0] 7]
            set mag1 [lindex [lindex $parameter 0] 8]
            set mag2 [lindex [lindex $parameter 0] 9]
            set mag3 [lindex [lindex $parameter 0] 10]
            set mag4 [lindex [lindex $parameter 0] 11]
            set npaths [lindex [lindex $parameter 0] 12]

            set energy [expr {((($mag1*(1+cos(((1*$teta-$po1)*($pi/180)))))/$npaths) + (($mag2*(1+cos(((2*$teta-$po2)*($pi/180)))))/$npaths) + (($mag3*(1+cos(((3*$teta-$po3)*($pi/180)))))/$npaths) + (($mag4*(1+cos(((4*$teta-$po4)*($pi/180)))))/$npaths)) / 627.50947415151515}]            
            
            # Check where the atoms $index0 $index1 and $index2 belong
            set logic0 [lsearch -all -regexp $fragList "(^| )${index0}($| )"]
            set logic1 [lsearch -all -regexp $fragList "(^| )${index1}($| )"]
            set logic2 [lsearch -all -regexp $fragList "(^| )${index2}($| )"]
            set logic3 [lsearch -all -regexp $fragList "(^| )${index3}($| )"]

            if {$logic0 == $logic1 && $logic1 == $logic2 && $logic2 == $logic3} {
                # Atoms belong to the same fragmennt
                if {$logic0 == ""} {
                    # All the atoms belong to the Fragment X
                    if {[lsearch -index 0 $dihedEnergyTable "none"] == -1} {
                        lappend dihedEnergyTable [list "none" 0]
                    }
                    set pos [lsearch -index 0 $dihedEnergyTable "none"]
                    set prevEnergy [lindex [lindex $dihedEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set dihedEnergyTable [lreplace $dihedEnergyTable $pos $pos [list "none" $finalEnergy]]
                } else {
                    # All atoms belong to a certain fragment
                    if {[lsearch -index 0 -exact $dihedEnergyTable "$logic0"] == -1} {
                        lappend dihedEnergyTable [list "$logic0" 0]
                    }
                    set pos [lsearch -index 0 $dihedEnergyTable "$logic0"]
                    set prevEnergy [lindex [lindex $dihedEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set dihedEnergyTable [lreplace $dihedEnergyTable $pos $pos [list "$logic0" $finalEnergy]]
                }
            } else {
                # Atoms belong to different fragments
                set listLogic [lsort -unique [list $logic0 $logic1 $logic2 $logic3]]
                set listLogicLength [llength $listLogic]
                if {$listLogicLength == 2} {
                    if {[lsearch -index 0 -exact $dihedEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]"] == -1} {
                        lappend dihedEnergyTable [list "[lindex $listLogic 0]-[lindex $listLogic 1]" 0]
                    }
                    set pos [lsearch -index 0 -exact $dihedEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]"]  
                    set prevEnergy [lindex [lindex $dihedEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set dihedEnergyTable [lreplace $dihedEnergyTable $pos $pos [list "[lindex $listLogic 0]-[lindex $listLogic 1]" $finalEnergy]]
                } elseif {$listLogicLength == 3} {
                    if {[lsearch -index 0 -exact $dihedEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]"] == -1} {
                        lappend dihedEnergyTable [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]" 0]
                    }
                    set pos [lsearch -index 0 -exact $dihedEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]"]  
                    set prevEnergy [lindex [lindex $dihedEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set dihedEnergyTable [lreplace $dihedEnergyTable $pos $pos [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]" $finalEnergy]]
                } elseif {$listLogicLength == 4} {
                    if {[lsearch -index 0 -exact $dihedEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]"] == -1} {
                        lappend dihedEnergyTable [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]" 0]
                    }
                    set pos [lsearch -index 0 -exact $dihedEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]"]  
                    set prevEnergy [lindex [lindex $dihedEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set dihedEnergyTable [lreplace $dihedEnergyTable $pos $pos [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]" $finalEnergy]]
                } else {
                    puts "Something went wrong."
                }
            }

        } else {
            puts "Missing parameter for dihedral $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
            puts $outputFile "Missing parameter for dihedral $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
        }

    }

    # Calculating the final energy for torsions
    puts "\tEnergy Decomposition:"
    puts $outputFile "\tEnergy Decomposition:"
    set fe 0
    set dihedEnergyTable [lsort -index 0 $dihedEnergyTable]
    foreach line $dihedEnergyTable {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t\t\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            set aLength [llength $a]

            if {$aLength == 1} {
                set description "Fragment [lindex $a 0]\t\t\t\t\t"
            } elseif {$aLength == 2} {
                if {[lindex $a 0] != "" && [lindex $a 1] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t\t\t\t"
                } else {
                    set b [lindex $a 0]
                    if {$b == ""} {
                        set b [lindex $a 1]
                    }
                    set description "Fragment $b + Fragment X\t\t\t\t"
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
                    set description "Fragment $b + Fragment $c + Fragment X\t\t"
                }
            } elseif {$aLength == 4} {
                if {[lindex $a 0] != "" && [lindex $a 1] != "" && [lindex $a 2] != "" && [lindex $a 3] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1] + Fragment [lindex $a 2] + Fragment [lindex $a 3]"
                } else {
                    set b [lindex $a 0]
                    set c [lindex $a 1]
                    set d [lindex $a 2]
                    set e [lindex $a 3]

                    if {$b == ""} {
                        set b "X"
                    } elseif {$c == ""} {
                        set c "X"
                    } elseif {$d == ""} {
                        set d "X"
                    } elseif {$e == ""} {
                        set e "X"
                    } else {
                        puts "Something went wrong!"
                        puts $outputFile "Something went wrong!"
                    }
                    set description "Fragment $b + Fragment $c + Fragment $d + Fragment $e"
                }
            } else {
                puts "Something went wrong!"
                puts $outputFile "Something went wrong!"
            }

        }
        puts "\t\t$description\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

    #Export variable
    set dihedralTotalList $dihedEnergyTable

}

proc impropers {inputFile outputFile} {
    global colori colorf impTotalList
    puts "####################\n##$colori Out-of-plane $colorf\n####################"
    puts $outputFile "####################\n## Out-of-plane \n####################"
    # Read the inputFile. The following variables are set: numberAtoms pi types charges parameters connectivity xyz fragList
    source "$inputFile"

    # Get impropers list
    set uniqueImp {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set atoms1 [lindex $connectivity $index]

        if {[llength $atoms1] == 3} {
            set newLine [list [lindex $atoms1 0] [lindex $atoms1 1] $index [lindex $atoms1 2]]
                            
            if {[lsearch $uniqueImp $newLine] == -1} {
                lappend uniqueImp $newLine
            }
        }

    }

    ## Get ImpTrs
    set impParam {}
    foreach line $parameters {
        set var [regexp -inline {ImpTrs\s+"(\S+)"\s+"(\S+)"\s+"(\S+)"\s+"(\S+)"\s+(\S+)\s+(\S+)\s+(\S+)} $line]
        if {$var == ""} {
            continue
        } else {
            lappend impParam [lrange $var 1 7]
        }
    }

    # Calculate energy

    set impEnergyTable {}
    foreach imp $uniqueImp {
        set parameter ""
        set index0 [lindex $imp 0]
        set index1 [lindex $imp 1]
        set index2 [lindex $imp 2]
        set index3 [lindex $imp 3]

        set xa [lindex [lindex $xyz $index0] 0]
        set ya [lindex [lindex $xyz $index0] 1]
        set za [lindex [lindex $xyz $index0] 2]
        set xb [lindex [lindex $xyz $index1] 0]
        set yb [lindex [lindex $xyz $index1] 1]
        set zb [lindex [lindex $xyz $index1] 2]
        set xc [lindex [lindex $xyz $index2] 0]
        set yc [lindex [lindex $xyz $index2] 1]
        set zc [lindex [lindex $xyz $index2] 2]
        set xd [lindex [lindex $xyz $index3] 0]
        set yd [lindex [lindex $xyz $index3] 1]
        set zd [lindex [lindex $xyz $index3] 2]

        set type0 [lindex $types $index0]
        set type1 [lindex $types $index1]
        set type2 [lindex $types $index2]
        set type3 [lindex $types $index3]

        set xyz0 [list $xa $ya $za]
        set xyz1 [list $xb $yb $zb]
        set xyz2 [list $xc $yc $zc]
        set xyz3 [list $xd $yd $zd]

        set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $impParam $type0] $type1] $type2] $type3]
        set listMeasure [list $index0 $index1 $index2 $index3]
        if {$parameter == ""} {
            set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $impParam $type0] $type3] $type2] $type1]
            set listMeasure [list $index0 $index3 $index2 $index1]
        }
        if {$parameter == ""} {
            set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $impParam $type3] $type0] $type2] $type1]
            set listMeasure [list $index3 $index0 $index2 $index1]
        }
        if {$parameter == ""} {
            set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $impParam $type3] $type1] $type2] $type0]
            set listMeasure [list $index3 $index1 $index2 $index0]
        }
        if {$parameter == ""} {
            set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $impParam $type1] $type0] $type2] $type3]
            set listMeasure [list $index1 $index0 $index2 $index3]
        }
        if {$parameter == ""} {
            set parameter [lsearch -index 3 -all -inline -exact [lsearch -index 2 -all -inline -exact [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $impParam $type1] $type3] $type2] $type0]
            set listMeasure [list $index1 $index3 $index2 $index0]
        }


        if {$parameter != ""} {
            # Calculate the cross product for the first three atoms
            set u [::math::linearalgebra::sub $xyz0 $xyz2]
            set v [::math::linearalgebra::sub $xyz1 $xyz2]
            set p1 [::math::linearalgebra::crossproduct $u $v]
            set u [::math::linearalgebra::sub $xyz3 $xyz2]
            set v [::math::linearalgebra::sub $xyz2 $xyz1]
            set p2 [::math::linearalgebra::crossproduct $u $v]

            set omega [expr [::math::linearalgebra::angle $p1 $p2] / ($pi/180)]

            set mag [lindex [lindex $parameter 0] 4]
            set po [lindex [lindex $parameter 0] 5]
            set period [lindex [lindex $parameter 0] 6]

            set energy [expr ($mag * (1 - cos($period*(($omega-$po)*($pi/180))))) / 627.50947415151515]

            set logic0 [lsearch -all -regexp $fragList "(^| )${index0}($| )"]
            set logic1 [lsearch -all -regexp $fragList "(^| )${index1}($| )"]
            set logic2 [lsearch -all -regexp $fragList "(^| )${index2}($| )"]
            set logic3 [lsearch -all -regexp $fragList "(^| )${index3}($| )"]
            
            if {$logic0 == $logic1 && $logic1 == $logic2 && $logic2 == $logic3} {
                # Atoms belong to the same fragmennt
                if {$logic0 == ""} {
                    # All the atoms belong to the Fragment X
                    if {[lsearch -index 0 $impEnergyTable "none"] == -1} {
                        lappend impEnergyTable [list "none" 0]
                    }
                    set pos [lsearch -index 0 $impEnergyTable "none"]
                    set prevEnergy [lindex [lindex $impEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set impEnergyTable [lreplace $impEnergyTable $pos $pos [list "none" $finalEnergy]]
                } else {
                    # All atoms belong to a certain fragment
                    if {[lsearch -index 0 -exact $impEnergyTable "$logic0"] == -1} {
                        lappend impEnergyTable [list "$logic0" 0]
                    }
                    set pos [lsearch -index 0 $impEnergyTable "$logic0"]
                    set prevEnergy [lindex [lindex $impEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set impEnergyTable [lreplace $impEnergyTable $pos $pos [list "$logic0" $finalEnergy]]
                }
            } else {
                # Atoms belong to different fragments
                set listLogic [lsort -unique [list $logic0 $logic1 $logic2 $logic3]]
                set listLogicLength [llength $listLogic]
                if {$listLogicLength == 2} {
                    if {[lsearch -index 0 -exact $impEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]"] == -1} {
                        lappend impEnergyTable [list "[lindex $listLogic 0]-[lindex $listLogic 1]" 0]
                    }
                    set pos [lsearch -index 0 -exact $impEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]"]  
                    set prevEnergy [lindex [lindex $impEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set impEnergyTable [lreplace $impEnergyTable $pos $pos [list "[lindex $listLogic 0]-[lindex $listLogic 1]" $finalEnergy]]
                } elseif {$listLogicLength == 3} {
                    if {[lsearch -index 0 -exact $impEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]"] == -1} {
                        lappend impEnergyTable [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]" 0]
                    }
                    set pos [lsearch -index 0 -exact $impEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]"]  
                    set prevEnergy [lindex [lindex $impEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set impEnergyTable [lreplace $impEnergyTable $pos $pos [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]" $finalEnergy]]
                } elseif {$listLogicLength == 4} {
                    if {[lsearch -index 0 -exact $impEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]"] == -1} {
                        lappend impEnergyTable [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]" 0]
                    }
                    set pos [lsearch -index 0 -exact $impEnergyTable "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]"]  
                    set prevEnergy [lindex [lindex $impEnergyTable $pos] 1]
                    set finalEnergy [expr $prevEnergy + $energy]
                    set impEnergyTable [lreplace $impEnergyTable $pos $pos [list "[lindex $listLogic 0]-[lindex $listLogic 1]-[lindex $listLogic 2]-[lindex $listLogic 3]" $finalEnergy]]
                } else {
                    puts "Something went wrong."
                }
            }

        } else {
            puts "/!\\ WARNING /!\\ Missing parameter for improper $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
            puts $outputFile "/!\\ WARNING /!\\ Missing parameter for improper $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
        }

    }

    # Calculating the final energy for impropers
    puts "\tEnergy Decomposition:"
    puts $outputFile "\tEnergy Decomposition:"
    set fe 0
    set impEnergyTable [lsort -index 0 $impEnergyTable]
    foreach line $impEnergyTable {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t\t\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            set aLength [llength $a]

            if {$aLength == 1} {
                set description "Fragment [lindex $a 0]\t\t\t\t\t"
            } elseif {$aLength == 2} {
                if {[lindex $a 0] != "" && [lindex $a 1] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t\t\t\t"
                } else {
                    set b [lindex $a 0]
                    if {$b == ""} {
                        set b [lindex $a 1]
                    }
                    set description "Fragment $b + Fragment X\t\t\t\t"
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
                    set description "Fragment $b + Fragment $c + Fragment X\t\t"
                }
            } elseif {$aLength == 4} {
                if {[lindex $a 0] != "" && [lindex $a 1] != "" && [lindex $a 2] != "" && [lindex $a 3] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1] + Fragment [lindex $a 2] + Fragment [lindex $a 3]"
                } else {
                    set b [lindex $a 0]
                    set c [lindex $a 1]
                    set d [lindex $a 2]
                    set e [lindex $a 3]

                    if {$b == ""} {
                        set b "X"
                    } elseif {$c == ""} {
                        set c "X"
                    } elseif {$d == ""} {
                        set d "X"
                    } elseif {$e == ""} {
                        set e "X"
                    } else {
                        puts "Something went wrong!"
                        puts $outputFile "Something went wrong!"
                    }
                    set description "Fragment $b + Fragment $c + Fragment $d + Fragment $e"
                }
            } else {
                puts "Something went wrong!"
                puts $outputFile "Something went wrong!"
            }

        }
        puts "\t\t$description\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

    # Export variable
    set impTotalList $impEnergyTable
}

proc nonbond {inputFile outputFile} {
    global colori colorf vdwTotalList coloumbTotalList
    # Read the inputFile. The following variables are set: numberAtoms pi types charges parameters connectivity xyz fragList
    source "$inputFile"

    puts "Calculating non-bond interactions..."

    set numberThreads [numberCPUs]

    set atomIncrement [format %.0f [expr $numberAtoms / (($numberThreads*(2+2*$numberThreads))/4) + 1]]

    set last 0
    for {set core 0} { $last < $numberAtoms } { incr core } {
        set first $last
        set last [expr $first + ($atomIncrement * ($core + 1))]
        
        if {$last > $numberAtoms} {
            set last $numberAtoms
        }

        # Create a temporary file to store the results from the threads
        file tempfile [subst $core]tempFile

        set [subst $core]cpu [thread::create {
                proc nonbond {atomStart atomEnd inputFile tempFile} {
                    source "$inputFile"


                    ## Get VDW
                    set vdwParam {}
                    foreach line $parameters {
                        set var [regexp -inline {VDW\s+"(\S+)"\s+(\S+)\s+(\S+)} $line]
                        if {$var == ""} {
                            continue
                        } else {
                            lappend vdwParam [lrange $var 1 3]
                        }
                    }

                    set b 0
                    set c 0

                    set coloumbEnergyTable {}
                    set vdwEnergyTable {}
                    for {set i $atomStart} { $i < $atomEnd } { incr i } {
                        # Get atoms 1-2 bonds away
                        set scale083Atoms {}
                        set scale0Atoms {}

                        set typei [lindex $types $i]
                        set chargei [lindex $charges $i]
                        set atoms1 [lindex $connectivity $i]
                        set scale0Atoms [concat $scale0Atoms $atoms1]

                        foreach atom $atoms1 {
                            set atoms2 [lindex $connectivity $atom]
                            set scale0Atoms [concat $scale0Atoms $atoms2]

                            foreach atom1 $atoms2 {
                                if {$atom1 != $i} {
                                    set atoms3 [lindex $connectivity $atom1]
                            
                                    foreach atom2 $atoms3 {
                                        if {$atom2 != $i && $atom2 != $atom} {
                                            lappend scale083Atoms $atom2
                                        }
                                    }
                                }
                            }
                        }

                        set scale0Atoms [lsort -unique $scale0Atoms]
                        set scale083Atoms [lsort -unique $scale083Atoms]

                        
                        for {set j [expr $i + 1]} { $j < $numberAtoms } { incr j } {
                            lassign [lindex $xyz $i] xi yi zi
                            lassign [lindex $xyz $j] xj yj zj

                            set distance [expr {( ($xi - $xj)**2 + ($yi - $yj)**2 + ($zi - $zj)**2 )**(0.5)}]
                            
                            set typej [lindex $types $j]
                            set chargej [lindex $charges $j]

                            ## VDW
                            unset xi yi zi xj yj zj

                            set atomi [lsearch -index 0 -inline -exact $vdwParam $typei]
                            if {$atomi != ""} {
                                set ri [lindex $atomi 1]
                                set ei [lindex $atomi 2]
                            } else {
                                puts "Atom $typei not found"
                                puts $outputFile "Atom $typei not found"
                            }
                            set atomj [lsearch -index 0 -inline -exact $vdwParam $typej]
                            if {$atomi != ""} {
                                set rj [lindex $atomj 1]
                                set ej [lindex $atomj 2]
                            } else {
                                puts "Atom $typej not found"
                                puts $outputFile "Atom $typej not found"
                            }

                            set rij [expr {$ri + $rj}]
                            set eij [expr {($ei * $ej)**(0.5)}]

                            unset ri ei rj ej

                            set energy [expr {((($eij * ($rij)**12) / (($distance)**12)) - ((2 * $eij * ($rij)**6) / (($distance)**6))) / 627.50947415151515}]


                            ## Coulomb

                            set energyCoulomb [expr {(332.063712827427 * (($chargei * $chargej) / $distance)) / 627.50947415151515}]
                            
                            set logic0 ""
                            set logic1 ""
                            set index 0
                            foreach frag $fragList {
                                if {[lsearch -exact $frag $i] != -1} {
                                    set logic0 $index
                                }
                                if {[lsearch -exact $frag $j] != -1} {
                                    set logic1 $index
                                }
                                incr index
                            }

                            if {[lsearch $scale0Atoms $j] != -1} {
                                # Do nothing
                            } elseif {[lsearch $scale083Atoms $j] != -1} {
                                set energy [expr $energy * 0.5]
                                set energyCoulomb [expr $energyCoulomb / 1.2]
                                
                                if {$logic0 == $logic1} {
                                    # Atoms belonging to the same fragment
                                    if {$logic0 == ""} {
                                        # Both atoms belonging to the Fragment X
                                        if {[lsearch -index 0 $vdwEnergyTable "none"] == -1} {
                                            lappend vdwEnergyTable [list "none" 0]
                                        }
                                        set pos [lsearch -index 0 $vdwEnergyTable "none"]
                                        set vdwEnergyTable [lreplace $vdwEnergyTable $pos $pos [list "none" [expr [lindex [lindex $vdwEnergyTable $pos] 1] + $energy]]]

                                        if {[lsearch -index 0 $coloumbEnergyTable "none"] == -1} {
                                            lappend coloumbEnergyTable [list "none" 0]
                                        }
                                        set pos [lsearch -index 0 $coloumbEnergyTable "none"]
                                        set coloumbEnergyTable [lreplace $coloumbEnergyTable $pos $pos [list "none" [expr [lindex [lindex $coloumbEnergyTable $pos] 1] + $energyCoulomb]]]
                                    } else {
                                        # Atoms belonging to a certain fragment
                                        if {[lsearch -index 0 $vdwEnergyTable "$logic0-$logic1"] == -1} {
                                             lappend vdwEnergyTable [list "$logic0-$logic1" 0]
                                        }
                                        set pos [lsearch -index 0 $vdwEnergyTable "$logic0-$logic1"]
                                        set vdwEnergyTable [lreplace $vdwEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $vdwEnergyTable $pos] 1] + $energy]]]


                                        if {[lsearch -index 0 $coloumbEnergyTable "$logic0-$logic1"] == -1} {
                                             lappend coloumbEnergyTable [list "$logic0-$logic1" 0]
                                        }
                                        set pos [lsearch -index 0 $coloumbEnergyTable "$logic0-$logic1"]
                                        set coloumbEnergyTable [lreplace $coloumbEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $coloumbEnergyTable $pos] 1] + $energyCoulomb]]]
                                    }
                                } else {
                                    # Atoms belonging to different fragments
                                    if {[lsearch -index 0 -exact $vdwEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $vdwEnergyTable             "$logic1-$logic0"] == -1} {
                                        lappend vdwEnergyTable [list "$logic0-$logic1" 0]
                                    }
                                    set pos [lsearch -index 0 -exact $vdwEnergyTable "$logic0-$logic1"]
                                    if {$pos == -1} {
                                        set pos [lsearch -index 0 -exact $vdwEnergyTable "$logic1-$logic0"]
                                    }
                                    set vdwEnergyTable [lreplace $vdwEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $vdwEnergyTable $pos] 1] + $energy]]]


                                    if {[lsearch -index 0 -exact $coloumbEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $coloumbEnergyTable             "$logic1-$logic0"] == -1} {
                                        lappend coloumbEnergyTable [list "$logic0-$logic1" 0]
                                    }
                                    set pos [lsearch -index 0 -exact $coloumbEnergyTable "$logic0-$logic1"]
                                    if {$pos == -1} {
                                        set pos [lsearch -index 0 -exact $coloumbEnergyTable "$logic1-$logic0"]
                                    }
                                    set coloumbEnergyTable [lreplace $coloumbEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $coloumbEnergyTable $pos] 1] + $energyCoulomb]]]
                                }

                            } else {
                                
                                if {$logic0 == $logic1} {
                                    # Atoms belonging to the same fragment
                                    if {$logic0 == ""} {
                                        # Both atoms belonging to the Fragment X
                                        if {[lsearch -index 0 $vdwEnergyTable "none"] == -1} {
                                            lappend vdwEnergyTable [list "none" 0]
                                        }
                                        set pos [lsearch -index 0 $vdwEnergyTable "none"]
                                        set vdwEnergyTable [lreplace $vdwEnergyTable $pos $pos [list "none" [expr [lindex [lindex $vdwEnergyTable $pos] 1] + $energy]]]

                                        if {[lsearch -index 0 $coloumbEnergyTable "none"] == -1} {
                                            lappend coloumbEnergyTable [list "none" 0]
                                        }
                                        set pos [lsearch -index 0 $coloumbEnergyTable "none"]
                                        set coloumbEnergyTable [lreplace $coloumbEnergyTable $pos $pos [list "none" [expr [lindex [lindex $coloumbEnergyTable $pos] 1] + $energyCoulomb]]]
                                    } else {
                                        # Atoms belonging to a certain fragment
                                        if {[lsearch -index 0 $vdwEnergyTable "$logic0-$logic1"] == -1} {
                                             lappend vdwEnergyTable [list "$logic0-$logic1" 0]
                                        }
                                        set pos [lsearch -index 0 $vdwEnergyTable "$logic0-$logic1"]
                                        set vdwEnergyTable [lreplace $vdwEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $vdwEnergyTable $pos] 1] + $energy]]]


                                        if {[lsearch -index 0 $coloumbEnergyTable "$logic0-$logic1"] == -1} {
                                             lappend coloumbEnergyTable [list "$logic0-$logic1" 0]
                                        }
                                        set pos [lsearch -index 0 $coloumbEnergyTable "$logic0-$logic1"]
                                        set coloumbEnergyTable [lreplace $coloumbEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $coloumbEnergyTable $pos] 1] + $energyCoulomb]]]
                                    }
                                } else {
                                    # Atoms belonging to different fragments
                                    if {[lsearch -index 0 -exact $vdwEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $vdwEnergyTable             "$logic1-$logic0"] == -1} {
                                        lappend vdwEnergyTable [list "$logic0-$logic1" 0]
                                    }
                                    set pos [lsearch -index 0 -exact $vdwEnergyTable "$logic0-$logic1"]
                                    if {$pos == -1} {
                                        set pos [lsearch -index 0 -exact $vdwEnergyTable "$logic1-$logic0"]
                                    }
                                    set vdwEnergyTable [lreplace $vdwEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $vdwEnergyTable $pos] 1] + $energy]]]


                                    if {[lsearch -index 0 -exact $coloumbEnergyTable "$logic0-$logic1"] == -1 && [lsearch -index 0 -exact $coloumbEnergyTable             "$logic1-$logic0"] == -1} {
                                        lappend coloumbEnergyTable [list "$logic0-$logic1" 0]
                                    }
                                    set pos [lsearch -index 0 -exact $coloumbEnergyTable "$logic0-$logic1"]
                                    if {$pos == -1} {
                                        set pos [lsearch -index 0 -exact $coloumbEnergyTable "$logic1-$logic0"]
                                    }
                                    set coloumbEnergyTable [lreplace $coloumbEnergyTable $pos $pos [list "$logic0-$logic1" [expr [lindex [lindex $coloumbEnergyTable $pos] 1] + $energyCoulomb]]]
                                }

                            }
 
                        }

                    }

                    set file [open "$tempFile" w]
                    puts $file [list $coloumbEnergyTable $vdwEnergyTable]
                    close $file
                }

                thread::wait
            }]
        
        ::thread::send -async [set ${core}cpu] [list nonbond $first $last $inputFile [subst $${core}tempFile]] result

        set maxCore $core
    }


    for {set core 0} { $core <= $maxCore } { incr core } {
        vwait result
    }    

    set vdwList {}
    set coloumbList {}
    for {set core 0} { $core <= $maxCore } { incr core } {

        set file [open "[subst $${core}tempFile]" r]
        set var [read $file]

        lappend vdwList [lindex $var 1]
        lappend coloumbList [lindex $var 0]
    }
    set vdwList [lsort -index 0 [join $vdwList]]
    set coloumbList [lsort -index 0 [join $coloumbList]]
   
   set vdwFinal {}
   set totalVDW 0
   foreach line $vdwList {
        set a [lindex $line 0]

        set pos [lsearch -index 0 -exact $vdwFinal "$a"]
        if {$pos == -1} {
           lappend vdwFinal $line
        } else {
            set b [lindex [lsearch -index 0 -exact -inline $vdwFinal $a] 1]
            set c [expr $b + [lindex $line 1]]
            set vdwFinal [lreplace $vdwFinal $pos $pos [list "$a" "$c"]]
        }
   }


   set coloumbFinal {}
   set a ""
   set b ""
   set c ""
   foreach line $coloumbList {
        set a [lindex $line 0]

        set pos [lsearch -index 0 -exact $coloumbFinal $a]
        if {$pos == -1} {
            lappend coloumbFinal $line
        } else {
            set b [lindex [lsearch -index 0 -exact -inline $coloumbFinal $a] 1]
            set c [expr $b + [lindex $line 1]]
            set coloumbFinal [lreplace $coloumbFinal $pos $pos [list "$a" "$c"]]
        }
   }

   # Calculating the final energy for VWD interactions
    puts "####################\n## $colori VDW $colorf\n####################"
    puts $outputFile "####################\n## VDW \n####################"
    puts "\tEnergy Decomposition:"
    puts $outputFile "\tEnergy Decomposition:"
    set fe 0
    set vdwFinal [lsort -index 0 $vdwFinal]
    foreach line $vdwFinal {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            if {[lindex $a 0] == [lindex $a 1]} {
                set description "Fragment [lindex $a 0]\t\t"
            } elseif {[lindex $a 0] == ""} {
                set description "Fragment [lindex $a 1] + Fragment X\t"
            } elseif {[lindex $a 1] == ""} {
                set description "Fragment [lindex $a 0] + Fragment X\t"
            } else {
                set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t"
            }
        }
        puts "\t\t$description\t\t\t\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t\t\t\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

    # Export variable
    set vdwTotalList $vdwFinal

    # Calculating the final energy for Coloumb interactions
    puts "####################\n## $colori Coloumb $colorf\n####################"
    puts $outputFile "####################\n## Coloumb \n####################"
    puts "\tEnergy Decomposition:"
    puts $outputFile "\tEnergy Decomposition:"
    set fe 0
    set coloumbFinal [lsort -index 0 $coloumbFinal]
    foreach line $coloumbFinal {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            if {[lindex $a 0] == [lindex $a 1]} {
                set description "Fragment [lindex $a 0]\t\t"
            } elseif {[lindex $a 0] == ""} {
                set description "Fragment [lindex $a 1] + Fragment X\t"
            } elseif {[lindex $a 1] == ""} {
                set description "Fragment [lindex $a 0] + Fragment X\t"
            } else {
                set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t"
            }
        }
        puts "\t\t$description\t\t\t\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t\t\t\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

    # Export variable
    set coloumbTotalList $coloumbFinal

}

proc totalEnergies {outputFile} {
    global colori colorf bondTotalList angleTotalList dihedralTotalList impTotalList vdwTotalList coloumbTotalList

    set list [lsort -index 0 [join [list $bondTotalList $angleTotalList $dihedralTotalList $impTotalList $vdwTotalList $coloumbTotalList]]]

    set finalEnergies {}
    foreach line $list {
        set a [lindex $line 0]
        set z [lsort [split $a "-"]]

        if {[llength $z] == 2} {
            if {[lindex $z 0] == [lindex $z 1]} {
                set a "[lindex $z 0]"
            } else {
                set a "[lindex $z 0]-[lindex $z 1]"
            }
        }

        set pos [lsearch -index 0 -exact $finalEnergies "$a"]
        if {$pos == -1} {
            lappend finalEnergies [list $a [lindex $line 1]]
        } else {
            set b [lindex [lsearch -index 0 -exact -inline $finalEnergies "$a"] 1]
            set c [expr $b + [lindex $line 1]]
            set finalEnergies [lreplace $finalEnergies $pos $pos [list "$a" "$c"]]
        }
   }

    # Calculating the final energy for Total Energy interactions
    puts "####################\n## $colori Total Energy $colorf\n####################"
    puts $outputFile "####################\n## Total Energy \n####################"
    puts "\tEnergy Decomposition:"
    puts $outputFile "\tEnergy Decomposition:"
    set fe 0
    set finalEnergies [lsort -index 0 $finalEnergies]
    foreach line $finalEnergies {
        set fe [expr $fe + [lindex $line 1]]
        if {[lindex $line 0] == "none"} {
            set description "Fragment X\t\t\t\t\t"
        } else {
            set a [split [lindex $line 0] "-"]
            set aLength [llength $a]

            if {$aLength == 1} {
                set description "Fragment [lindex $a 0]\t\t\t\t\t"
            } elseif {$aLength == 2} {
                if {[lindex $a 0] != "" && [lindex $a 1] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1]\t\t\t\t"
                } else {
                    set b [lindex $a 0]
                    if {$b == ""} {
                        set b [lindex $a 1]
                    }
                    set description "Fragment $b + Fragment X\t\t\t\t"
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
                    set description "Fragment $b + Fragment $c + Fragment X\t\t"
                }
            } elseif {$aLength == 4} {
                if {[lindex $a 0] != "" && [lindex $a 1] != "" && [lindex $a 2] != "" && [lindex $a 3] != ""} {
                    set description "Fragment [lindex $a 0] + Fragment [lindex $a 1] + Fragment [lindex $a 2] + Fragment [lindex $a 3]"
                } else {
                    set b [lindex $a 0]
                    set c [lindex $a 1]
                    set d [lindex $a 2]
                    set e [lindex $a 3]

                    if {$b == ""} {
                        set b "X"
                    } elseif {$c == ""} {
                        set c "X"
                    } elseif {$d == ""} {
                        set d "X"
                    } elseif {$e == ""} {
                        set e "X"
                    } else {
                        puts "Something went wrong!"
                        puts $outputFile "Something went wrong!"
                    }
                    set description "Fragment $b + Fragment $c + Fragment $d + Fragment $e"
                }
            } else {
                puts "Something went wrong!"
                puts $outputFile "Something went wrong!"
            }

        }
        puts "\t\t$description\t[lindex $line 1]\tHartree"
        puts $outputFile "\t\t$description\t[lindex $line 1]\tHartree"
    }

    puts "$colori\tTotal Energy: $fe Hartree\n$colorf"
    puts $outputFile "\tTotal Energy: $fe Hartree\n"

}

##### Get Number of Processors
##### From https://stackoverflow.com/questions/29482303/how-to-find-the-number-of-cpus-in-tcl
proc numberCPUs {} {
    # Windows puts it in an environment variable
    global tcl_platform env
    if {$tcl_platform(platform) eq "windows"} {
        return $env(NUMBER_OF_PROCESSORS)
    }

    # Check for sysctl (OSX, BSD)
    set sysctl [auto_execok "sysctl"]
    if {[llength $sysctl]} {
        if {![catch {exec {*}$sysctl -n "hw.ncpu"} cores]} {
            return $cores
        }
    }

    # Assume Linux, which has /proc/cpuinfo, but be careful
    if {![catch {open "/proc/cpuinfo"} f]} {
        set cores [regexp -all -line {^processor\s} [read $f]]
        close $f
        if {$cores > 0} {
            return $cores
        }
    }

    # No idea what the actual number of cores is; exhausted all our options
    # Fall back to returning 1; there must be at least that because we're running on it!
    return 1
}

### Start procedure
set dir [file rootname [file normalize "$argv"]]
set outputFilePath "${dir}_OUTPUT.txt"
set outputFile [open "$outputFilePath" w]

puts "################################################################################\n### Energy Split\n################################################################################\n### Output\n################################################################################\n"
puts $outputFile "################################################################################\n### Energy Split\n################################################################################\n### Output\n################################################################################\n"
puts "### Energy Split detected [numberCPUs] cpu(s) available.\n"
puts $outputFile "### Energy Split detected [numberCPUs] cpu(s) available.\n"
puts "####################\n### $colori Fragments $colorf\n####################"
puts $outputFile "####################\n### Fragments \n####################"
source $argv
set i 0
foreach line $fragList {
    puts "$colori\tFragment $i $colorf: $line"
    puts $outputFile "\tFragment $i: $line"
    incr i
}
puts "\n"


bond $argv $outputFile
angleDihedral $argv $outputFile
impropers $argv $outputFile
nonbond $argv $outputFile
totalEnergies $outputFile


puts $outputFile "################################################################################\n### The calculation finished successfully on [clock format [clock seconds] -format %Y_%b_%d] at [clock format [clock seconds] -format %H:%M:%S]. \n################################################################################\n### Developed by Henrique S. Fernandes (henriquefer11@gmail.com)\n################################################################################"

close $outputFile

puts "################################################################################\n### The calculation finished successfully on [clock format [clock seconds] -format %Y_%b_%d] at [clock format [clock seconds] -format %H:%M:%S]. \n### Output file: $outputFilePath\n################################################################################\n### Developed by Henrique S. Fernandes (henriquefer11@gmail.com)\n################################################################################"