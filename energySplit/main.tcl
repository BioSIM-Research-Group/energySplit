package provide mainEnergySplit 0.1

proc energySplit::initialProcedure {listFragments fileName} {
    set nFragments [llength $listFragments]

    # Check if all the fragments are unique. Each atom should belong to a unique fragment. Two fragments cannot include thhe same atom.
    set error 0
    for {set i 0} {$i < $nFragments} {incr i} {
        for {set index [expr $i + 1]} {$index < $nFragments} {incr index} {
            set test [[atomselect top "([lindex $listFragments $i]) and ([lindex $listFragments $index])"] get index]
            if {$test == ""} {
                continue
            } else {
                set error 1
                puts "Error: Two fragments share one or more atoms. Each atom should only belong to an unique fragment.\nAn error was found between fragments: \"[lindex $listFragments $i]\" and \"[lindex $listFragments $index]\""
            }
        }
    }

    # If no error was reported on the previous test, the calculation can start.
    if {$error == 1} {
        puts "\nThe calculation abort because of the previous error(s)."
        return
    } else {
        puts "\nPreparing to start the energy calculation...\n"
    }

    # Convert VMD selections to index lists
    set fragList {}
    foreach frag $listFragments {
        set a [atomselect top "$frag"]
        lappend fragList [$a get index]
        $a delete

    }


    # Prepare a input file
    #set fileName "input-energySplit-[clock format [clock seconds] -format %Y%b%d_%H%M%S].tcl"
    set inputFile [open "$fileName" w]

    # Writting the header of the input file
    puts $inputFile "###\n### Energy Split v$energySplit::version\n###\n### Developed by: Henrique S. Fernandes (henriquefer11@gmail.com)\n###\n### Input file\n## Date: [clock format [clock seconds] -format %Y%b%d]\n## Hour: [clock format [clock seconds] -format %H:%M:%S]\n## Name: [molinfo top get name]"

    # Setting Variables
    puts $inputFile "\n## Variables\n"
    puts $inputFile "set pi [format %.30f [expr {acos(-1)}]]"
    set all [atomselect top "all"]
    puts $inputFile "set numberAtoms [$all num]"
    puts $inputFile "set types [list [$all get type]]"
    puts $inputFile "set charges [list [$all get charge]]"
    puts $inputFile "set xyz [list [$all get [list x y z]]]"
    puts $inputFile "set connectivity [list [energySplit::connectivityFromVMD [$all num]]]"
    if {$energySplit::loadedParameters != ""} {
        set parameters [split $energySplit::loadedParameters "\n"]
    } else {
        puts "Trying to get parameters from molUP..."
        catch {set parameters [split [.molUP.frame0.major.mol[molinfo top].tabs.tabInput.param get 1.0 end] "\n"]}
    }
    puts $inputFile "set parameters [list $parameters]"
    puts $inputFile "set fragList [list $fragList]"
    $all delete

    close $inputFile

    # Calculation
    #exec tclsh $::energySplitpath/energySplitCalculation.tcl "$fileName"
}

proc energySplit::connectivityFromVMD {numberAtoms} {
    set connect {}

    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set a [$sel getbonds]
        lappend connect [lindex $a 0]
        $sel delete
    }

    return $connect
}


### The following procedures were provided by molUP
proc energySplit::loadPrmtop {} {

    set fileTypes {
                {{AMBER prmtop file (.prmtop)}       {.prmtop}        }
    }
    
    set path [tk_getOpenFile -filetypes $fileTypes -defaultextension ".prmtop" -title "Choose a AMBER prmtop file..."]

    set parametersInfoFromPrmtop ""

    if {$path != ""} {
        append parametersInfoFromPrmtop "NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2\n"

        ### Get Atom Types
        set listAtomTypes [energySplit::prmtopGetDataFromFlag "AMBER_ATOM_TYPE" $path]
            ## Uppercase all atom types
            set newList {}
            set lowerCaseList {}
            set alternativeAtomTypes [list J K X Y Z 8 9 I V 5 6 7 F G H Q R S T U L W]
            set i 0
            foreach atomType $listAtomTypes {
                if {[string is lower [string range $atomType 0 0]] == 1} {
                    if {[lsearch $lowerCaseList $atomType] == -1} {
                        lappend lowerCaseList $atomType

                        #Change type on VMD lists
                        set selection [atomselect top "type \"$atomType\""]
                        set type "[string toupper [string range $atomType 0 0][lindex $alternativeAtomTypes $i]]"
                        $selection set type $type

                        set atomType "[string toupper [string range $atomType 0 0][lindex $alternativeAtomTypes $i]]"
                        
                        incr i
                    } else {
                        set index [lsearch $lowerCaseList $atomType]
                        set atomType "[string toupper [string range $atomType 0 0][lindex $alternativeAtomTypes $index]]"
                    }

                    lappend newList $atomType
                } else {
                    lappend newList $atomType
                }
            }
            set listAtomTypes $newList


        set listAtomTypesUnique [lsort -unique $listAtomTypes]


        ### Get VDW Information
        set lennardJonesA [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "LENNARD_JONES_ACOEF" $path]]
        set lennardJonesB [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "LENNARD_JONES_BCOEF" $path]]
        set icoList [energySplit::prmtopGetDataFromFlag "NONBONDED_PARM_INDEX" $path]
        set atomTypesList [energySplit::prmtopGetDataFromFlag "ATOM_TYPE_INDEX" $path]

        set nTypes [lindex [lsort -real -decreasing $atomTypesList] 0]
        foreach index $listAtomTypesUnique {
            set atomTypeIndex [lindex [lsearch $listAtomTypes $index] 0]

            if {$atomTypeIndex != -1} {
                set valueAtomTypeIndex [lindex $atomTypesList $atomTypeIndex]

                set ico [expr $nTypes*($valueAtomTypeIndex-1)+$valueAtomTypeIndex - 1]
                set coefIndex [lindex $icoList $ico]

                set coefA [lindex $lennardJonesA [expr $coefIndex - 1]]
                set coefB [lindex $lennardJonesB [expr $coefIndex - 1]]

                if {$coefA != 0 && $coefB != 0 && $coefA != "" && $coefB != ""} {
                    set r [expr ((2*$coefA/$coefB)**(1/6.000000000000000))/2]
                    set e [expr $coefB / (4* $coefA/$coefB)]
                } else {
                    set r 0.0000
                    set e 0.0000
                }
                append parametersInfoFromPrmtop "VDW \"$index\" [format %8.4f $r] [format %8.4f $e]\n"
            }
        }


        ### Get Bonds Information
        set bondsIncH [energySplit::prmtopGetDataFromFlag "BONDS_INC_HYDROGEN" $path]
        set bondsNotH [energySplit::prmtopGetDataFromFlag "BONDS_WITHOUT_HYDROGEN" $path]
        set bondsList [concat $bondsIncH $bondsNotH]
        set bondForceList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "BOND_FORCE_CONSTANT" $path]]
        set bondEquiList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "BOND_EQUIL_VALUE" $path]]

        set uniqueList {}
        set numberBonds [expr [llength $bondsList] / 3]
        for {set index 0} { $index < $numberBonds } { incr index } {
            set atom1Index [expr [lindex $bondsList [expr ($index*3+0)]] /3]
            set atom2Index [expr [lindex $bondsList [expr ($index*3+1)]] /3]
            set valueIndexes [expr [lindex $bondsList [expr ($index*3+2)]] - 1]

            set atom1 [lindex $listAtomTypes $atom1Index]
            set atom2 [lindex $listAtomTypes $atom2Index]
            set force [lindex $bondForceList $valueIndexes]
            set equil [lindex $bondEquiList $valueIndexes]

            if {[string match {*\**} $atom1] == 1} {
                set atom1Search "[string range $atom1 0 0]\[\*\]"
            } else {
                set atom1Search $atom1
            }

             if {[string match {*\**} $atom2] == 1} {
                set atom2Search "[string range $atom2 0 0]\[\*\]"
            } else {
                set atom2Search $atom2
            }

            if {[lsearch -exact $uniqueList "[string trim $atom1] [string trim $atom2]"] == -1} {
                lappend uniqueList "[string trim $atom1] [string trim $atom2]"
                lappend uniqueList "[string trim $atom2] [string trim $atom1]"
                append parametersInfoFromPrmtop "HrmStr1 \"$atom1\" \"$atom2\" [format %6.2f $force] [format %6.4f $equil]\n"
            } else {
                #Do nothing
            }
        }


        ### Get Angles Information
        set anglesIncH [energySplit::prmtopGetDataFromFlag "ANGLES_INC_HYDROGEN" $path]
        set anglesNotH [energySplit::prmtopGetDataFromFlag "ANGLES_WITHOUT_HYDROGEN" $path]
        set anglesList [concat $anglesIncH $anglesNotH]
        set angleForceList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "ANGLE_FORCE_CONSTANT" $path]]
        set angleEquiList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "ANGLE_EQUIL_VALUE" $path]]

        set uniqueList {}
        set numberAngles [expr [llength $anglesList] / 4]
        for {set index 0} { $index < $numberAngles } { incr index } {
            set atom1Index [expr [lindex $anglesList [expr ($index*4+0)]] /3]
            set atom2Index [expr [lindex $anglesList [expr ($index*4+1)]] /3]
            set atom3Index [expr [lindex $anglesList [expr ($index*4+2)]] /3]
            set valueIndexes [expr [lindex $anglesList [expr ($index*4+3)]] - 1]

            set atom1 [lindex $listAtomTypes $atom1Index]
            set atom2 [lindex $listAtomTypes $atom2Index]
            set atom3 [lindex $listAtomTypes $atom3Index]
            set force [lindex $angleForceList $valueIndexes]
            set equil [expr [lindex $angleEquiList $valueIndexes]*180/3.1415926535897931]


            if {[string match {*\**} $atom1] == 1} {
                set atom1Search "[string range $atom1 0 0]\[\*\]"
            } else {
                set atom1Search $atom1
            }

            if {[string match {*\**} $atom2] == 1} {
                set atom2Search "[string range $atom2 0 0]\[\*\]"
            } else {
                set atom2Search $atom2
            }

            if {[string match {*\**} $atom3] == 1} {
                set atom3Search "[string range $atom3 0 0]\[\*\]"
            } else {
                set atom3Search $atom3
            }

            if {[lsearch -exact $uniqueList "[string trim $atom1] [string trim $atom2] [string trim $atom3]"] == -1} {
                lappend uniqueList "[string trim $atom1] [string trim $atom2] [string trim $atom3]"
                lappend uniqueList "[string trim $atom3] [string trim $atom2] [string trim $atom1]"
                append parametersInfoFromPrmtop "HrmBnd1 \"$atom1\" \"$atom2\" \"$atom3\"  [format %4.2f $force] [format %7.4f $equil]\n"
            } else {
                #Do nothing
            }
        }
                ## Add the angle parameters for TIP3P water molecules
                append parametersInfoFromPrmtop "HrmBnd1 \"HW\" \"HW\" \"OW\"   0.00   0.0000\n"
                append parametersInfoFromPrmtop "HrmBnd1 \"HW\" \"OW\" \"HW\"   0.00   0.0000\n"





       ### Get Dihedral Angles Information
        set dihedIncH [energySplit::prmtopGetDataFromFlag "DIHEDRALS_INC_HYDROGEN" $path]
        set dihedNotH [energySplit::prmtopGetDataFromFlag "DIHEDRALS_WITHOUT_HYDROGEN" $path]
        set dihedList [concat $dihedIncH $dihedNotH]
        set dihedForceList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "DIHEDRAL_FORCE_CONSTANT" $path]]
        set dihedPeriodList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "DIHEDRAL_PERIODICITY" $path]]
        set dihedPhaseList [energySplit::prmtopTableExtractValues [energySplit::prmtopGetDataFromFlag "DIHEDRAL_PHASE" $path]]


        set uniqueList {}
        set uniqueListDihed {}
        set uniqueListDihed1 {}
        set dihedralAnglesList {}
        set numberDihed [expr [llength $dihedList] / 5]
        for {set index 0} { $index < $numberDihed } { incr index } {
            set atom1Index [expr [lindex $dihedList [expr ($index*5+0)]] /3]
            set atom2Index [expr [lindex $dihedList [expr ($index*5+1)]] /3]
            set atom3Index [expr [lindex $dihedList [expr ($index*5+2)]] /3]
            set atom4Index [expr [lindex $dihedList [expr ($index*5+3)]] /3]
            set valueIndexes [expr [lindex $dihedList [expr ($index*5+4)]] - 1]

            if {$atom4Index < 0} {
                ## Improper Dihedral Angle
                set atom1 [lindex $listAtomTypes $atom1Index]
                set atom2 [lindex $listAtomTypes $atom2Index]
                set atom3 [lindex $listAtomTypes [expr abs($atom3Index)]]
                set atom4 [lindex $listAtomTypes [expr abs($atom4Index)]]
                set force [lindex $dihedForceList $valueIndexes]
                set period [lindex $dihedPeriodList $valueIndexes]
                set phase [expr [lindex $dihedPhaseList $valueIndexes]*180/3.1415926535897931]

                 if {[lsearch -exact $uniqueList "$atom1 $atom2 $atom3 $atom4"] == -1} {
                     lappend uniqueList "$atom1 $atom3 $atom2 $atom4"
                     lappend uniqueList "$atom3 $atom1 $atom2 $atom4"
                     lappend uniqueList "$atom2 $atom1 $atom3 $atom4"
                     lappend uniqueList "$atom1 $atom2 $atom3 $atom4"
                     lappend uniqueList "$atom3 $atom2 $atom1 $atom4"
                     lappend uniqueList "$atom2 $atom3 $atom1 $atom4"

                     lappend uniqueList "$atom4 $atom1 $atom3 $atom2"
                     lappend uniqueList "$atom4 $atom3 $atom1 $atom2"
                     lappend uniqueList "$atom4 $atom2 $atom1 $atom3"
                     lappend uniqueList "$atom4 $atom1 $atom2 $atom3"
                     lappend uniqueList "$atom4 $atom3 $atom2 $atom1"
                     lappend uniqueList "$atom4 $atom2 $atom3 $atom1"
                     append parametersInfoFromPrmtop "ImpTrs \"$atom1\" \"$atom2\" \"$atom3\" \"$atom4\" [format %4.1f $force] [format %5.1f $phase] 2.0\n"
                 } else {
                     #Do nothing
                 }

                
            } elseif {$atom4Index > 0} {
                set atom1 [lindex $listAtomTypes $atom1Index]
                set atom2 [lindex $listAtomTypes $atom2Index]
                set atom3 [lindex $listAtomTypes [expr abs($atom3Index)]]
                set atom4 [lindex $listAtomTypes [expr abs($atom4Index)]]
                set force [lindex $dihedForceList $valueIndexes]
                set period [lindex $dihedPeriodList $valueIndexes]
                set phase [expr [lindex $dihedPhaseList $valueIndexes]*180/3.1415926535897931]

                if {[lsearch $dihedralAnglesList "$atom1 $atom2 $atom3 $atom4 $force $phase $period"] == -1} {
                    lappend dihedralAnglesList "$atom1 $atom2 $atom3 $atom4 $force $phase $period"
                } else {

                }

                if {[lsearch -exact $uniqueListDihed "$atom1 $atom2 $atom3 $atom4"] == -1} {
                     lappend uniqueListDihed "$atom1 $atom2 $atom3 $atom4"
                     lappend uniqueListDihed "$atom4 $atom3 $atom2 $atom1"

                     lappend uniqueListDihed1 "$atom1 $atom2 $atom3 $atom4"
                 } else {
                     #Do nothing
                 }

            }
        }


        foreach dihed $uniqueListDihed1 {
            set values [lsearch -all $dihedralAnglesList "[subst $dihed] *"]

            set phase1 0
            set phase2 0
            set phase3 0
            set phase4 0
            set force1 0
            set force2 0
            set force3 0
            set force4 0
            set atom1 X
            set atom2 X
            set atom3 X
            set atom4 X

            foreach value $values {

                if {[lindex [lindex $dihedralAnglesList $value] 6] == 1} {

                    set force1 [lindex [lindex $dihedralAnglesList $value] 4]
                    set phase1 [lindex [lindex $dihedralAnglesList $value] 5]

                } elseif {[lindex [lindex $dihedralAnglesList $value] 6] == 2} {

                    set force2 [lindex [lindex $dihedralAnglesList $value] 4]
                    set phase2 [lindex [lindex $dihedralAnglesList $value] 5]

                } elseif {[lindex [lindex $dihedralAnglesList $value] 6] == 3} {

                    set force3 [lindex [lindex $dihedralAnglesList $value] 4]
                    set phase3 [lindex [lindex $dihedralAnglesList $value] 5]

                } elseif {[lindex [lindex $dihedralAnglesList $value] 6] == 4} {

                    set force4 [lindex [lindex $dihedralAnglesList $value] 4]
                    set phase4 [lindex [lindex $dihedralAnglesList $value] 5]

                } else {     
                }
            }
                
            append parametersInfoFromPrmtop "AmbTrs \"[lindex $dihed 0]\" \"[lindex $dihed 1]\" \"[lindex $dihed 2]\" \"[lindex $dihed 3]\" [format %3.0f $phase1] [format %3.0f $phase2] [format %3.0f $phase3] [format %3.0f $phase4]   [format %5.3f $force1] [format %5.3f $force2] [format %5.3f $force3] [format %5.3f $force4] 1.0\n"
            set phase1 0
            set phase2 0
            set phase3 0
            set phase4 0
            set force1 0
            set force2 0
            set force3 0
            set force4 0
            set atom1 X
            set atom2 X
            set atom3 X
            set atom4 X
        }

    } else {}

    set energySplit::loadedParameters $parametersInfoFromPrmtop

    tk_messageBox -icon info -type ok -title "Energy Split" -message "The parameters were successfully loaded."
}

proc energySplit::prmtopGetDataFromFlag {flag path} {
    catch {exec $energySplit::grep -n "%FLAG" $path} listFlagLines

    set listFlagLines [split $listFlagLines "\n"]
    
    set pos [lsearch -index 1 -all $listFlagLines $flag]

    set firstLine [expr [lindex [split [lindex $listFlagLines $pos] ":" ] 0] + 2]
    set lastLine [expr [lindex [split [lindex $listFlagLines [expr $pos + 1]] ":" ] 0] -1]

    catch {exec $energySplit::sed -n "$firstLine,$lastLine p" $path} data

    return $data

}

proc energySplit::prmtopTableExtractValues {data} {
    set data [split $data "\n"]
    set finalList {}
    foreach line $data {
        for {set i 0} {$i < 5} {incr i} {
            set value [string trim [string range $line [expr 16 * $i] [expr 16 * $i + 15]]]
            if {$value != ""} {
                lappend finalList $value
            }
        }
    }
    return $finalList
}
