
# Global variables
global pi numberAtoms parameters types charges connectivity xyz

set ::tcl_precision 17

# Collect Connectivity from VMD
proc connectivityFromVMD {} {
    global numberAtoms
    set connect {}

    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set a [$sel getbonds]
        lappend connect [lindex $a 0]
        $sel delete
    }

    return $connect
}

# Defining initial variables
proc variables {} {
    global pi numberAtoms parameters types charges connectivity xyz totalEnergy
    # set pi
    set pi [format %.30f [expr {acos(-1)}]]

    # total number of atoms
    set all [atomselect top "all"]
    set numberAtoms [$all num]
    set types [$all get type]
    set charges [$all get charge]
    set xyz [$all get [list x y z]]
    set connectivity [connectivityFromVMD]
    $all delete

    # get parameters list from molUP
    set parameters [split [.molUP.frame0.major.mol[molinfo top].tabs.tabInput.param get 1.0 end] "\n"]

    set totalEnergy 0
}


#########################################################################
################# DIHEDRAL ##############################################
#########################################################################

proc anglesDihedrals {} {
    global pi numberAtoms parameters types charges totalEnergy

    set uniqueAngles {}
    set uniqueDihedral {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set atoms1 [lindex [$sel getbonds] 0]
        $sel delete

        foreach atom $atoms1 {
            set sel1 [atomselect top "index $atom"]
            set atoms2 [lindex [$sel1 getbonds] 0]
            $sel1 delete

            foreach atom1 $atoms2 {
                if {$atom1 != $index} {
                    set sel2 [atomselect top "index $atom1"]
                    set atoms3 [lindex [$sel2 getbonds] 0]
                    $sel2 delete

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
            set teta [measure angle [list $index0 $index1 $index2]]
            set kangle [lindex [lindex $parameter 0] 3]
            set teta0 [lindex [lindex $parameter 0] 4]

            set energyAngle [expr {((($teta - $teta0)*($pi/180))**2) * $kangle}]
            set angleEnergy [expr {$angleEnergy + $energyAngle}]
        } else {
            puts "Missing parameter for angle $type0-$type1-$type2 (indexes $index0 $index1 $index2)"
        }

    }
    puts "Angles Energy: [expr $angleEnergy / 627.50947415151515] Hartree | $angleEnergy kcal/mol"

    set totalEnergy [expr $totalEnergy + ($angleEnergy / 627.50947415151515)]

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
    set dihedralEnergy 0
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
            set dihedralEnergy [expr {$dihedralEnergy + $energyDihed}]
            
        } else {
            puts "Missing parameter for dihedral $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
        }

    }
    puts "Dihedral Energy: [expr $dihedralEnergy / 627.50947415151515] Hartree | $dihedralEnergy kcal/mol"

    set totalEnergy [expr $totalEnergy + ($dihedralEnergy / 627.50947415151515)]
}


proc bonds {} {
    global pi numberAtoms parameters types charges totalEnergy

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

    set bondEnergy 0
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set atomType [$sel get type]
        set bonds [lindex [$sel getbonds] 0]
        foreach bond $bonds {

            if {$bond > $index} {
                set parameter ""
                set sel1 [atomselect top "index $bond"]
                set atomType1 [$sel1 get type]
                set l [measure bond [list $index $bond]]

                set parameter [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $bondParam $atomType] $atomType1]
                if {$parameter == ""} {
                    set parameter [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $bondParam $atomType1] $atomType]
                }

                if {$parameter != ""} {
                    set l0 [lindex [lindex $parameter 0] 3]
                    set kb [lindex [lindex $parameter 0] 2]

                    set energy [expr {(($l - $l0)**2) * $kb}]
                    set bondEnergy [expr {$bondEnergy + $energy}]
                } else {
                    puts "Missing parameter for bond $atomType-$atomType1 (indexes $index $bond)"
                }

                $sel1 delete
            }
        }
        $sel delete
    }

    puts "Bonds Energy: [expr $bondEnergy / 627.50947415151515] Hartree | $bondEnergy kcal/mol"

    set totalEnergy [expr $totalEnergy + ($bondEnergy / 627.50947415151515)]
}

proc impropers {} {
    global pi parameters numberAtoms types charges totalEnergy

    # Get impropers list
    set uniqueImp {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set atoms1 [lsort [lindex [$sel getbonds] 0]]
        $sel delete

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
    set impEnergy 0
    foreach imp $uniqueImp {

        set parameter ""
        set index0 [lindex $imp 0]
        set index1 [lindex $imp 1]
        set index2 [lindex $imp 2]
        set index3 [lindex $imp 3]

        set sel [atomselect top "index $index0"]
        set type0 [$sel get type]
        set xyz0 [$sel get {x y z}]
        $sel delete
        set sel [atomselect top "index $index1"]
        set type1 [$sel get type]
        set xyz1 [$sel get {x y z}]
        $sel delete
        set sel [atomselect top "index $index2"]
        set type2 [$sel get type]
        set xyz2 [$sel get {x y z}]
        $sel delete
        set sel [atomselect top "index $index3"]
        set type3 [$sel get type]
        set xyz3 [$sel get {x y z}]
        $sel delete


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
            set u [vecsub [lindex $xyz2 0] [lindex $xyz1 0]]
            set v [vecsub [lindex $xyz2 0] [lindex $xyz0 0]]
            set p1 [veccross $u $v]
            set u [vecsub [lindex $xyz2 0] [lindex $xyz3 0]]
            set v [vecsub [lindex $xyz2 0] [lindex $xyz1 0]]
            set p2 [veccross $u $v]

            set omega [expr (acos([vecdot $p1 $p2]/([veclength $p1] * [veclength $p2])))/($pi/180)]


            set mag [lindex [lindex $parameter 0] 4]
            set po [lindex [lindex $parameter 0] 5]
            set period [lindex [lindex $parameter 0] 6]

            set energyImp [expr $mag * (1 - cos($period*(($omega-$po)*($pi/180))))]
            set impEnergy [expr $impEnergy + $energyImp]


            # puts "$omega $listMeasure $type0 $type1 $type2 $type3 [lindex [lindex $parameter 0] 0] [lindex [lindex $parameter 0] 1] [lindex [lindex $parameter 0] 2] [lindex [lindex $parameter 0] 3] $mag $po $period $energyImp [expr $energyImp / 627.50947415151515]"
            
        } else {
            puts "/!\\ Missing parameter for improper $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
        }

    }
    puts "Impropers Energy: [expr $impEnergy / 627.50947415151515] Hartree | $impEnergy kcal/mol"

    set totalEnergy [expr $totalEnergy + ($impEnergy / 627.50947415151515)]
}


proc nonBondMultiCore {} {
    global pi parameters numberAtoms types charges xyz connectivity totalEnergy

    set tmpFile [open "tmpFile.tcl" w]
    puts $tmpFile "set parameters [list $parameters]"
    puts $tmpFile "set types [list $types]"
    puts $tmpFile "set charges [list $charges]"
    puts $tmpFile "set connectivity [list $connectivity]"
    puts $tmpFile "set xyz [list $xyz]"
    puts $tmpFile "set pi $pi"
    puts $tmpFile "set numberAtoms $numberAtoms"
    close $tmpFile

    set a [exec tclsh nonbond.tcl $numberAtoms]
    set a [split $a "\n"]

    file delete -force "tmpFile.tcl"
    
    set VDW 0
    set elec 0
    foreach b $a {
        if {[lsearch -index 0 $b "ENERGY"] == 0} {
            set VDW [expr $VDW + [lindex $b 1]]
            set elec [expr $elec + [lindex $b 2]]
        }
    }

    puts "VDW Energy: $VDW Hartree | [expr $VDW * 627.50947415151515] kcal/mol"
    puts "Coulomb Energy: $elec Hartree | [expr $elec * 627.50947415151515] kcal/mol"

    set totalEnergy [expr $totalEnergy + $VDW + $elec]
}


#### Time #####
proc duration { secs } {
 set timeatoms [ list ]
 if { [ catch {
    foreach div { 86400 3600 60 1 } \
            mod { 0 24 60 60 } \
           name { day hour min sec } {
       set n [ expr {$secs / $div} ]
       if { $mod > 0 } { set n [ expr {$n % $mod} ] }
       if { $n > 1 } {
          lappend timeatoms "$n ${name}s"
       } elseif { $n == 1 } {
         lappend timeatoms "$n ${name}"
       }
    }
 } err ] } {
    return -code error "duration: $err"
 }
 return [ join $timeatoms ]
}

proc tempo {time0 time1} {
   ## transformar o tempo num integer
   scan "$time0" "%dh %dm %ds   %s %s %s" h0 m0 s0 mt0 d0 y0
   scan "$time1" "%dh %dm %ds   %s %s %s" h1 m1 s1 mt1 d1 y1
   set time0 [clock scan "$h0:$m0:$s0 $mt0 $d0 $y0"]
   set time1 [clock scan "$h1:$m1:$s1 $mt1 $d1 $y1"]
   ## contas de diferença do tempo
   set timeD [expr abs ($time0-$time1)]
   set timeDiff "1 secs"
   if {$timeD!=0} {set timeDiff [duration $timeD]}
   return $timeDiff
}

# ############### Start
# set time0 [clock format [clock seconds] -format "%Hh %Mm %Ss   %d %b %y"]
# puts "## Preparing system..."
# variables
# puts "Done!\n"
# puts "## Calculating energy of bonds..."
# bonds
# puts "Done!\n"
# puts "## Calculating energy of angles and dihedral angles..."
# anglesDihedrals
# puts "Done!\n"
# puts "## Calculating energy of impropers..."
# impropers
# puts "Done!\n"
# puts "## Calculating energy of non-bond interactions (Coulomb and VDW)..."
# nonBondMultiCore
# puts "Done!\n"
# set time1 [clock format [clock seconds] -format "%Hh %Mm %Ss   %d %b %y"]
# puts "#####\n#####\n"
# puts "Total Energy: $totalEnergy hartree | [expr $totalEnergy * 627.50947415151515] kcal/mol\nCalculation time: [tempo $time0 $time1]"
# puts "\n#####\n#####"
# ################ End

proc bondsFrag {selection} {
    global pi numberAtoms parameters types charges totalEnergy

    #fragList is the variable containing a list of all atoms belonging to the fragment
    set fragList [[atomselect top "$selection"] get index]

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

    set bondEnergyFrag 0
    set bondEnergyInteraction 0
    set bondEnergyNotFrag 0
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set atomType [$sel get type]
        set bonds [lindex [$sel getbonds] 0]
        foreach bond $bonds {

            if {$bond > $index} {
                set parameter ""
                set sel1 [atomselect top "index $bond"]
                set atomType1 [$sel1 get type]
                set l [measure bond [list $index $bond]]

                set parameter [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $bondParam $atomType] $atomType1]
                if {$parameter == ""} {
                    set parameter [lsearch -index 1 -all -inline -exact [lsearch -index 0 -all -inline -exact $bondParam $atomType1] $atomType]
                }

                if {$parameter != ""} {
                    set l0 [lindex [lindex $parameter 0] 3]
                    set kb [lindex [lindex $parameter 0] 2]

                    set energy [expr {(($l - $l0)**2) * $kb}]

                    set logic0 [lsearch $fragList $bond]
                    set logic1 [lsearch $fragList $index]

                    if {($logic0 == -1) && ($logic1 == -1)} {
                        set bondEnergyNotFrag [expr {$bondEnergyNotFrag + $energy}]
                    } elseif {($logic0 != -1) && ($logic1 != -1)} {
                        set bondEnergyFrag [expr {$bondEnergyFrag + $energy}]
                    } elseif {($logic0 == -1) && ($logic1 != -1)} {
                        set bondEnergyInteraction [expr {$bondEnergyInteraction + $energy}]
                    } elseif {($logic0 != -1) && ($logic1 == -1)} {
                        set bondEnergyInteraction [expr {$bondEnergyInteraction + $energy}]
                    } else {
                        puts "An error occurred."
                    }

                } else {
                    puts "Missing parameter for bond $atomType-$atomType1 (indexes $index $bond)"
                }

                $sel1 delete
            }
        }
        $sel delete
    }

    puts "::Bonds Energy:\n - Fragment: [expr $bondEnergyFrag / 627.50947415151515] Hartree | $bondEnergyFrag kcal/mol\n - Boundary: [expr $bondEnergyInteraction / 627.50947415151515] Hartree | $bondEnergyInteraction kcal/mol\n - Non-Fragment: [expr $bondEnergyNotFrag / 627.50947415151515] Hartree | $bondEnergyNotFrag kcal/mol\n TOTAL: [expr ($bondEnergyNotFrag + $bondEnergyFrag + $bondEnergyInteraction) / 627.50947415151515] Hartree | [expr $bondEnergyNotFrag + $bondEnergyFrag + $bondEnergyInteraction] kcal/mol\n"

}

proc anglesDihedralsFrag {selection} {
    global pi numberAtoms parameters types charges totalEnergy

    #fragList is the variable containing a list of all atoms belonging to the fragment
    set fragList [[atomselect top "$selection"] get index]

    set uniqueAngles {}
    set uniqueDihedral {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set atoms1 [lindex [$sel getbonds] 0]
        $sel delete

        foreach atom $atoms1 {
            set sel1 [atomselect top "index $atom"]
            set atoms2 [lindex [$sel1 getbonds] 0]
            $sel1 delete

            foreach atom1 $atoms2 {
                if {$atom1 != $index} {
                    set sel2 [atomselect top "index $atom1"]
                    set atoms3 [lindex [$sel2 getbonds] 0]
                    $sel2 delete

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
    set angleEnergyFrag 0
    set angleEnergyInteraction 0
    set angleEnergyNotFrag 0
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
            set teta [measure angle [list $index0 $index1 $index2]]
            set kangle [lindex [lindex $parameter 0] 3]
            set teta0 [lindex [lindex $parameter 0] 4]

            set energyAngle [expr {((($teta - $teta0)*($pi/180))**2) * $kangle}]

            set logic0 [lsearch $fragList $index0]
            set logic1 [lsearch $fragList $index1]
            set logic2 [lsearch $fragList $index2]

            if {[llength [lsearch -all -inline [list $logic0 $logic1 $logic2] "-1"]] == 3} {
                set angleEnergyNotFrag [expr {$angleEnergyNotFrag + $energyAngle}]
            } elseif {[llength [lsearch -all -inline [list $logic0 $logic1 $logic2] "-1"]] == 0} {
                set angleEnergyFrag [expr {$angleEnergyFrag + $energyAngle}]
            } else {
                set angleEnergyInteraction [expr {$angleEnergyInteraction + $energyAngle}]
            }

            
        } else {
            puts "Missing parameter for angle $type0-$type1-$type2 (indexes $index0 $index1 $index2)"
        }

    }
    puts "::Angles Energy:\n - Fragment: [expr $angleEnergyFrag / 627.50947415151515] Hartree | $angleEnergyFrag kcal/mol\n - Boundary: [expr $angleEnergyInteraction / 627.50947415151515] Hartree | $angleEnergyInteraction kcal/mol\n - Non-Fragment: [expr $angleEnergyNotFrag / 627.50947415151515] Hartree | $angleEnergyNotFrag kcal/mol\n TOTAL: [expr ($angleEnergyNotFrag + $angleEnergyFrag + $angleEnergyInteraction) / 627.50947415151515] Hartree | [expr $angleEnergyNotFrag + $angleEnergyFrag + $angleEnergyInteraction] kcal/mol\n"


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
    puts "::Dihedral Energy:\n - Fragment: [expr $dihedralEnergyFrag / 627.50947415151515] Hartree | $dihedralEnergyFrag kcal/mol\n - Boundary: [expr $dihedralEnergyInteraction / 627.50947415151515] Hartree | $dihedralEnergyInteraction kcal/mol\n - Non-Fragment: [expr $dihedralEnergyNotFrag / 627.50947415151515] Hartree | $dihedralEnergyNotFrag kcal/mol\n TOTAL: [expr ($dihedralEnergyNotFrag + $dihedralEnergyFrag + $dihedralEnergyInteraction) / 627.50947415151515] Hartree | [expr $dihedralEnergyNotFrag + $dihedralEnergyFrag + $dihedralEnergyInteraction] kcal/mol\n"

}

proc impropersFrag {selection} {
    global pi parameters numberAtoms types charges totalEnergy

    #fragList is the variable containing a list of all atoms belonging to the fragment
    set fragList [[atomselect top "$selection"] get index]

    # Get impropers list
    set uniqueImp {}
    for {set index 0} { $index < $numberAtoms } { incr index } {
        set sel [atomselect top "index $index"]
        set atoms1 [lsort [lindex [$sel getbonds] 0]]
        $sel delete

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
    set impEnergyFrag 0
    set impEnergyInteraction 0
    set impEnergyNotFrag 0
    foreach imp $uniqueImp {

        set parameter ""
        set index0 [lindex $imp 0]
        set index1 [lindex $imp 1]
        set index2 [lindex $imp 2]
        set index3 [lindex $imp 3]

        set sel [atomselect top "index $index0"]
        set type0 [$sel get type]
        set xyz0 [$sel get {x y z}]
        $sel delete
        set sel [atomselect top "index $index1"]
        set type1 [$sel get type]
        set xyz1 [$sel get {x y z}]
        $sel delete
        set sel [atomselect top "index $index2"]
        set type2 [$sel get type]
        set xyz2 [$sel get {x y z}]
        $sel delete
        set sel [atomselect top "index $index3"]
        set type3 [$sel get type]
        set xyz3 [$sel get {x y z}]
        $sel delete


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
            set u [vecsub [lindex $xyz2 0] [lindex $xyz1 0]]
            set v [vecsub [lindex $xyz2 0] [lindex $xyz0 0]]
            set p1 [veccross $u $v]
            set u [vecsub [lindex $xyz2 0] [lindex $xyz3 0]]
            set v [vecsub [lindex $xyz2 0] [lindex $xyz1 0]]
            set p2 [veccross $u $v]

            set omega [expr (acos([vecdot $p1 $p2]/([veclength $p1] * [veclength $p2])))/($pi/180)]


            set mag [lindex [lindex $parameter 0] 4]
            set po [lindex [lindex $parameter 0] 5]
            set period [lindex [lindex $parameter 0] 6]

            set energyImp [expr $mag * (1 - cos($period*(($omega-$po)*($pi/180))))]

            set logic0 [lsearch $fragList $index0]
            set logic1 [lsearch $fragList $index1]
            set logic2 [lsearch $fragList $index2]
            set logic3 [lsearch $fragList $index3]

            if {[llength [lsearch -all -inline [list $logic0 $logic1 $logic2 ] "-1"]] == 4} {
                set impEnergyNotFrag [expr {$impEnergyNotFrag + $energyImp}]
            } elseif {[llength [lsearch -all -inline [list $logic0 $logic1 $logic2 $logic3] "-1"]] == 0} {
                set impEnergyFrag [expr {$impEnergyFrag + $energyImp}]
            } else {
                set impEnergyInteraction [expr {$impEnergyInteraction + $energyImp}]
            }
            
        } else {
            puts "/!\\ Missing parameter for improper $type0-$type1-$type2-$type3 (indexes $index0 $index1 $index2 $index3)"
        }

    }
    puts "::Impropers Energy:\n - Fragment: [expr $impEnergyFrag / 627.50947415151515] Hartree | $impEnergyFrag kcal/mol\n - Boundary: [expr $impEnergyInteraction / 627.50947415151515] Hartree | $impEnergyInteraction kcal/mol\n - Non-Fragment: [expr $impEnergyNotFrag / 627.50947415151515] Hartree | $impEnergyNotFrag kcal/mol\n TOTAL: [expr ($impEnergyNotFrag + $impEnergyFrag + $impEnergyInteraction) / 627.50947415151515] Hartree | [expr $impEnergyNotFrag + $impEnergyFrag + $impEnergyInteraction] kcal/mol\n"
}

proc nonBondMultiCoreFrag {selection} {
    global pi parameters numberAtoms types charges xyz connectivity totalEnergy

    set tmpFile [open "tmpFileFrag.tcl" w]
    puts $tmpFile "set parameters [list $parameters]"
    puts $tmpFile "set types [list $types]"
    puts $tmpFile "set charges [list $charges]"
    puts $tmpFile "set connectivity [list $connectivity]"
    puts $tmpFile "set xyz [list $xyz]"
    puts $tmpFile "set pi $pi"
    puts $tmpFile "set numberAtoms $numberAtoms"
    set list [[atomselect top "$selection"] get index]
    puts $tmpFile "set fragList [list $list]"
    close $tmpFile

    set a [exec tclsh nonbondFrag.tcl $numberAtoms]
    set a [split $a "\n"]

    file delete -force "tmpFileFrag.tcl"
    
    set vdwEnergyFrag 0
    set vdwEnergyInteraction 0
    set vdwEnergyNotFrag 0
    set coloumbEnergyFrag 0
    set coloumbEnergyInteraction 0
    set coloumbEnergyNotFrag 0
    foreach b $a {
        if {[lsearch -index 0 $b "ENERGY"] == 0} {
            set vdwEnergyFrag [expr $vdwEnergyFrag + [lindex $b 1]]
            set vdwEnergyInteraction [expr $vdwEnergyInteraction + [lindex $b 2]]
            set vdwEnergyNotFrag [expr $vdwEnergyNotFrag + [lindex $b 3]]
            set coloumbEnergyFrag [expr $coloumbEnergyFrag + [lindex $b 4]]
            set coloumbEnergyInteraction [expr $coloumbEnergyInteraction + [lindex $b 5]]
            set coloumbEnergyNotFrag [expr $coloumbEnergyNotFrag + [lindex $b 6]]
        }
    }

    puts "::VDW Energy:\n - Fragment: $vdwEnergyFrag Hartree | [expr $vdwEnergyFrag * 627.50947415151515] kcal/mol\n - Boundary: $vdwEnergyInteraction Hartree | [expr $vdwEnergyInteraction * 627.50947415151515] kcal/mol\n - Non-Fragment: $vdwEnergyNotFrag Hartree | [expr $vdwEnergyNotFrag * 627.50947415151515] kcal/mol\n TOTAL: [expr ($vdwEnergyNotFrag + $vdwEnergyFrag + $vdwEnergyInteraction)] Hartree | [expr ($vdwEnergyNotFrag + $vdwEnergyFrag + $vdwEnergyInteraction) * 627.50947415151515] kcal/mol\n"

    puts "::Coloumb Energy:\n - Fragment: $coloumbEnergyFrag Hartree | [expr $coloumbEnergyFrag * 627.50947415151515] kcal/mol\n - Boundary: $coloumbEnergyInteraction Hartree | [expr $coloumbEnergyInteraction * 627.50947415151515] kcal/mol\n - Non-Fragment: $coloumbEnergyNotFrag Hartree | [expr $coloumbEnergyNotFrag * 627.50947415151515] kcal/mol\n TOTAL: [expr ($coloumbEnergyNotFrag + $coloumbEnergyFrag + $coloumbEnergyInteraction)] Hartree | [expr ($coloumbEnergyNotFrag + $coloumbEnergyFrag + $coloumbEnergyInteraction) * 627.50947415151515] kcal/mol\n"


}

# variables
# bondsFrag "index 1 to 10"
# anglesDihedralsFrag "index 1 to 10"
# impropersFrag "index 1 to 10"
# nonBondMultiCoreFrag "index 1 to 10"


