package require Thread

set numberAtoms [lindex $argv 0]
global numberAtoms

proc nonBondThreads {} {
    global pi parameters numberAtoms types charges connectivity xyz fragList

    set numberThreads [numberCPUs]
    set atomIncrement [format %.0f [expr $numberAtoms / $numberThreads + 1]]


    for {set core 0} { [expr $core * $atomIncrement] < $numberAtoms } { incr core } {
        set first [expr $core * $atomIncrement]
        set last [expr $first + $atomIncrement]

        puts "$first $last"
        
        if {$last > $numberAtoms} {
            set last $numberAtoms
        }

        set [subst $core]cpu [thread::create {
                source "tmpFileFrag.tcl"
                global pi parameters numberAtoms types charges connectivity xyz fragList
                proc nonbond {atomStart atomEnd} {
                    global pi parameters numberAtoms types charges connectivity xyz fragList

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


                    set vdwEnergyFrag 0
                    set vdwEnergyInteraction 0
                    set vdwEnergyNotFrag 0
                    set coloumbEnergyFrag 0
                    set coloumbEnergyInteraction 0
                    set coloumbEnergyNotFrag 0

                    for {set i $atomStart} { $i < $atomEnd } { incr i } {
                        # Get atoms 1-2 bonds away
                        set scale083Atoms {}
                        set scale0Atoms {}

                        # set sel [atomselect top "index $i"]
                        set typei [lindex $types $i]
                        set chargei [lindex $charges $i]
                        set atoms1 [lindex $connectivity $i]
                        set scale0Atoms [concat $scale0Atoms $atoms1]
                        # $sel delete

                        foreach atom $atoms1 {
                            # set sel1 [atomselect top "index $atom"]
                            set atoms2 [lindex $connectivity $atom]
                            set scale0Atoms [concat $scale0Atoms $atoms2]
                            # $sel1 delete

                            foreach atom1 $atoms2 {
                                if {$atom1 != $i} {
                                    # set sel2 [atomselect top "index $atom1"]
                                    set atoms3 [lindex $connectivity $atom1]
                                    # $sel2 delete
                            
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
                            set xi [lindex [lindex $xyz $i] 0]
                            set yi [lindex [lindex $xyz $i] 1]
                            set zi [lindex [lindex $xyz $i] 2]

                            set xj [lindex [lindex $xyz $j] 0]
                            set yj [lindex [lindex $xyz $j] 1]
                            set zj [lindex [lindex $xyz $j] 2]

                            set distance [expr {( ($xi - $xj)**2 + ($yi - $yj)**2 + ($zi - $zj)**2 )**(0.5)}]
                            
                            set typej [lindex $types $j]
                            set chargej [lindex $charges $j]

                            ## VDW
                            set ri 0
                            set ei 0
                            set rj 0
                            set ej 0
                            set xi 0
                            set yi 0
                            set zi 0
                            set xj 0
                            set yj 0
                            set zj 0

                            set atomi [lsearch -index 0 -inline -exact $vdwParam $typei]
                            if {$atomi == ""} {
                                puts "Atom $typei not found"
                            } else {
                                set ri [lindex $atomi 1]
                                set ei [lindex $atomi 2]
                            }
                            set atomj [lsearch -index 0 -inline -exact $vdwParam $typej]
                            if {$atomi == ""} {
                                puts "Atom $typej not found"
                            } else {
                                set rj [lindex $atomj 1]
                                set ej [lindex $atomj 2]
                            }

                            set rij [expr {$ri + $rj}]
                            set eij [expr {($ei * $ej)**(0.5)}]

                            set energy [expr {(($eij * ($rij)**12) / (($distance)**12)) - ((2 * $eij * ($rij)**6) / (($distance)**6))}]


                            ## Coulomb

                            set energyCoulomb [expr {332.063712827427 * (($chargei * $chargej) / $distance)}]

                            set logic0 [lsearch $fragList $i]
                            set logic1 [lsearch $fragList $j]

                            if {[lsearch $scale0Atoms $j] != -1} {
                                # Do nothing
                            } elseif {[lsearch $scale083Atoms $j] != -1} {
                                if {($logic0 == -1) && ($logic1 == -1)} {
                                    set vdwEnergyNotFrag [expr {$vdwEnergyNotFrag + ($energy * 0.5)}]
                                    set coloumbEnergyNotFrag [expr {$coloumbEnergyNotFrag + ($energyCoulomb / 1.2)}]
                                } elseif {($logic0 != -1) && ($logic1 != -1)} {
                                    set vdwEnergyFrag [expr {$vdwEnergyFrag + ($energy * 0.5)}]
                                    set coloumbEnergyFrag [expr {$coloumbEnergyFrag + ($energyCoulomb / 1.2)}]
                                } elseif {($logic0 == -1) && ($logic1 != -1)} {
                                    set vdwEnergyInteraction [expr {$vdwEnergyInteraction + ($energy * 0.5)}]
                                    set coloumbEnergyInteraction [expr {$coloumbEnergyInteraction + ($energyCoulomb / 1.2)}]
                                } elseif {($logic0 != -1) && ($logic1 == -1)} {
                                    set vdwEnergyInteraction [expr {$vdwEnergyInteraction + ($energy * 0.5)}]
                                    set coloumbEnergyInteraction [expr {$coloumbEnergyInteraction + ($energyCoulomb / 1.2)}]
                                } else {
                                    puts "An error occurred."
                                }

                            } else {
                                
                                if {($logic0 == -1) && ($logic1 == -1)} {
                                    set vdwEnergyNotFrag [expr {$vdwEnergyNotFrag + $energy}]
                                    set coloumbEnergyNotFrag [expr {$coloumbEnergyNotFrag + $energyCoulomb}]
                                } elseif {($logic0 != -1) && ($logic1 != -1)} {
                                    set vdwEnergyFrag [expr {$vdwEnergyFrag + $energy}]
                                    set coloumbEnergyFrag [expr {$coloumbEnergyFrag + $energyCoulomb}]
                                } elseif {($logic0 == -1) && ($logic1 != -1)} {
                                    set vdwEnergyInteraction [expr {$vdwEnergyInteraction + $energy}]
                                    set coloumbEnergyInteraction [expr {$coloumbEnergyInteraction + $energyCoulomb}]
                                } elseif {($logic0 != -1) && ($logic1 == -1)} {
                                    set vdwEnergyInteraction [expr {$vdwEnergyInteraction + $energy}]
                                    set coloumbEnergyInteraction [expr {$coloumbEnergyInteraction + $energyCoulomb}]
                                } else {
                                    puts "An error occurred."
                                }

                            }
                        }
                    }

                    puts [list "ENERGY" [expr {$vdwEnergyFrag / 627.50947415151515}] [expr {$vdwEnergyNotFrag / 627.50947415151515}] [expr {$vdwEnergyInteraction / 627.50947415151515}] [expr {$coloumbEnergyFrag / 627.50947415151515}] [expr {$coloumbEnergyNotFrag / 627.50947415151515}] [expr {$coloumbEnergyInteraction / 627.50947415151515}]]

                    # puts "VDW Energy: [expr $vdwEnergy / 627.50947415151515] Hartree | $vdwEnergy kcal/mol"
                    # puts "Coulomb Energy: [expr $coloumbEnergy / 627.50947415151515] Hartree | $coloumbEnergy kcal/mol"
                }
                thread::wait
            }]
        
        ::thread::send -async [set ${core}cpu] [list nonbond $first $last] result

        set maxCore $core
    }

    for {set core 0} { $core <= $maxCore } { incr core } {
        vwait result
    }

    for {set core 0} { $core <= $maxCore } { incr core } {
        ::thread::release [set ${core}cpu]
    }
    
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

puts [nonBondThreads]