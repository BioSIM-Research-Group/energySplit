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
    set parameters [split [.molUP.frame0.major.mol[molinfo top].tabs.tabInput.param get 1.0 end] "\n"]
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
