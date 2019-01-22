package provide mainEnergySplit 0.1

proc energySplit::initialProcedure {listFragments} {
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

    
}