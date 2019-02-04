package provide guiEnergySplit  0.1

proc energySplit::topGui {} {
    #### Check if the window exists
	if {[winfo exists $::energySplit::topGui]} {wm deiconify $::energySplit::topGui ;return $::energySplit::topGui}
	toplevel $::energySplit::topGui

	#### Title of the windows
	wm title $energySplit::topGui "Energy Split v$energySplit::version " ;

    #### Screen Size
    set sWidth [expr [winfo vrootwidth  $::energySplit::topGui] -0]
	set sHeight [expr [winfo vrootheight $::energySplit::topGui] -100]

    #### Window Size and Position
	wm geometry $::energySplit::topGui 700x400+[expr $sWidth / 2 - 700 / 2]+[expr $sHeight / 2 - 400 / 2]
	$::energySplit::topGui configure -background {white}


    ####################################################################################################################
	####################################################################################################################
	####################################################################################################################


    #### Pack background frame
    grid columnconfigure $energySplit::topGui  0   -weight 1
    grid rowconfigure $energySplit::topGui     0   -weight 1
    
    set f0 $energySplit::topGui.frame0
    grid [ttk::frame $f0] -in $energySplit::topGui -sticky news
    grid columnconfigure $f0  0   -weight 1
    grid rowconfigure $f0     1   -weight 1

    #### Header
    ## Left
    grid [ttk::frame $f0.h1 \
        -height 100 \
        ] -in $f0 -row 0 -column 0 -sticky news

        grid [ttk::label $f0.h1.nameApp \
        -text "Energy Split" \
        ] -in $f0.h1 -row 0 -column 0 -sticky w -padx [list 20 5] -pady [list 10 10]


    


    ## Right
    grid [ttk::frame $f0.h2 \
        -height 100 \
        ] -in $f0 -row 0 -column 1 -sticky news

        grid [ttk::button $f0.h2.removeButton \
        -command {energySplit::deleteLine} \
        -text "Delete" \
        ] -in $f0.h2 -row 0 -column 0 -sticky ens -pady [list 5 0] -padx [list 0 0]

        set imgAddButton [image create photo -height 36 -width 36 -format gif -data {R0lGODdhJAAkAPEAAAAAAP91GyZFySZFySH5BAEAAAIALAAAAAAkACQAAAJ5lI+py40Bo3TUyItDVbnfbXkitI2m5oic14wr+8KPm5hxV9s5HfLHOfP1VENikHZCHYkfUGsiaDoXTel0R4pCr1iUlVsMg5dF2Q+o1S1l6HXGjbnhzsK0sRzvvvV5+NafBagk94cXQZF095T4lcJ4eMU4Rrc3aYlQAAA7
}]

        grid [ttk::button $f0.h2.addButton \
        -command {energySplit::addLine} \
        -text "Add" \
        ] -in $f0.h2 -row 0 -column 1 -sticky ens -pady [list 5 0] -padx [list 5 20]

    grid columnconfigure $f0.h1     0   -weight 2
    grid columnconfigure $f0.h2     0   -weight 1
    grid columnconfigure $f0.h2     1   -weight 1


    ### Main
    grid [ttk::frame $f0.l1 \
        ] -in $f0 -row 1 -column 0 -columnspan 2 -sticky news

    grid columnconfigure $f0.l1     0   -weight 1

    grid [tablelist::tablelist $f0.l1.table \
		-showeditcursor true \
		-columns {0 "Fragment" center 0 "Selection" center} \
		-stretch 1 \
		-background white \
		-selectmode extended \
        -yscrollcommand [list $f0.l1.yscb set] \
		-xscrollcommand [list $f0.l1.xscb set] \
		-state normal \
		-borderwidth 2 \
		-relief flat \
		] -in $f0.l1 -row 0 -column 0 -sticky news -padx [list 10 0]

    grid [ttk::scrollbar $f0.l1.yscb \
			-orient vertical \
			-command [list $f0.l1.table  yview]\
			] -in $f0.l1 -row 0 -column 1 -sticky ns -padx [list 0 10]

	grid [ttk::scrollbar $f0.l1.xscb \
			-orient horizontal \
			-command [list $f0.l1.table  xview]\
			] -in $f0.l1 -row 1 -column 0 -sticky we -padx [list 10 0]

    $f0.l1.table insert 0 [list "X" "all"]
    $f0.l1.table configcolumns 1 -editable true
    $f0.l1.table configcells 0,1 -editable false

    bind $energySplit::topGui.frame0.l1.table <<TablelistCellUpdated>> {energySplit::applyTableList}

    grid rowconfigure $f0.l1     0   -weight 1

    ### Footer
    grid [ttk::frame $f0.f1 \
        ] -in $f0 -row 2 -column 0 -columnspan 2 -sticky news
        
    grid columnconfigure $f0.f1     0   -weight 1

    grid [ttk::button $f0.f1.saveButton \
        -command {energySplit::save} \
        -text "Save" \
        ] -in $f0.f1 -row 0 -column 0 -sticky ens -pady [list 5 5] -padx [list 0 0]

    grid [ttk::button $f0.f1.saveRunButton \
        -command {energySplit::saveRun} \
        -text "Save & Run" \
        ] -in $f0.f1 -row 0 -column 1 -sticky ens -pady [list 5 5] -padx [list 5 20]

}

proc energySplit::deleteLine {} {
    set selection [$energySplit::topGui.frame0.l1.table curselection]
    $energySplit::topGui.frame0.l1.table delete $selection
    set size [$energySplit::topGui.frame0.l1.table size]
    for {set index 1} { $index < $size } { incr index } {
        $energySplit::topGui.frame0.l1.table configcells [subst $index],0 -text [subst $index]
    }
    set allText "all"
    set currentSelected [$energySplit::topGui.frame0.l1.table getcolumns 1]
    set i 0
    foreach text $currentSelected {
        if {$i != 0} {
            append allText " and not ($text)"
            #$energySplit::topGui.frame0.l1.table configcells [subst $i],0 -background "[lindex $energySplit::colorsbg [expr $i -1]]"
        }
        incr i
    }
    $energySplit::topGui.frame0.l1.table configcells 0,1 -text "$allText"
}

proc energySplit::addLine {} {
    $energySplit::topGui.frame0.l1.table insert end [list "[$energySplit::topGui.frame0.l1.table size]" ""]
    set allText "all"
    set currentSelected [$energySplit::topGui.frame0.l1.table getcolumns 1]
    set i 0
    foreach text $currentSelected {
        if {$i != 0} {
            append allText " and not ($text)"
            #$energySplit::topGui.frame0.l1.table configcells [subst $i],0 -background "[lindex $energySplit::colorsbg [expr $i -1]]"
        }
        incr i
    }
    $energySplit::topGui.frame0.l1.table configcells 0,1 -text "$allText"
}

proc energySplit::applyTableList {} {
    set allText "all"
    set currentSelected [$energySplit::topGui.frame0.l1.table getcolumns 1]
    set i 0
    foreach text $currentSelected {
        if {$i != 0} {
            append allText " and not ($text)"
            #$energySplit::topGui.frame0.l1.table configcells [subst $i],0 -background "[lindex $energySplit::colorsbg [expr $i -1]]"
        }
        incr i
    }
    $energySplit::topGui.frame0.l1.table configcells 0,1 -text "$allText"
}

proc energySplit::save {} {
    # Get fragments
    set fragment [lrange [$energySplit::topGui.frame0.l1.table getcolumns 1] 1 end]

    # Get the destination file
    set outputFile [tk_getSaveFile -defaultextension "tcl" -initialfile "input-energySplit-[clock format [clock seconds] -format %Y%b%d_%H%M%S].tcl" -title "Energy Split - Saving input file..."]

    # Save the input file
    energySplit::initialProcedure $fragment $outputFile
}

proc energySplit::saveRun {} {
    # Get fragments
    set fragment [lrange [$energySplit::topGui.frame0.l1.table getcolumns 1] 1 end]

    # Get the destination file
    set outputFile [tk_getSaveFile -defaultextension "tcl" -initialfile "input-energySplit-[clock format [clock seconds] -format %Y%b%d_%H%M%S].tcl" -title "Energy Split - Saving input file..."]
    
    # Save the input file
    energySplit::initialProcedure $fragment "$outputFile"

    # Run the calculation

    exec tclsh $::energySplitpath/energySplitCalculation.tcl "$outputFile"
 
    tk_messageBox -icon info -type ok -title "Energy Split" -message "The calculation finished successfully.\n\nOutput:\n[file rootname [file normalize ${outputFile}]]_OUTPUT.txt"
}
