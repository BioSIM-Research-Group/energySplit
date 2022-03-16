package provide energySplit 0.1

#### Initial procedures ####
namespace eval energySplit:: {
    namespace export energySplit

    # Packages
    package require mainEnergySplit                 0.1
    package require guiEnergySplit                  0.1

    
    package require tablelist


    # Variables
    # The installation path of this plugin is stored in the energySplitPath variable
    variable version    "0.2.0"
    variable colorsbg [list blue red gray orange yellow tan green white pink cyan purple]
    variable loadedParameters ""

    # GUI
    variable topGui         ".energySplitGui"

    ## sed
		if {[string first "Windows" $::tcl_platform(os)] != -1} {
			variable sed "$::energySplitpath/windowsDependencies/sed.exe"
            variable grep "$::energySplitpath/windowsDependencies/grep.exe"
		} else {
			variable sed "sed"
            variable grep "grep"
		}

}

proc energySplit::start {} {
    # Launch GUI
    energySplit::topGui

    # Print success message
    puts "Ennergy Split v$energySplit::version was loaded sucessfully."
}