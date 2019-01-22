package provide energySplit 0.1

#### Initial procedures ####
namespace eval energySplit:: {
    namespace export energySplit

    # Packages
    package require mainEnergySplit                 0.1


    # Variables
    # The installation path of this plugin is stored in the energySplitPath variable
    variable version    "0.1"

}

proc energySplit::start {} {
    puts "Ennergy Split v$energySplit::version was loaded sucessfully."
}