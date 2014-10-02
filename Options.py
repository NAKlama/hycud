#!/usr/bin/python3

# DO NOT CHANGE THE NEXT LINE!
options_ver = 'v3.2.2'

##################################################
#          Default values and options            #
#------------------------------------------------#
# Please copy this file to ~/.hycud to keep      #
# your custom options between updates.           #
##################################################

# Default path to hydropro executable
default_calcHydr            = "hydropro10-lnx.exe"

# Default path to the REMO directory
default_REMOdir             = "/nmr5/nare/programs/REMO"

# Default number of threads to use for REMO execution
#
# NOTE: This will have the progress indicator for REMO skip around if
#       increased. Progress indicator is only a rough approximation when
#       multithreaded REMO execution is used
default_threads             = 1

# Default directory to create temporary storage directories in
#
# NOTE: For a machine with sufficient RAM, putting this onto a RAM-disk
#       will speed up the splitting step significantly
default_temporaryStorage    = "."

# Default niceness
#
# NOTE: A high value decreases the priority of the programs HYCUD executes,
#       allowing them to run in the background without disrupting other tasks
#       Lowering this value will increase priority
default_niceness            = 20

# Default fragment size
#
# The input proteins will be split into equal fragments of this size, if no
# Fragmentation was specified on the command line
default_fragmentSize        = 14

# Default template file
#
# This is the default name of the HydroPro Datafile used as a template
# A simple filename will search for a file of that name in the directory
# HYCUD was executed in, while an absolute path will make a specific file
# the default template
default_templateFile        = "hydropro.dat"

# Since HydroPro runs multi threded, this isn't usually needed, and may
# get you in trouble with your administator, since it will make the load
# average of the machine increase a lot
#
# On machines with many cores, activating this will give you an option to
# run several instances of HydroPro in parallel
# allow_HydroPro_MultiTreaded = True
allow_HydroPro_MultiTreaded = False

# HYCUD calculates the values displayed in the seciont 'SUMMARY' weighted by
# molecular weight or non exchangable protons. By default these values are not
# displayed, and since this isn't often needed, the option to display these values
# is hidden by default
#
# By setting the option to True (like in the commented line below) an option will
# be shown that will then display these weighted values
# allow_option_weighted_averages = True
allow_option_weighted_averages = False

# By setting this option to True, these weighted averages will be shown by default
# The previous option will be ignored in that case.
# default_show_weighted_averages = True
default_show_weighted_averages = False

# HYCUD can work in daemon mode. This mode is normally deactivates, but can be
# activated by setting this option to true
allow_option_daemon_mode = False

# The following options are needed for daemon mode, they are ignored otherwise
daemon_db_user = "User"
daemon_db_pass = "Password"
daemon_db_db   = "DB_name"
