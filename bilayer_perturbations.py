#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#=======================================================================
# create parser
#=======================================================================
version_nb="0.1.0"
parser = argparse.ArgumentParser(prog='bilayer_perturbations', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
********************************************************
v''' + version_nb + '''
author: Jean Helie
git: https://github.com/jhelie/bilayer_perturbations.git
DOI: 
********************************************************

[ Description ]

This script calculates the evolution of a lipid bilayer thickness and of the lipids tails
order parameters.

If proteins are present in the system the script can calculate the local perturbation of
those two metrics induced by transmembrane proteins.

The statistics and graphs on metrics and their perturbations are broken down by leaflet, 
lipid species and transmembrane protein clusters size.

The perturbations calculated can also be visualised in VMD by loading the text files
produced by the script (see https://github.com/jhelie/set_users_fields.git).

Visit the DOI for related paper and full information or see below for a summary of how the
script works.

bilayer thickness
-----------------

A thickness is associated to each lipids headgroups for each frame. It is calculated as 
the average geometric distance between the head group and its closest headgroup neighbours
in the opposite leaflet. This means that the shape of the bilayer does not prevent
meaningful calculations of its thickness.
The number of neighbours to take into account in the opposite leaflet can be set with the
--neighbours option.

lipids tails order parameter
----------------------------

This script computes the second rank order parameter as defined by:
 P2 = 0.5*(3*<cos**2(theta)> - 1)

where theta is the angle between the bond and the bilayer normal.
 P2 = 1      perfect alignement with the bilayer normal
 P2 = 0      random orientation
 P2 = -0.5   anti-alignement

The bilayer normal is considered to be z axis. This means that, unlike for the calculation
of the thickness, the more your bilayer deforms the further from the actual local normal
the z axis will be and the less meaningful the calculated order parameters will be.

detection of transmembrane protein clusters
-------------------------------------------

Two clustering algorithms can be used to identify protein clusters.
->Connectivity based (relies on networkX module):
  A protein is considered in a cluster if it is within a distance less than --nx_cutoff
  from another protein. This means that a single protein can act as a connector between
  two otherwise disconnected protein clusters.
  This algorithm can be ran using either the minimum distante between proteins (default, 
  --algorithm 'min') or the distance between their center of geometry (--algorithm 'cog').
  The 'min' option scales as the square of the number of proteins and can thus be very
  slow for large systems.

->Density based (relies on the sklearn module and its implementation of DBSCAN):
  A protein is considered in a cluster if is surrounded by at least --db_neighbours other
  proteins within a radius of --db_radius.
  This density based approach is usually less suited to the detection of protein
  clusters but as a general rule the more compact the clusters, the smaller --db_radius
  the higher --db_neighbours can be - for details on this algorithm see its online
  documentation.
  This algorithm is selected by setting the --algorithm option to 'density'.

The identified protein clusters are considered to be transmembrane only if the closest
lipid headgroup neighbours to the cluster particles are all within the same leaflet.

[ Requirements ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - networkX (if option --radial is used and option --algorithm set to 'min' or 'cog')
 - sklearn (if option --radial is used and option --algorithm set to 'density')

[ Notes ]

1. It's a good idea to pre-process the trajectory first and to only output the relevant
   particles (e.g. no water and no cholesterol).

2. By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
   the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
   routine.
   This optimisation process can take time in large systems and you can specify your own
   cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
    --leaflet_cutoff 15
   
   In very large systems (more then ~50,000 phospholipids) LeafletFinder can fail (or
   rather networkX, that it relies on) can fail , to avoid this
   you can choose not to use this routine:
    --leaflet_cutoff large
   In this case lipids whose headgroups z value is above the average lipids z value will be
   considered to make up the upper leaflet and those whos headgroups z value is below the average
   will be considered to make the lower leaflet. This means that your bilayer should be as flat as
   possible in the gro file that you supply in order to get a meaningful outcome.

3. In case lipids flipflop during the trajectory, a file listing them can be supplied via the -i flag.
   This file can be the output of the ff_detect script and should follow the format:
   'resname,resid,starting_leaflet' format on each line e.g. 'POPC,145,lower'
   If flipflopping lipids are not identified they may add some significant noise to the results

4. The code can easily be updated to add more lipids, for now the following tails can be dealt with:
    - Martini: DHPC,DHPE,DLPC,DLPE,DAPC,DUPC,DPPC,DPPE,DPPS,DPPG,DSPC,DSPE,POPC,POPE,POPS,POPG,PPCS,PIP2,PIP3,GM3

5. The order parameter calculated for each (handled) lipd specie can be visualised with VMD.
   This can be done either with pdb files (output frequency controled via -w flag) or with the 
   xtc trajectory.
     - pdb file: the order parameter info is stored in the beta factor column. Just open
                 the pdb with VMD and choose Draw Style > Coloring Method > Beta 
     - xtc file: the order parameter info is stored in a .txt file in /3_VMD/ and you can load it into
                 the user field in the xtc by sourcing the script 'set_user_fields.tcl' and running the
                 procedure 'set_order_param'

6. The colour associated to each lipid specie can be defined by supplying a colour file containing
   'resname,colour' on each line (a line with a colour MUST be defined for all species).
   Colours can be specified using single letter code (e.g. 'r'), hex code  or the name of colormap.
   In case a colormap is used, its name must be specified as the colour for each lipid specie - see
   the matplotlib website for a list of standard colour maps.
   If no colour is used the 'jet' colour map is used by default.

[ Usage ]
	
Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro]
-x			: trajectory file [.xtc] (optional)
-c			: colour definition file, see note 6
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		10	: process every t-frames
-w			: write annotated pdbs every [w] processed frames (optional, see note 5)
--radial		: calculate radial perturbations around protein clusters
--smooth		: nb of points to use for data smoothing (optional)
--perturb 	0	: perturbation to calculate (see note XX)
	       [1]	  thickness
		2	  order parameters
		3	  thickness + order parameters

Lipids identification  
-----------------------------------------------------
--flipflops		: input file with flipflopping lipids, see note 3
--forcefield		: forcefield options, see note 3
--leaflet_cutoff	: leaflet identification technique, see note 2

radial perturbations and protein clusters identification
-----------------------------------------------------
-g			: cluster groups definition file, see note XXX
-p			: protein selection file, (optional, see note XXX)
--radial_radius 	100	: max radius to which represent the radial perturbations (Angtrom)
--radial_bins 	25	: number of bins into which discretise data for radial perturbations
--algorithm	min	: 'cog','min' or 'density', see note XXX
--nx_cutoff 	8	: networkX cutoff distance for lipid-lipid contact (Angtrom)
--db_radius 	20	: DBSCAN search radius (Angtrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --radius	
 
Other options
-----------------------------------------------------
--neighbours	5	: nb of nearest opposite neighbours to use for thickness calculation
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-c', nargs=1, dest='colour_file', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[10000000000000], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('-w', nargs=1, dest='frames_write_dt', default=[1000000000000000], type=int, help=argparse.SUPPRESS)
parser.add_argument('--radial', dest='radial', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--smooth', nargs=1, dest='nb_smoothing', default=[0], type=int, help=argparse.SUPPRESS)
parser.add_argument('--perturb', dest='perturb', choices=['0','1','2','3'], default='1', help=argparse.SUPPRESS)

#lipids identification
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--forcefield', dest='forcefield_opt', choices=['martini'], default='martini', help=argparse.SUPPRESS)
parser.add_argument('--leaflet_cutoff', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#radial and protein clusters options
parser.add_argument('-g', nargs=1, dest='cluster_groups_file', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-p', nargs=1, dest='selection_file_prot', default=['auto'], help=argparse.SUPPRESS)
parser.add_argument('--radial_radius', nargs=1, dest='radial_radius', default=[100], type=float, help=argparse.SUPPRESS)
parser.add_argument('--radial_bins', nargs=1, dest='radial_nb_bins', default=[25], type=int, help=argparse.SUPPRESS)
parser.add_argument('--algorithm', dest='m_algorithm', choices=['cog','min','density'], default='min', help=argparse.SUPPRESS)
parser.add_argument('--nx_cutoff', nargs=1, dest='cutoff_connect', default=[8], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_radius', nargs=1, dest='dbscan_dist', default=[20], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_neighbours', nargs=1, dest='dbscan_nb', default=[3], type=int, help=argparse.SUPPRESS)

#other options
parser.add_argument('--neighbours', nargs=1, dest='thick_nb_neighbours', default=[5], type=int, help=argparse.SUPPRESS)
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=======================================================================
# store inputs
#=======================================================================
args=parser.parse_args()
#data options
args.grofilename=args.grofilename[0]
args.xtcfilename=args.xtcfilename[0]
args.colour_file=args.colour_file[0]
args.output_folder=args.output_folder[0]
args.t_start=args.t_start[0]
args.t_end=args.t_end[0]
args.frames_dt=args.frames_dt[0]
args.frames_write_dt=args.frames_write_dt[0]
args.nb_smoothing=args.nb_smoothing[0]
args.perturb=int(args.perturb[0])
#lipids identification
args.selection_file_ff=args.selection_file_ff[0]
#radial and protein clusters options
args.cluster_groups_file=args.cluster_groups_file[0]
args.selection_file_prot=args.selection_file_prot[0]
args.radial_radius=args.radial_radius[0]
args.radial_nb_bins=args.radial_nb_bins[0]
args.cutoff_connect=args.cutoff_connect[0]
args.dbscan_dist=args.dbscan_dist[0]
args.dbscan_nb=args.dbscan_nb[0]
#other options
args.cutoff_leaflet=args.cutoff_leaflet[0]
args.thick_nb_neighbours=args.thick_nb_neighbours[0]

#process options
if args.perturb == 0:
	args.radial=True
if args.radial:
	global radial_sizes
	global radial_step
	radial_sizes = {}
	radial_sizes["all frames"] = []
	radial_step = args.radial_radius/float(args.radial_nb_bins)
global lipids_ff_nb
global nb_frames_processed
lipids_ff_nb = 0
nb_frames_processed = 0
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflet_cutoff option should be a number or 'large', see note 2"
		sys.exit(1)

#=======================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=======================================================================

#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.colour_file!="no" and not os.path.isfile(args.colour_file):
	print "Error: file " + str(args.colour_file) + " not found."
	sys.exit(1)
if args.selection_file_ff!="no" and not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.xtcfilename=="no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
	elif '--smooth' in sys.argv:
		print "Error: --smooth option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=======================================================================
# create folders and log file
#=======================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "bilayer_perturbations_" + args.grofilename[:-4]
	else:
		args.output_folder = "bilayer_perturbations_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	
	#case: thickness to be calculated
	#--------------------------------
	if args.perturb == 1 or args.perturb == 3:
		os.mkdir(args.output_folder + "/thickness")
		#1 species
		os.mkdir(args.output_folder + "/thickness/1_species")
		if args.xtcfilename!="no":
			os.mkdir(args.output_folder + "/thickness/1_species/xvg")
			os.mkdir(args.output_folder + "/thickness/1_species/png")
			if args.nb_smoothing>1:
				os.mkdir(args.output_folder + "/thickness/1_species/smoothed")
				os.mkdir(args.output_folder + "/thickness/1_species/smoothed/png")
				os.mkdir(args.output_folder + "/thickness/1_species/smoothed/xvg")
		#2 snapshots
		os.mkdir(args.output_folder + "/thickness/2_snapshots")
		#3 vmd
		if args.xtcfilename!="no":
			os.mkdir(args.output_folder + "/thickness/3_VMD")
	
	#case: order parameters to be calculated
	#---------------------------------------
	if args.perturb == 2 or args.perturb == 3:
		os.mkdir(args.output_folder + "/order_param")
		#1 non flipfloppping lipids
		os.mkdir(args.output_folder + "/order_param/1_nff")
		if args.xtcfilename!="no":
			os.mkdir(args.output_folder + "/order_param/1_nff/xvg")
			os.mkdir(args.output_folder + "/order_param/1_nff/png")
			if args.nb_smoothing>1:
				os.mkdir(args.output_folder + "/order_param/1_nff/smoothed")
				os.mkdir(args.output_folder + "/order_param/1_nff/smoothed/png")
				os.mkdir(args.output_folder + "/order_param/1_nff/smoothed/xvg")
		#2 snapshots
		os.mkdir(args.output_folder + "/order_param/2_snapshots")
		#3 vmd
		if args.xtcfilename!="no":
			os.mkdir(args.output_folder + "/order_param/3_VMD")
		#4 flipflopping lipids
		if args.selection_file_ff!="no":
			os.mkdir(args.output_folder + "/order_param/4_ff")
			if args.xtcfilename!="no":
				os.mkdir(args.output_folder + "/order_param/4_ff/xvg")
				os.mkdir(args.output_folder + "/order_param/4_ff/png")
				if args.nb_smoothing>1:
					os.mkdir(args.output_folder + "/order_param/4_ff/smoothed")
					os.mkdir(args.output_folder + "/order_param/4_ff/smoothed/png")
					os.mkdir(args.output_folder + "/order_param/4_ff/smoothed/xvg")
	
	#case: radial
	#------------
	if args.radial:
		os.mkdir(args.output_folder + "/radial")
		#density
		os.mkdir(args.output_folder + "/radial/density")
		os.mkdir(args.output_folder + "/radial/density/sizes")
		os.mkdir(args.output_folder + "/radial/density/sizes/xvg")
		os.mkdir(args.output_folder + "/radial/density/sizes/png")
		os.mkdir(args.output_folder + "/radial/density/species")
		os.mkdir(args.output_folder + "/radial/density/species/xvg")
		os.mkdir(args.output_folder + "/radial/density/species/png")
		if args.xtcfilename != "no":
			os.mkdir(args.output_folder + "/radial/density/snapshots")
			os.mkdir(args.output_folder + "/radial/density/snapshots/sizes")
			os.mkdir(args.output_folder + "/radial/density/snapshots/sizes/xvg")
			os.mkdir(args.output_folder + "/radial/density/snapshots/sizes/png")
			os.mkdir(args.output_folder + "/radial/density/snapshots/species")
			os.mkdir(args.output_folder + "/radial/density/snapshots/species/xvg")
			os.mkdir(args.output_folder + "/radial/density/snapshots/species/png")

		#thickness
		if args.perturb == 1 or args.perturb == 3:
			os.mkdir(args.output_folder + "/radial/thickness")
			os.mkdir(args.output_folder + "/radial/thickness/sizes")
			os.mkdir(args.output_folder + "/radial/thickness/sizes/xvg")
			os.mkdir(args.output_folder + "/radial/thickness/sizes/png")
			os.mkdir(args.output_folder + "/radial/thickness/species")
			os.mkdir(args.output_folder + "/radial/thickness/species/xvg")
			os.mkdir(args.output_folder + "/radial/thickness/species/png")
			if args.xtcfilename != "no":
				os.mkdir(args.output_folder + "/radial/thickness/snapshots")			
				os.mkdir(args.output_folder + "/radial/thickness/snapshots/sizes")
				os.mkdir(args.output_folder + "/radial/thickness/snapshots/sizes/xvg")
				os.mkdir(args.output_folder + "/radial/thickness/snapshots/sizes/png")
				os.mkdir(args.output_folder + "/radial/thickness/snapshots/species")
				os.mkdir(args.output_folder + "/radial/thickness/snapshots/species/xvg")
				os.mkdir(args.output_folder + "/radial/thickness/snapshots/species/png")

		#order parameters
		if args.perturb == 2 or args.perturb == 3:
			os.mkdir(args.output_folder + "/radial/order_param")
			os.mkdir(args.output_folder + "/radial/order_param/sizes")
			os.mkdir(args.output_folder + "/radial/order_param/sizes/xvg")
			os.mkdir(args.output_folder + "/radial/order_param/sizes/png")
			os.mkdir(args.output_folder + "/radial/order_param/species")
			os.mkdir(args.output_folder + "/radial/order_param/species/xvg")
			os.mkdir(args.output_folder + "/radial/order_param/species/png")
			if args.xtcfilename != "no":
				os.mkdir(args.output_folder + "/radial/order_param/snapshots")
				os.mkdir(args.output_folder + "/radial/order_param/snapshots/sizes")
				os.mkdir(args.output_folder + "/radial/order_param/snapshots/sizes/xvg")
				os.mkdir(args.output_folder + "/radial/order_param/snapshots/sizes/png")
				os.mkdir(args.output_folder + "/radial/order_param/snapshots/species")
				os.mkdir(args.output_folder + "/radial/order_param/snapshots/species/xvg")
				os.mkdir(args.output_folder + "/radial/order_param/snapshots/species/png")
		
	#create log
	#----------
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/bilayer_perturbations.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[bilayer_perturbations v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python bilayer_perturbations.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	#copy input files
	#----------------
	if args.colour_file!="no":
		shutil.copy2(args.colour_file,args.output_folder + "/")
	if args.selection_file_ff!="no":
		shutil.copy2(args.selection_file_ff,args.output_folder + "/")

################################################################################################################################################
# DATA LOADING
################################################################################################################################################

def set_forcefield_lipids_info():
	#leaflet selection
	global leaflet_beads
	global leaflet_sele_string
	leaflet_beads = {}
	leaflet_beads['martini'] = "name PO4 or name PO3 or name B1A"
	leaflet_sele_string = leaflet_beads[args.forcefield_opt]

	#order parameter: lipids handled
	global op_lipids_possible
	op_lipids_possible = {}
	op_lipids_possible["martini"]=['DHPC','DHPE','DLPC','DLPE','DAPC','DUPC','DPPC','DPPE','DPPS','DPPG','DSPC','DSPE','POPC','POPE','POPS','POPG','PPCS','PIP2','PIP3','GM3']
	
	#order parameter: lipid tails
	global bond_names
	global tail_boundaries
	bond_names = {}
	tail_boundaries = {}
	if args.forcefield_opt=="martini":
		#define possible head, link and tail sequences
		head_PC=" NC3-PO4"
		head_PE=" NH3-PO4"
		head_PS=" CNO-PO4"
		head_PG=" GL0-PO4"
		head_I2=" PO1-RP1 PO2-RP1 PO2-RP2 RP1-RP2 RP1-RP3 RP2-RP3 RP3-PO3"
		head_I3=" PO0-RP2 PO1-RP1 PO2-RP1 PO2-RP2 RP1-RP2 RP1-RP3 RP2-RP3 RP3-PO3"
		head_GM=" very-big"
		link1=" PO4-GL1 GL1-GL2 GL1-C1A GL2-C1B"
		link2=" PO4-GL1 GL1-GL2 GL1-D1A GL2-D1B"
		link3=" PO4-AM1 AM1-AM2 AM1-C1A AM2-D1B"
		link4=" PO3-GL1 GL1-GL2 GL1-C1A GL2-C1B"
		link5=" B1A-AM2 AM1-AM2 AM1-C1A AM2-D1B"
		tail_DH=" C1A-C2A C1B-C2B"
		tail_DL=" C1A-C2A C2A-C3A C1B-C2B C2B-C3B"
		tail_DP=" C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B"
		tail_DU=" C1A-D2A D2A-D3A D3A-C4A C1B-D2B D2B-D3B D3B-C4B"
		tail_DS=" C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B"
		tail_DO=" C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B"
		tail_DA=" D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B"
		tail_PO=" C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B"
		tail_PI=" C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B C4B-C5B"
		tail_PP=" C1A-C2A C2A-C3A C3A-C4A D1B-C2B C2B-C3B C3B-C4B"
		
		#define lipids composition
		bond_names['DHPC']=head_PC + link1 + tail_DH
		bond_names['DHPE']=head_PE + link1 + tail_DH
		bond_names['DLPC']=head_PC + link1 + tail_DL
		bond_names['DLPE']=head_PE + link1 + tail_DL
		bond_names['DAPC']=head_PC + link2 + tail_DA
		bond_names['DUPC']=head_PC + link1 + tail_DU
		bond_names['DPPC']=head_PC + link1 + tail_DP
		bond_names['DPPE']=head_PE + link1 + tail_DP
		bond_names['DPPS']=head_PS + link1 + tail_DP
		bond_names['DPPG']=head_PG + link1 + tail_DP
		bond_names['DSPC']=head_PC + link1 + tail_DS
		bond_names['DSPE']=head_PE + link1 + tail_DS
		bond_names['POPC']=head_PC + link1 + tail_PO
		bond_names['POPE']=head_PE + link1 + tail_PO
		bond_names['POPS']=head_PS + link1 + tail_PO
		bond_names['POPG']=head_PG + link1 + tail_PO
		bond_names['PPCS']=head_PC + link3 + tail_PP
		bond_names['PIP2']=head_I2 + link4 + tail_PI
		bond_names['PIP3']=head_I3 + link4 + tail_PI
		bond_names['GM3']=head_GM + link5 + tail_PP
		
		#define tail boundaries (start_A, length_A, start_B, length_B)
		tail_boundaries['DHPC']=[5,1,6,1]
		tail_boundaries['DHPE']=[5,1,6,1]	
		tail_boundaries['DLPC']=[5,2,7,2]
		tail_boundaries['DLPE']=[5,2,7,2]	
		tail_boundaries['DAPC']=[5,4,9,4]
		tail_boundaries['DUPC']=[5,4,9,4]	
		tail_boundaries['DPPC']=[5,3,8,3]
		tail_boundaries['DPPG']=[5,3,8,3]	
		tail_boundaries['DPPE']=[5,3,8,3]	
		tail_boundaries['DPPS']=[5,3,8,3]		
		tail_boundaries['DOPC']=[5,4,9,4]
		tail_boundaries['DOPG']=[5,4,9,4]	
		tail_boundaries['DOPE']=[5,4,9,4]	
		tail_boundaries['DOPS']=[5,4,9,4]		
		tail_boundaries['DSPC']=[5,4,9,4]
		tail_boundaries['DSPE']=[5,4,9,4]
		tail_boundaries['POPC']=[5,3,8,4]
		tail_boundaries['POPE']=[5,3,8,4]
		tail_boundaries['POPS']=[5,3,8,4]
		tail_boundaries['POPG']=[5,3,8,4]
		tail_boundaries['PPCS']=[5,3,8,3]	
		tail_boundaries['PIP2']=[11,3,14,4]
		tail_boundaries['PIP3']=[12,3,15,4]
		tail_boundaries['GM3']=[5,3,8,3]

	return
def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global nb_frames_processed

	if args.xtcfilename=="no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		nb_frames_processed = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes
		nb_frames_processed =0
		U.trajectory.rewind()
	return
def identify_ff():
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index=[]
	lipids_ff_l2u_index=[]
	lipids_sele_ff={}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_beads[args.forcefield_opt] + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		try:
			#read current lipid details
			lip_resname = line.split(',')[0]
			lip_resnum = int(line.split(',')[1])
			lip_leaflet = line.split(',')[2][0:-1]
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + lines[l_index][:-1] + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_proteins():	
	print "\nIdentifying proteins..."
	
	#import modules
	if args.m_algorithm == "density":
		global DBSCAN
		from sklearn.cluster import DBSCAN
	else:
		global nx
		import networkx as nx

	#declare variables
	global proteins_nb
	global proteins_sele
	global proteins_sele_string
	global proteins_sele_string_VMD
	global proteins_boundaries
	global proteins_nb_atoms
	proteins_nb = 0
	proteins_sele = {}
	proteins_sele_string = {}
	proteins_sele_string_VMD = {}
	proteins_boundaries = {}
	
	#check for protein presence
	if U.selectAtoms("protein").numberOfAtoms() == 0:
		print "Error: no protein detected."
		sys.exit(1)
	
	#case: selection file provided
	if args.selection_file_prot!="auto":
		print " -reading protein selection file..."
		with open(args.selection_file_prot) as f:
			lines = f.readlines()
		proteins_nb=len(lines)
		proteins_sele["all"] = MDAnalysis.core.AtomGroup.AtomGroup([])
		for p_index in range(0,proteins_nb):
			line = lines[p_index]
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			try:
				print " p[" + str(p_index) + "]=U.selectAtoms(" + line[0:-1] + ")"
				proteins_sele[p_index] = U.selectAtoms(line[1:-2])
				proteins_sele["all"] += proteins_sele[p_index]
				proteins_boundaries[p_index] = [proteins_sele[p_index].indices()[0] + 1, proteins_sele[p_index].indices()[proteins_sele[p_index].numberOfAtoms()]+1]
				proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
				proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			except:
				print "Error:"
				print line
				print "->invalid selection string."
				sys.exit(1)
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
	
	#case: automatic detection
	else:
		#declare local variables
		proteins_ca_nb = {}
		proteins_ca_nmax = 0
		proteins_ca_group = {}
		proteins_boundaries = {}
	
		#retrieve 1st atom info
		proteins_sele["all"] = U.selectAtoms("protein")
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
		prec_resnum = proteins_sele["all"][0].resnum
		prec_segid = proteins_sele["all"][0].segid
		prec_atnum = proteins_sele["all"][0].number+1
		prev_atnum = proteins_sele["all"][0].number+1						#atom corresponding to the beginning of the current protein
		#browse following atoms
		for a in proteins_sele["all"][1:]:
			delta_res = a.resnum-prec_resnum
			delta_atm = a.number+1-prec_atnum
			if delta_res < 0 or a.segid != prec_segid or delta_atm > 1:
				proteins_boundaries[proteins_nb] = [prev_atnum,prec_atnum]
				proteins_nb += 1
				prev_atnum = a.number + 1
			prec_resnum = a.resnum
			prec_atnum = a.number + 1
			prec_segid = a.segid		
		#add last protein section
		if prev_atnum < proteins_sele["all"][proteins_nb_atoms-1].number:
			proteins_boundaries[proteins_nb] = [prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
			proteins_nb += 1
		
		#display results
		print " -protein found:", proteins_nb
		print " -protein boundaries (atom numbers): see protein.sele file"
		#create protein selections and save into a txt file
		filename_sele=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
		output_stat = open(filename_sele, 'w')	
		output_stat.write("#This file was generated by the script bilayer_perturbations v" + str(version_nb) +"\n")
		output_stat.write("#The lines below correspond to MDAnalysis section string, e.g. U.selectAtoms(LINE)\n")
		output_stat.write("\n")	
		for p_index in range(0, proteins_nb):
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
			proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			proteins_sele[p_index] = U.selectAtoms(proteins_sele_string[p_index])
			output_stat.write(proteins_sele_string[p_index] + "\n")
		output_stat.close()
	print ""

	return
def identify_leaflets():
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	global leaflet_sele_total_nb_atoms
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if numpy.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			leaflet_sele["upper"]["all species"] = L.group(0)
			leaflet_sele["lower"]["all species"] = L.group(1)
		else:
			leaflet_sele["upper"]["all species"] = L.group(1)
			leaflet_sele["lower"]["all species"] = L.group(0)
		leaflet_sele["both"]["all species"] = leaflet_sele["lower"]["all species"] + leaflet_sele["lower"]["all species"]
		if numpy.shape(L.groups())[0]==2:
			print " -found 2 leaflets: ", leaflet_sele["upper"]["all species"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"]["all species"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, numpy.shape(L.groups())[0]):
				other_lipids+=L.group(g).numberOfResidues()
			print " -found " + str(numpy.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"]["all species"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"]["all species"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	else:
		leaflet_sele["both"]["all species"] = U.selectAtoms(leaflet_sele_string)
		tmp_lipids_avg_z = leaflet_sele["both"]["all species"].centerOfGeometry()[2]
		leaflet_sele["upper"]["all species"] = leaflet_sele["both"]["all species"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		leaflet_sele["lower"]["all species"] = leaflet_sele["both"]["all species"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"]["all species"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"]["all species"].numberOfResidues(), '(lower) lipids'
	
	#store full selections
	for l in ["lower","upper","both"]:
		leaflet_sele_atoms[l]["all species"] = leaflet_sele[l]["all species"].residues.atoms
		
	#overal lipids selection
	leaflet_sele_total_nb_atoms = leaflet_sele["both"]["all species"].numberOfAtoms()

	return
def identify_species():
	print "\nIdentifying membrane composition..."
	
	#declare variables
	global op_bonds
	global op_lipids_handled
	global membrane_comp
	global leaflet_ratio
	global leaflet_species
	global lipids_sele_nff
	global lipids_sele_nff_VMD_string
	global lipids_resnums_list
	op_bonds = {}
	op_lipids_handled = {}
	membrane_comp = {}
	leaflet_ratio = {}
	leaflet_species = {}
	lipids_sele_nff = {}
	lipids_sele_nff_VMD_string = {}
	lipids_resnums_list = {}
	
	#specie identification
	for l in ["lower","upper","both"]:
		leaflet_species[l] = list(numpy.unique(leaflet_sele[l]["all species"].resnames()))

	#create specie selections
	for l in ["lower","upper"]:
		for s in leaflet_species[l]:		
			leaflet_sele[l][s] = leaflet_sele[l]["all species"].selectAtoms("resname " + str(s))
			leaflet_sele_atoms[l][s] = leaflet_sele[l][s].residues.atoms
	
	#calculate membrane composition
	for l in ["lower","upper"]:
		membrane_comp[l] = " -" + str(l) + ":"
		leaflet_ratio[l] = {}
		for s in leaflet_species[l]:
			leaflet_ratio[l][s] = round(leaflet_sele[l][s].numberOfResidues()/float(leaflet_sele[l]["all species"].numberOfResidues())*100,1)
			membrane_comp[l] += " " + s + " (" + str(leaflet_ratio[l][s]) + "%)"
		print membrane_comp[l]
	
	#create individual lipid selections
	for l in ["lower","upper"]:
		lipids_sele_nff[l] = {}
		lipids_sele_nff_VMD_string[l] = {}
		for s in leaflet_species[l]:
			lipids_sele_nff[l][s] = {}
			lipids_sele_nff_VMD_string[l][s] = {}
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				lipids_sele_nff[l][s][r_index] = leaflet_sele[l][s].selectAtoms("resnum " + str(leaflet_sele[l][s].resnums()[r_index])).residues.atoms
				lipids_sele_nff_VMD_string[l][s][r_index] = "resname " + str(s) + " and resid " + str(leaflet_sele[l][s].resnums()[r_index])
	
	#case: need to calculate order parameters
	if args.perturb == 2 or args.perturb == 3:
		#create list of lipids to take into account for order parameters calculation
		for l in ["lower","upper"]:
			op_lipids_handled[l] = []
			leaflet_sele[l]["op"] = MDAnalysis.core.AtomGroup.AtomGroup([])
			for s in leaflet_species[l]:
				if s in op_lipids_possible[args.forcefield_opt]:
					op_lipids_handled[l].append(s)
					leaflet_sele[l]["op"] += leaflet_sele[l][s]
		op_lipids_handled["both"]= list(numpy.unique(op_lipids_handled["lower"] + op_lipids_handled["upper"]))
		leaflet_sele["both"]["op"] = leaflet_sele["lower"]["op"] + leaflet_sele["upper"]["op"]
		if len(op_lipids_handled["both"]) == 0:
			print "Error: none of the lipid species can be taken into account - double check the forcefield option, see order_param -h."
			sys.exit(1)
		else:
			membrane_comp["op"] = "\nLipids handled for order parameters calculations: "
			for s in op_lipids_handled["both"]:
				membrane_comp["op"] += str(s) + ", "
			membrane_comp["op"] = membrane_comp["op"][:-2]
			print membrane_comp["op"]
			
		#create list of bonds for each lipid specie
		for s in op_lipids_handled["both"]:
			op_bonds[s]=[]
			for bond_name in bond_names[s].split():
				op_bonds[s].append(bond_name.split("-"))

		#check that all the necessary lipids tails particles are present in order to calculate order parameter
		for l in ["lower","upper"]:
			for s in op_lipids_handled[l]:
				tail_A_start = tail_boundaries[s][0]
				for bond in op_bonds[s][tail_A_start:]:
					tmp_sele_b0 = leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[0]))
					tmp_sele_b1 = leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[1]))
					if tmp_sele_b0.numberOfAtoms() == 0:
						print "Error: missing particles " + str(bond[0]) + " in " + str(s) + " lipids in the " + str(l) + " leaflet."
						print "->cannot calculate tail order parameters."
						sys.exit(1)
					if tmp_sele_b1.numberOfAtoms() == 0:
						print "Error: missing particles " + str(bond[1]) + " in " + str(s) + " lipids in the " + str(l) + " leaflet."
						print "->cannot calculate tail order parameters, choose a different --perturb option or use different inputs."
						sys.exit(1)

	#create list of residues for each lipid specie		
	for l in ["lower","upper"]:
		lipids_resnums_list[l] = {}
		lipids_resnums_list[l]["all species"] = list(leaflet_sele[l]["all species"].resnums())
		for s in leaflet_species[l]:
			lipids_resnums_list[l][s] = list(leaflet_sele[l][s].resnums())
		if args.perturb == 2 or args.perturb == 3:
			lipids_resnums_list[l]["op"] = list(leaflet_sele[l]["op"].resnums())
	
	return
def initialise_colours_and_groups():
	#declare variables
	global colours_lipids
	global colours_lipids_nb
	global colours_lipids_map
	colours_lipids = {}
	colours_lipids_nb = 0
	colours_lipids_map = "jet"
	colormaps_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']

	#colours: lipid species
	#======================
	#process colour file if provided
	#-------------------------------
	if args.colour_file != "no":
		#read file
		print "\nReading lipids colour definition file..."
		with open(args.colour_file) as f:
			lines = f.readlines()
		colours_lipids_nb = len(lines)
		for line_index in range(0,colours_lipids_nb):
			l_content = lines[line_index].split(',')
			colours_lipids[l_content[0]] = l_content[1][:-1]
		print " -found the following colours definition:"
		for s in colours_lipids.keys():
			print " -" + str(s) + ": " + str(colours_lipids[s])
		
		#check whether a colour map was specified
		if colours_lipids_nb > 1 and len(numpy.unique(colours_lipids.values())) == 1:
			if numpy.unique(colours_lipids.values())[0] in colormaps_possible:
				colours_lipids_map = numpy.unique(colours_lipids.values())[0]
			else:
				print "Error: either the same color was specified for all species or the color map '" + str(numpy.unique(colours_lipids.values())[0]) + "' is not valid."
				sys.exit(1)
		else:
			colours_lipids_map = "custom"
			
		#check that all detected species have a colour specified
		for s in leaflet_species["both"]:
			if s not in colours_lipids.keys():
				print "Error: no colour specified for " + str(s) + "."
				sys.exit(1)
	
	#generate colours from colour map if necessary
	#---------------------------------------------
	if colours_lipids_map != "custom":
		tmp_cmap = cm.get_cmap(colours_lipids_map)
		colours_lipids_value = tmp_cmap(numpy.linspace(0, 1, len(leaflet_species["both"])))
		for s_index in range(0, len(leaflet_species["both"])):
			colours_lipids[leaflet_species["both"][s_index]]=colours_lipids_value[s_index]

	#colours: size groups
	#====================
	if args.cluster_groups_file != "no":
		global colours_groups
		global groups_number
		global groups_boundaries
		global groups_sizes_dict
		colours_groups = {}
		colours_groups_map = "custom"
				
		#read file
		print "\nReading cluster groups definition file..."
		with open(args.cluster_groups_file) as f:
			lines = f.readlines()
		groups_number = len(lines)
		for g_index in range(0,groups_number):
			l_content=lines[g_index].split(',')
			tmp_beg = int(l_content[0])
			tmp_end = l_content[1]
			colours_groups[g_index] = l_content[2][:-1]					#attribute colour to group
			if tmp_end == "max":
				tmp_end = 100000										#put a stupidly big size to cap the open ended group
			else:
				tmp_end = int(tmp_end)
			groups_boundaries.append([tmp_beg,tmp_end])
			
		#display results
		print " -found " + str(groups_number) + " cluster groups:"
		for g_index in range(0,groups_number):
			if groups_boundaries[g_index][1]==100000:
				print "   g" + str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "+, " + str(colours_groups[g_index])
			else:
				print "   g" + str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + ", " + str(colours_groups[g_index])

		#check for boundaries overlapping
		prev_beg = groups_boundaries[0][0]
		prev_end = groups_boundaries[0][1]
		if prev_end < prev_beg:
			print "Error: the max size is smaller than the min size for specified cluster groups " + str(groups_boundaries[0]) + "."
			sys.exit(1)
		for g in groups_boundaries[1:]:
			if g[1] < g[0]:
				print "Error: the max size is smaller than the min size for group " + str(g) + "."
				sys.exit(1)
			if g[0] <= prev_end:
				print "Error: specified cluster groups [" + str(prev_beg) + "," + str(prev_end) + "] and " + str(g) + " overlap or are not in increasing order."
				sys.exit(1)
			prev_beg=g[0]
			prev_end=g[1]
		
		#create equivalency table between groups and sizes
		for g_index in range(0,groups_number):
			bb = groups_boundaries[g_index]
			tmp_beg = bb[0]
			tmp_end = bb[1]
			for tmp_size in range(tmp_beg, tmp_end+1):
				groups_sizes_dict[tmp_size] = g_index
		for tmp_size in list(set(range(1,max(groups_sizes_dict.keys()))) - set(groups_sizes_dict.keys())): 	#this handles potentially unaccounted for sizes up to the maximum specified by the user
			groups_sizes_dict[tmp_size] = groups_number
		if max(groups_sizes_dict.keys())!= 100000:															    #this handles potentially unaccounted for sizes above the maximum specified by the user (in case it's not an open group)
			for tmp_size in range(max(groups_sizes_dict.keys())+1,100001):
				groups_sizes_dict[tmp_size] = groups_number

		#check whether a colour map was specified
		if groups_number > 1 and len(numpy.unique(colours_groups.values())) == 1:
			if numpy.unique(colours_groups.values())[0] in colormaps_possible:
				colours_groups_map = numpy.unique(colours_groups.values())[0]
			else:
				print "Error: either the same color was specified for all groups or the color map '" + str(numpy.unique(colours_groups.values())[0]) + "' is not valid."
				sys.exit(1)

		#generate colours from colour map if necessary
		if groups_colors_map != "custom":
			tmp_cmap=cm.get_cmap(groups_colors_map)
			groups_colors_value = tmp_cmap(numpy.linspace(0, 1, groups_number))
			for g_index in range(0, groups_number):
				colours_groups[g_index] = groups_colors_value[g_index]
		
		#add the group "other" 
		colours_groups[groups_number] = "#C0C0C0"
		
	return

################################################################################################################################################
# FUNCTIONS: data structures
################################################################################################################################################

def data_struct_time_coords():

	global time_stamp
	global z_lower
	global z_upper
	global z_ff

	time_stamp = {}
	z_lower = {}
	z_upper = {}
	z_ff = {}
	
	for l_index in range(0,lipids_ff_nb):
		z_ff[l_index] = {}

	return
def data_struct_op_ff():

	#order parameter: full data
	global lipids_op_ff_tailA
	global lipids_op_ff_tailB
	global lipids_op_ff_tails

	#order paramater: sorted data
	global lipids_op_ff_tailA_sorted
	global lipids_op_ff_tailB_sorted
	global lipids_op_ff_tails_sorted

	#order paramater: smoothed data
	global lipids_op_ff_tailA_smooth
	global lipids_op_ff_tailB_smooth
	global lipids_op_ff_tails_smooth
	
	lipids_op_ff_tailA = {}
	lipids_op_ff_tailB = {}
	lipids_op_ff_tails = {}
	lipids_op_ff_tailA_sorted = {}
	lipids_op_ff_tailB_sorted = {}
	lipids_op_ff_tails_sorted = {}
	lipids_op_ff_tailA_smooth = {}
	lipids_op_ff_tailB_smooth = {}
	lipids_op_ff_tails_smooth = {}

	for l_index in range(0,lipids_ff_nb):
		lipids_op_ff_tailA[l_index] = {}
		lipids_op_ff_tailB[l_index] = {}
		lipids_op_ff_tails[l_index] = {}
		lipids_op_ff_tailA_sorted[l_index] = []
		lipids_op_ff_tailB_sorted[l_index] = []
		lipids_op_ff_tails_sorted[l_index] = []
		lipids_op_ff_tailA[l_index]["all frames"] = []
		lipids_op_ff_tailB[l_index]["all frames"] = []
		lipids_op_ff_tails[l_index]["all frames"] = []

	return
def data_struct_op_nff():
	
	global lipids_op_nff
	lipids_op_nff = {}

	#op: each lipid
	for l in ["lower","upper"]:
		lipids_op_nff[l] = {}
		for s in op_lipids_handled[l]:
			lipids_op_nff[l][s] = {}
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				lipids_op_nff[l][s][r_index] = {}
	
	#op: each specie
	lipids_op_nff["data"] = {}
	lipids_op_nff["sorted"] = {}
	lipids_op_nff["sorted"]["avg"] = {}
	lipids_op_nff["sorted"]["std"] = {}
	lipids_op_nff["smoothed"] = {}
	lipids_op_nff["smoothed"]["avg"] = {}
	lipids_op_nff["smoothed"]["std"] = {}
	for tail in ["tailA", "tailB", "tails"]:
		lipids_op_nff["data"][tail] = {}
		lipids_op_nff["sorted"]["avg"][tail] = {}	
		lipids_op_nff["sorted"]["std"][tail] = {}		
		lipids_op_nff["smoothed"]["avg"][tail] = {}
		lipids_op_nff["smoothed"]["std"][tail] = {}
		for l in ["lower", "upper"]:
			lipids_op_nff["data"][tail][l] = {}
			lipids_op_nff["sorted"]["avg"][tail][l] = {}	
			lipids_op_nff["sorted"]["std"][tail][l] = {}	
			lipids_op_nff["smoothed"]["avg"][tail][l] = {}
			lipids_op_nff["smoothed"]["std"][tail][l] = {}
			for s in op_lipids_handled[l] + ["all species"]:
				lipids_op_nff["data"][tail][l][s] = {}
				lipids_op_nff["data"][tail][l][s]["all frames"] = []
				lipids_op_nff["sorted"]["avg"][tail][l][s] = []
				lipids_op_nff["sorted"]["std"][tail][l][s] = []	
	
	return
def data_struct_thick():

	global lipids_thick_nff
	lipids_thick_nff = {}
	
	#thickness: each lipid
	for l in ["lower","upper"]:
		lipids_thick_nff[l] = {}
		for s in leaflet_species[l]:
			lipids_thick_nff[l][s] = {}
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				lipids_thick_nff[l][s][r_index] = {}

	#thickness: each specie
	lipids_thick_nff["data"] = {}
	lipids_thick_nff["sorted"] = {}
	lipids_thick_nff["sorted"]["avg"] = {}
	lipids_thick_nff["sorted"]["std"] = {}
	lipids_thick_nff["smoothed"] = {}
	lipids_thick_nff["smoothed"]["avg"] = {}
	lipids_thick_nff["smoothed"]["std"] = {}
	for s in leaflet_species["both"] + ["all species"]:
		lipids_thick_nff["data"][s] = {}
		lipids_thick_nff["data"][s]["all frames"] = []
		lipids_thick_nff["sorted"]["avg"][s] = []
		lipids_thick_nff["sorted"]["std"][s] = []
	
	return
def data_struct_radial():
		
	#density
	#-------
	global radial_density
	radial_density = {}	
	for l in ["lower","upper"]:
		radial_density[l] = {}
		#all species
		radial_density[l]["all species"] = {}
		radial_density[l]["all species"]["all sizes"] = {}
		radial_density[l]["all species"]["all sizes"]["nb"] = {}
		radial_density[l]["all species"]["all sizes"]["pc"] = {}
		radial_density[l]["all species"]["all sizes"]["nb"]["all frames"] = numpy.zeros(args.radial_nb_bins)
		radial_density[l]["all species"]["all sizes"]["pc"]["all frames"] = numpy.zeros(args.radial_nb_bins)
		#each specie
		for s in leaflet_species[l]:
			radial_density[l][s] = {}
			radial_density[l][s]["all sizes"] = {}
			radial_density[l][s]["all sizes"]["nb"] = {}
			radial_density[l][s]["all sizes"]["pc"] = {}
			radial_density[l][s]["all sizes"]["nb"]["all frames"] = numpy.zeros(args.radial_nb_bins)
			radial_density[l][s]["all sizes"]["pc"]["all frames"] = numpy.zeros(args.radial_nb_bins)

	#thickness
	#---------
	if args.perturb == 1 or args.perturb == 3:
		global radial_thick
		radial_thick = {}	
		#all species
		radial_thick["all species"] = {}
		radial_thick["all species"]["avg"] = {}
		radial_thick["all species"]["std"] = {}
		radial_thick["all species"]["avg"]["all sizes"] = {}
		radial_thick["all species"]["std"]["all sizes"] = {}
		radial_thick["all species"]["avg"]["all sizes"]["all frames"] = {}
		radial_thick["all species"]["std"]["all sizes"]["all frames"] = {}
		#each specie
		for s in leaflet_species["both"]:
			radial_thick[s] = {}
			radial_thick[s]["avg"] = {}
			radial_thick[s]["std"] = {}
			radial_thick[s]["avg"]["all sizes"] = {}
			radial_thick[s]["std"]["all sizes"] = {}
			radial_thick[s]["avg"]["all sizes"]["all frames"] = {}
			radial_thick[s]["std"]["all sizes"]["all frames"] = {}
		for n in range(0,args.radial_nb_bins):
			radial_thick["all species"]["avg"]["all sizes"]["all frames"][n] = []
			radial_thick["all species"]["std"]["all sizes"]["all frames"][n] = []
			for s in leaflet_species["both"]:
				radial_thick[s]["avg"]["all sizes"]["all frames"][n] = []
				radial_thick[s]["std"]["all sizes"]["all frames"][n] = []

	#order parameter
	#---------------
	if args.perturb == 2 or args.perturb == 3:
		global radial_op
		radial_op = {}
		for l in ["lower","upper"]:
			radial_op[l] = {}
			#all species
			radial_op[l]["all species"] = {}
			radial_op[l]["all species"]["avg"] = {}
			radial_op[l]["all species"]["std"] = {}
			radial_op[l]["all species"]["avg"]["all sizes"] = {}
			radial_op[l]["all species"]["std"]["all sizes"] = {}
			radial_op[l]["all species"]["avg"]["all sizes"]["all frames"] = {}
			radial_op[l]["all species"]["std"]["all sizes"]["all frames"] = {}
			#each specie
			for s in op_lipids_handled[l]:
				radial_op[l][s] = {}
				radial_op[l][s]["avg"] = {}
				radial_op[l][s]["std"] = {}
				radial_op[l][s]["avg"]["all sizes"] = {}
				radial_op[l][s]["std"]["all sizes"] = {}
				radial_op[l][s]["avg"]["all sizes"]["all frames"] = {}
				radial_op[l][s]["std"]["all sizes"]["all frames"] = {}
			for n in range(0,args.radial_nb_bins):
				radial_op[l]["all species"]["avg"]["all sizes"]["all frames"][n] = []
				radial_op[l]["all species"]["std"]["all sizes"]["all frames"][n] = []
				for s in op_lipids_handled[l]:
					radial_op[l][s]["avg"]["all sizes"]["all frames"][n] = []
					radial_op[l][s]["std"]["all sizes"]["all frames"][n] = []

	return

################################################################################################################################################
# FUNCTIONS: core
################################################################################################################################################

def get_z_coords(f_nb):

	#to do include relevant beads selection string
	z_middle_instant = leaflet_sele_atoms["lower"]["all species"].selectAtoms("name PO4").centerOfGeometry()[2]+(leaflet_sele_atoms["upper"]["all species"].selectAtoms("name PO4").centerOfGeometry()[2]-leaflet_sele_atoms["lower"]["all species"].selectAtoms("name PO4").centerOfGeometry()[2])/float(2)
	z_lower[f_nb] = leaflet_sele_atoms["lower"]["all species"].selectAtoms("name PO4").centerOfGeometry()[2]-z_middle_instant
	z_upper[f_nb] = leaflet_sele_atoms["upper"]["all species"].selectAtoms("name PO4").centerOfGeometry()[2]-z_middle_instant
	for l in range(0,lipids_ff_nb):	
		z_ff[l][f_nb] = lipids_sele_ff[l].selectAtoms("name PO4").centerOfGeometry()[2]-z_middle_instant

	return
def get_distances(box_dim):
	
	#method: use minimum distance between proteins
	#---------------------------------------------
	if args.m_algorithm=="min":
		dist_matrix=100000*numpy.ones((proteins_nb,proteins_nb))
		for n in range(proteins_nb,1,-1):
			p_index=proteins_nb-n
			for pp in range(p_index+1,proteins_nb):
				dist_matrix[p_index,pp]=numpy.min(MDAnalysis.analysis.distances.distance_array(numpy.float32(proteins_sele[p_index].coordinates()),numpy.float32(proteins_sele[pp].coordinates()),box_dim))	
				dist_matrix[pp,p_index]=dist_matrix[p_index,pp]
											
	#method: use distance between cog
	#--------------------------------
	else:
		#store cog coord
		tmp_proteins_coords=numpy.zeros((proteins_nb,3))
		for p_index in range(0,proteins_nb):
			tmp_proteins_coords[p_index,:]=proteins_sele[p_index].centerOfGeometry()
		#calculate distances
		dist_matrix=MDAnalysis.analysis.distances.distance_array(numpy.float32(tmp_proteins_coords), numpy.float32(tmp_proteins_coords), box_dim)

	return dist_matrix
def calculate_cog(sele):				#TO UPDATE
	
	#to do:
	# move all particles into the same quadrant (translate by the appropriate box dim)
	# calc center of geometry
	# check it's in the box and put it back in if necessary
	
	tmp_cog=sele.centerOfGeometry()
	
	return tmp_cog
def calculate_radial(f_nb):
	
	global radial_step
	
	#detect protein clusters
	#=======================
	#identify clusters
	if args.m_algorithm!="density":
		tmp_groups = detect_clusters_connectivity(get_distances(U.trajectory.ts.dimensions), U.trajectory.ts.dimensions)
	else:
		tmp_groups = detect_clusters_density(get_distances(U.trajectory.ts.dimensions), U.trajectory.ts.dimensions)
	
	#store cluster size and selection
	radial_sizes[f_nb] = []
	tmp_cluster_selections = {}
	for g in tmp_groups:		
		#create selection for current cluster
		tmp_cluster_sele = MDAnalysis.core.AtomGroup.AtomGroup([])	
		for p_index in g:
			tmp_cluster_sele += proteins_sele[p_index]
		
		#check whether the cluster is TM
		#find closest PO4 particles for each particles of clusters, if all are in the same leaflet then it's surfacic [NB: this is done at the CLUSTER level (the same criteria at the protein level would probably fail)]
		dist_min_lower = numpy.min(MDAnalysis.analysis.distances.distance_array(tmp_cluster_sele.coordinates(), leaflet_sele["lower"]["all species"].coordinates(), U.trajectory.ts.dimensions),axis=1)
		dist_min_upper = numpy.min(MDAnalysis.analysis.distances.distance_array(tmp_cluster_sele.coordinates(), leaflet_sele["upper"]["all species"].coordinates(), U.trajectory.ts.dimensions),axis=1)
		dist = dist_min_upper-dist_min_lower

		#store current cluster details if it is a TM cluster
		if numpy.size(dist[dist>0]) != numpy.size(dist) and numpy.size(dist[dist>0]) !=0:
			if numpy.size(g) not in radial_sizes[f_nb]:
				radial_sizes[f_nb].append(numpy.size(g))
				tmp_cluster_selections[numpy.size(g)]=[]
			tmp_cluster_selections[numpy.size(g)].append(tmp_cluster_sele)
		
	#update data structures
	#======================
	radial_sizes[f_nb] = sorted(radial_sizes[f_nb])
	for c_size in radial_sizes[f_nb]:
		#add size if necessary
		#---------------------
		if c_size not in radial_sizes["all frames"]:
			radial_sizes["all frames"].append(c_size)
			radial_sizes["all frames"] = sorted(radial_sizes["all frames"])
			#density
			for l in ["lower","upper"]:	
				radial_density[l]["all species"][c_size] = {}
				radial_density[l]["all species"][c_size]["nb"] = {}
				radial_density[l]["all species"][c_size]["pc"] = {}
				radial_density[l]["all species"][c_size]["nb"]["all frames"] = numpy.zeros(args.radial_nb_bins)
				radial_density[l]["all species"][c_size]["pc"]["all frames"] = numpy.zeros(args.radial_nb_bins)
				for s in leaflet_species[l]:
					radial_density[l][s][c_size] = {}
					radial_density[l][s][c_size]["nb"] = {}
					radial_density[l][s][c_size]["pc"] = {}
					radial_density[l][s][c_size]["nb"]["all frames"] = numpy.zeros(args.radial_nb_bins)
					radial_density[l][s][c_size]["pc"]["all frames"] = numpy.zeros(args.radial_nb_bins)
				
			#thickness
			if args.perturb == 1 or args.perturb == 3:
				radial_thick["all species"]["avg"][c_size] = {}
				radial_thick["all species"]["std"][c_size] = {}
				radial_thick["all species"]["avg"][c_size]["all frames"] = {}
				radial_thick["all species"]["std"][c_size]["all frames"] = {}
				for s in leaflet_species["both"]:
					radial_thick[s]["avg"][c_size] = {}
					radial_thick[s]["std"][c_size] = {}
					radial_thick[s]["avg"][c_size]["all frames"] = {}
					radial_thick[s]["std"][c_size]["all frames"] = {}
				for n in range(0,args.radial_nb_bins):
					radial_thick["all species"]["avg"][c_size]["all frames"][n] = []
					radial_thick["all species"]["std"][c_size]["all frames"][n] = []
					for s in leaflet_species["both"]:
						radial_thick[s]["avg"][c_size]["all frames"][n] = []
						radial_thick[s]["std"][c_size]["all frames"][n] = []
				
			#order parameter
			if args.perturb == 2 or args.perturb == 3:
				for l in ["lower","upper"]:					
					radial_op[l]["all species"]["avg"][c_size] = {}
					radial_op[l]["all species"]["std"][c_size] = {}
					radial_op[l]["all species"]["avg"][c_size]["all frames"] = {}
					radial_op[l]["all species"]["std"][c_size]["all frames"] = {}
					for s in op_lipids_handled[l]:
						radial_op[l][s]["avg"][c_size] = {}
						radial_op[l][s]["std"][c_size] = {}
						radial_op[l][s]["avg"][c_size]["all frames"] = {}
						radial_op[l][s]["std"][c_size]["all frames"] = {}
					for n in range(0,args.radial_nb_bins):
						radial_op[l]["all species"]["avg"][c_size]["all frames"][n] = []
						radial_op[l]["all species"]["std"][c_size]["all frames"][n] = []
						for s in op_lipids_handled[l]:
							radial_op[l][s]["avg"][c_size]["all frames"][n] = []
							radial_op[l][s]["std"][c_size]["all frames"][n] = []
	
		#add frame entry
		#---------------
		#density
		for l in ["lower","upper"]:
			radial_density[l]["all species"][c_size]["nb"][f_nb] = numpy.zeros(args.radial_nb_bins)
			radial_density[l]["all species"][c_size]["pc"][f_nb] = numpy.zeros(args.radial_nb_bins)
			radial_density[l]["all species"]["all sizes"]["nb"][f_nb] = numpy.zeros(args.radial_nb_bins)
			radial_density[l]["all species"]["all sizes"]["pc"][f_nb] = numpy.zeros(args.radial_nb_bins)
			for s in leaflet_species[l]:
				radial_density[l][s][c_size]["nb"][f_nb] = numpy.zeros(args.radial_nb_bins)
				radial_density[l][s][c_size]["pc"][f_nb] = numpy.zeros(args.radial_nb_bins)
				radial_density[l][s]["all sizes"]["nb"][f_nb] = numpy.zeros(args.radial_nb_bins)
				radial_density[l][s]["all sizes"]["pc"][f_nb] = numpy.zeros(args.radial_nb_bins)
		#thickness
		if args.perturb == 1 or args.perturb == 3:
			radial_thick["all species"]["avg"][c_size][f_nb] = {}
			radial_thick["all species"]["std"][c_size][f_nb] = {}
			radial_thick["all species"]["avg"]["all sizes"][f_nb] = {}
			radial_thick["all species"]["std"]["all sizes"][f_nb] = {}
			for s in leaflet_species["both"]:
				radial_thick[s]["avg"][c_size][f_nb] = {}
				radial_thick[s]["std"][c_size][f_nb] = {}
				radial_thick[s]["avg"]["all sizes"][f_nb] = {}
				radial_thick[s]["std"]["all sizes"][f_nb] = {}
			for n in range(0,args.radial_nb_bins):
				radial_thick["all species"]["avg"][c_size][f_nb][n] = []
				radial_thick["all species"]["std"][c_size][f_nb][n] = []
				radial_thick["all species"]["avg"]["all sizes"][f_nb][n] = []
				radial_thick["all species"]["std"]["all sizes"][f_nb][n] = []
				for s in leaflet_species["both"]:
					radial_thick[s]["avg"][c_size][f_nb][n] = []
					radial_thick[s]["std"][c_size][f_nb][n] = []
					radial_thick[s]["avg"]["all sizes"][f_nb][n] = []
					radial_thick[s]["std"]["all sizes"][f_nb][n] = []
		#order param
		if args.perturb == 2 or args.perturb == 3:
			for l in ["lower","upper"]:
				radial_op[l]["all species"]["avg"][c_size][f_nb] = {}
				radial_op[l]["all species"]["std"][c_size][f_nb] = {}
				radial_op[l]["all species"]["avg"]["all sizes"][f_nb] = {}
				radial_op[l]["all species"]["std"]["all sizes"][f_nb] = {}
				for s in op_lipids_handled[l]:
					radial_op[l][s]["avg"][c_size][f_nb] = {}
					radial_op[l][s]["std"][c_size][f_nb] = {}
					radial_op[l][s]["avg"]["all sizes"][f_nb] = {}
					radial_op[l][s]["std"]["all sizes"][f_nb] = {}
				for n in range(0,args.radial_nb_bins):
					radial_op[l]["all species"]["avg"][c_size][f_nb][n] = []
					radial_op[l]["all species"]["std"][c_size][f_nb][n] = []
					radial_op[l]["all species"]["avg"]["all sizes"][f_nb][n] = []
					radial_op[l]["all species"]["std"]["all sizes"][f_nb][n] = []
					for s in op_lipids_handled[l]:
						radial_op[l][s]["avg"][c_size][f_nb][n] = []
						radial_op[l][s]["std"][c_size][f_nb][n] = []
						radial_op[l][s]["avg"]["all sizes"][f_nb][n] = []
						radial_op[l][s]["std"]["all sizes"][f_nb][n] = []
	
	#deal with lipids around them
	#============================
	for l in ["lower","upper"]:
		tmp_resnames = leaflet_sele[l]["all species"].resnames()
		tmp_resnums = numpy.zeros((1,leaflet_sele[l]["all species"].numberOfResidues()))
		tmp_resnums[0,:] = leaflet_sele[l]["all species"].resnums()
		for c_size in radial_sizes[f_nb]:
			for c_sele in tmp_cluster_selections[c_size]:
				#retrieve COG coords of current cluster
				tmp_c_cog = numpy.zeros((1,3))
				tmp_c_cog[0,:] = calculate_cog(c_sele)

				#calculate distance matrix between lipids and cluster cog and retrieve index of lipids within cutoff
				lip_dist = MDAnalysis.analysis.distances.distance_array(numpy.float32(tmp_c_cog), leaflet_sele[l]["all species"].selectAtoms(leaflet_beads[args.forcefield_opt]).coordinates(), U.trajectory.ts.dimensions)
				r_num_within = list(tmp_resnums[lip_dist<args.radial_radius])

				#drop data into the right radial bin for those lipids
				#----------------------------------------------------
				#case: density only
				if args.perturb == 0:
					for r_numf in r_num_within:
						r_num = int(r_numf)
						r_index = lipids_resnums_list[l]["all species"].index(r_num)
						r_specie = tmp_resnames[r_index]
						r_index_specie = lipids_resnums_list[l][r_specie].index(r_num)
						r_bin = int(numpy.floor(lip_dist[0,r_index]/float(radial_step)))
						
						#density
						radial_density[l][r_specie][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"]["all frames"][r_bin] += 1
				
				#case: density + thickness
				elif args.perturb == 1:
					for r_numf in r_num_within:
						r_num = int(r_numf)
						r_index = lipids_resnums_list[l]["all species"].index(r_num)
						r_specie = tmp_resnames[r_index]
						r_index_specie = lipids_resnums_list[l][r_specie].index(r_num)
						r_bin = int(numpy.floor(lip_dist[0,r_index]/float(radial_step)))
						
						#density
						radial_density[l][r_specie][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"]["all frames"][r_bin] += 1

						#thickness
						radial_thick["all species"]["avg"][c_size][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick["all species"]["avg"]["all sizes"][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"][c_size][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"]["all sizes"][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick["all species"]["avg"][c_size]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick["all species"]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"][c_size]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])

				#case: density + order parameters
				elif args.perturb == 2:
					for r_numf in r_num_within:
						r_num = int(r_numf)
						r_index = lipids_resnums_list[l]["all species"].index(r_num)
						r_specie = tmp_resnames[r_index]
						r_index_specie = lipids_resnums_list[l][r_specie].index(r_num)
						r_bin = int(numpy.floor(lip_dist[0,r_index]/float(radial_step)))

						#density
						radial_density[l][r_specie][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"]["all frames"][r_bin] += 1
					
						#order parameters
						if r_specie in op_lipids_handled[l]:
							radial_op[l][r_specie]["avg"][c_size][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l][r_specie]["avg"]["all sizes"][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"][c_size][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"]["all sizes"][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l][r_specie]["avg"][c_size]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l][r_specie]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"][c_size]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
				
				#case: density + thickness + order parameters
				else:
					for r_numf in r_num_within:
						r_num = int(r_numf)
						r_index = lipids_resnums_list[l]["all species"].index(r_num)
						r_specie = tmp_resnames[r_index]
						r_index_specie = lipids_resnums_list[l][r_specie].index(r_num)
						r_bin = int(numpy.floor(lip_dist[0,r_index]/float(radial_step)))
						
						#density
						radial_density[l][r_specie][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"][f_nb][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"][f_nb][r_bin] += 1
						radial_density[l][r_specie][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l][r_specie]["all sizes"]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"][c_size]["nb"]["all frames"][r_bin] += 1
						radial_density[l]["all species"]["all sizes"]["nb"]["all frames"][r_bin] += 1
						
						#thickness
						radial_thick["all species"]["avg"][c_size][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick["all species"]["avg"]["all sizes"][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"][c_size][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"]["all sizes"][f_nb][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick["all species"]["avg"][c_size]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick["all species"]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"][c_size]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])
						radial_thick[r_specie]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_thick_nff[l][r_specie][r_index_specie][f_nb])

						#order parameters
						if r_specie in op_lipids_handled[l]:
							radial_op[l][r_specie]["avg"][c_size][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l][r_specie]["avg"]["all sizes"][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"][c_size][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"]["all sizes"][f_nb][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l][r_specie]["avg"][c_size]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l][r_specie]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"][c_size]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
							radial_op[l]["all species"]["avg"]["all sizes"]["all frames"][r_bin].append(lipids_op_nff[l][r_specie][r_index_specie][f_nb])
			
	return
def calculate_radial_data():
	
	#case: density:
	for l in ["lower","upper"]:
		for s in leaflet_species[l]:
			for c_size in radial_sizes["all frames"] + ["all sizes"]:
				for f_nb in radial_density[l][s][c_size]["nb"].keys():			#f_nb can be "all frames"...
					for n in range(0,args.radial_nb_bins):								
						if radial_density[l]["all species"][c_size]["nb"][f_nb][n] != 0:
							radial_density[l][s][c_size]["pc"][f_nb][n] = radial_density[l][s][c_size]["nb"][f_nb][n] / float(radial_density[l]["all species"][c_size]["nb"][f_nb][n])*100
											
	#case: thickness
	#===============
	if args.perturb == 1 or args.perturb == 3:		
		for s in leaflet_species["both"] + ["all species"]:
			for c_size in radial_sizes["all frames"] + ["all sizes"]:
				for f_nb in radial_thick[s]["avg"][c_size].keys():			#f_nb can be "all frames"...
					for n in range(0,args.radial_nb_bins):			
						if s != "all species":
							if (s in leaflet_species["lower"] and s in leaflet_species["upper"]) and (radial_density["lower"][s][c_size]["nb"][f_nb][n] == 0 and radial_density["upper"][s][c_size]["nb"][f_nb][n] == 0):
								radial_thick[s]["avg"][c_size][f_nb][n] = numpy.nan
								radial_thick[s]["std"][c_size][f_nb][n] = numpy.nan
							elif (s in leaflet_species["lower"] and s not in leaflet_species["upper"]) and radial_density["lower"][s][c_size]["nb"][f_nb][n] == 0 :
								radial_thick[s]["avg"][c_size][f_nb][n] = numpy.nan
								radial_thick[s]["std"][c_size][f_nb][n] = numpy.nan
							elif (s not in leaflet_species["lower"] and s in leaflet_species["upper"]) and radial_density["upper"][s][c_size]["nb"][f_nb][n] == 0 :
								radial_thick[s]["avg"][c_size][f_nb][n] = numpy.nan
								radial_thick[s]["std"][c_size][f_nb][n] = numpy.nan
							else:
								tmp_avg = numpy.average(radial_thick[s]["avg"][c_size][f_nb][n])
								tmp_std = numpy.std(radial_thick[s]["avg"][c_size][f_nb][n])
								radial_thick[s]["avg"][c_size][f_nb][n] = tmp_avg
								radial_thick[s]["std"][c_size][f_nb][n] = tmp_std
						else:
							if radial_density["lower"]["all species"][c_size]["nb"][f_nb][n] == 0 and radial_density["upper"]["all species"][c_size]["nb"][f_nb][n] == 0:
								radial_thick[s]["avg"][c_size][f_nb][n] = numpy.nan
								radial_thick[s]["std"][c_size][f_nb][n] = numpy.nan
							else:
								tmp_avg = numpy.average(radial_thick[s]["avg"][c_size][f_nb][n])
								tmp_std = numpy.std(radial_thick[s]["avg"][c_size][f_nb][n])
								radial_thick[s]["avg"][c_size][f_nb][n] = tmp_avg
								radial_thick[s]["std"][c_size][f_nb][n] = tmp_std
							
	
	#case: order parameters
	#======================
	if args.perturb == 2 or args.perturb == 3:
		for l in ["lower","upper"]:
			for s in op_lipids_handled[l] + ["all species"]:
				for c_size in radial_sizes["all frames"] + ["all sizes"]:
					for f_nb in radial_op[l][s]["avg"][c_size].keys():			#f_nb can be "all frames"...
						for n in range(0,args.radial_nb_bins):			
							if radial_density[l][s][c_size]["nb"][f_nb][n] == 0:
								radial_op[l][s]["avg"][c_size][f_nb][n] = numpy.nan
								radial_op[l][s]["std"][c_size][f_nb][n] = numpy.nan
							else:
								tmp_avg = numpy.average(radial_op[l][s]["avg"][c_size][f_nb][n])
								tmp_std = numpy.std(radial_op[l][s]["avg"][c_size][f_nb][n])
								radial_op[l][s]["avg"][c_size][f_nb][n] = tmp_avg
								radial_op[l][s]["std"][c_size][f_nb][n] = tmp_std
				
	return
def calculate_thickness(f_nb):
	
	#array of associated thickness
	tmp_dist_t2b_dist = MDAnalysis.analysis.distances.distance_array(leaflet_sele["upper"]["all species"].coordinates(), leaflet_sele["lower"]["all species"].coordinates(), U.dimensions)
	tmp_dist_b2t_dist = MDAnalysis.analysis.distances.distance_array(leaflet_sele["lower"]["all species"].coordinates(), leaflet_sele["upper"]["all species"].coordinates(), U.dimensions)	
	tmp_dist_t2b_dist.sort()
	tmp_dist_b2t_dist.sort()
	tmp_dist_t2b_dist = tmp_dist_t2b_dist[:,:args.thick_nb_neighbours]
	tmp_dist_b2t_dist = tmp_dist_b2t_dist[:,:args.thick_nb_neighbours]
	tmp_dist_t2b_avg = numpy.zeros((leaflet_sele["upper"]["all species"].numberOfResidues(),1))
	tmp_dist_b2t_avg = numpy.zeros((leaflet_sele["lower"]["all species"].numberOfResidues(),1))
	tmp_dist_t2b_avg[:,0] = numpy.average(tmp_dist_t2b_dist, axis=1)
	tmp_dist_b2t_avg[:,0] = numpy.average(tmp_dist_b2t_dist, axis=1)
	
	#add frame entry
	for s in leaflet_species["both"] + ["all species"]:
		lipids_thick_nff["data"][s][f_nb] = []
		
	#upper leaflet
	for r_index in range(0,leaflet_sele["upper"]["all species"].numberOfResidues()):
		r_num = leaflet_sele["upper"]["all species"].resnums()[r_index]
		r_specie = leaflet_sele["upper"]["all species"].resnames()[r_index]
		r_index_specie = lipids_resnums_list["upper"][r_specie].index(r_num)		
		#lipids
		lipids_thick_nff["upper"][r_specie][r_index_specie][f_nb] = tmp_dist_t2b_avg[r_index,0]
		#species
		lipids_thick_nff["data"][r_specie][f_nb].append(tmp_dist_t2b_avg[r_index,0])
		lipids_thick_nff["data"][r_specie]["all frames"].append(tmp_dist_t2b_avg[r_index,0])
		lipids_thick_nff["data"]["all species"][f_nb].append(tmp_dist_t2b_avg[r_index,0])
		lipids_thick_nff["data"]["all species"]["all frames"].append(tmp_dist_t2b_avg[r_index,0])
		
	#lower leaflet
	for r_index in range(0,leaflet_sele["lower"]["all species"].numberOfResidues()):
		r_num = leaflet_sele["lower"]["all species"].resnums()[r_index]
		r_specie = leaflet_sele["lower"]["all species"].resnames()[r_index]
		r_index_specie = lipids_resnums_list["lower"][r_specie].index(r_num)
		#lipids
		lipids_thick_nff["lower"][r_specie][r_index_specie][f_nb] = tmp_dist_b2t_avg[r_index,0]
		#species
		lipids_thick_nff["data"][r_specie][f_nb].append(tmp_dist_b2t_avg[r_index,0])
		lipids_thick_nff["data"][r_specie]["all frames"].append(tmp_dist_b2t_avg[r_index,0])
		lipids_thick_nff["data"]["all species"][f_nb].append(tmp_dist_b2t_avg[r_index,0])
		lipids_thick_nff["data"]["all species"]["all frames"].append(tmp_dist_b2t_avg[r_index,0])
		
	return
def calculate_order_parameters(f_nb):
	
	#add frame entry
	for tail in ["tailA", "tailB", "tails"]:
		for l in ["lower", "upper"]:
			lipids_op_nff["data"][tail][l]["all species"][f_nb] = []	
	
	#non flipflopping lipids
	#-----------------------
	for l in ["lower","upper"]:
		for s in op_lipids_handled[l]:
			#retrieve tail boundaries for current lipid type
			tail_A_start = tail_boundaries[s][0]
			tail_B_start = tail_boundaries[s][2]
			tail_A_length = tail_boundaries[s][1]
			tail_B_length = tail_boundaries[s][3]
										
			#calculate 'op' for each bond in lipid tails
			b_index = 0
			tmp_bond_array = numpy.zeros((leaflet_sele[l][s].numberOfResidues(),tail_A_length + tail_B_length))
			for bond in op_bonds[s][tail_A_start:]:
				v = numpy.zeros((leaflet_sele_atoms[l][s].numberOfResidues(),3))
				v_norm2 = numpy.zeros((leaflet_sele_atoms[l][s].numberOfResidues(),1))
				v[:,0] = leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[0])).coordinates()[:,0] - leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[1])).coordinates()[:,0]
				v[:,1] = leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[0])).coordinates()[:,1] - leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[1])).coordinates()[:,1]
				v[:,2] = leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[0])).coordinates()[:,2] - leaflet_sele_atoms[l][s].selectAtoms("name " + str(bond[1])).coordinates()[:,2]
				v_norm2[:,0] = v[:,0]**2 + v[:,1]**2 + v[:,2]**2
				tmp_bond_array[:,b_index] = 0.5*(3*(v[:,2]**2)/v_norm2[:,0]-1)
				b_index += 1

			#calculate order parameters
			tmp_op_tails_array = numpy.zeros((leaflet_sele[l][s].numberOfResidues(),3))  #tail A / tail B / both
			tmp_op_tails_array[:,0] = numpy.average(tmp_bond_array[:,0:tail_A_length],axis=1)
			tmp_op_tails_array[:,1] = numpy.average(tmp_bond_array[:,tail_A_length:tail_A_length+tail_B_length],axis=1)
			tmp_op_tails_array[:,2] = numpy.average(tmp_op_tails_array[:,:2], axis=1)
			
			#store value: each lipid
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				lipids_op_nff[l][s][r_index][f_nb] = tmp_op_tails_array[r_index, 2]

			#store value: species
			tail_names = ["tailA", "tailB", "tails"]
			for tail_index in range(0,3):
				tail = tail_names[tail_index]				
				lipids_op_nff["data"][tail][l][s][f_nb] = list(tmp_op_tails_array[:, tail_index])
				lipids_op_nff["data"][tail][l][s]["all frames"] += list(tmp_op_tails_array[:, tail_index])
				lipids_op_nff["data"][tail][l]["all species"][f_nb] += list(tmp_op_tails_array[:, tail_index])
				lipids_op_nff["data"][tail][l]["all species"]["all frames"] += list(tmp_op_tails_array[:, tail_index])
				
	#flipflopping lipids
	#-------------------
	if args.selection_file_ff!="no":
		for l_index in range(0, lipids_ff_nb):
			#retrieve lipid info
			tmp_s = lipids_ff_info[l_index][0]
			tail_A_start = tail_boundaries[tmp_s][0]
			tail_B_start = tail_boundaries[tmp_s][2]
			tail_A_length = tail_boundaries[tmp_s][1]
			tail_B_length = tail_boundaries[tmp_s][3]

			#calculate order parameters
			tmp_bond_array = []
			for bond in op_bonds[tmp_s][tail_A_start:]:
				vx = lipids_sele_ff[l_index].selectAtoms("name " + str(bond[0])).coordinates()[0,0] - lipids_sele_ff[l_index].selectAtoms("name " + str(bond[1])).coordinates()[0,0]
				vy = lipids_sele_ff[l_index].selectAtoms("name " + str(bond[0])).coordinates()[0,1] - lipids_sele_ff[l_index].selectAtoms("name " + str(bond[1])).coordinates()[0,1]
				vz = lipids_sele_ff[l_index].selectAtoms("name " + str(bond[0])).coordinates()[0,2] - lipids_sele_ff[l_index].selectAtoms("name " + str(bond[1])).coordinates()[0,2]
				v_norm2 = vx**2 + vy**2 + vz**2
				tmp_bond_array.append(0.5*(3*(vz**2)/float(v_norm2)-1))
			lipids_op_ff_tailA[l_index][f_nb] = numpy.average(tmp_bond_array[0:tail_A_length])
			lipids_op_ff_tailB[l_index][f_nb] = numpy.average(tmp_bond_array[tail_A_length:tail_A_length+tail_B_length])
			lipids_op_ff_tails[l_index][f_nb] = numpy.average(lipids_op_ff_tailA[l_index][f_nb], lipids_op_ff_tailB[l_index][f_nb])
			if args.xtcfilename == "no":
				lipids_op_ff_tailA[l_index]["all frames"] = numpy.average(tmp_bond_array[0:tail_A_length])
				lipids_op_ff_tailB[l_index]["all frames"] = numpy.average(tmp_bond_array[tail_A_length:tail_A_length+tail_B_length])
				lipids_op_ff_tails[l_index]["all frames"]= numpy.average(lipids_op_ff_tailA[l_index][f_nb], lipids_op_ff_tailB[l_index][f_nb])
			else:
				lipids_op_ff_tailA[l_index]["all frames"].append(numpy.average(tmp_bond_array[0:tail_A_length]))
				lipids_op_ff_tailB[l_index]["all frames"].append(numpy.average(tmp_bond_array[tail_A_length:tail_A_length+tail_B_length]))
				lipids_op_ff_tails[l_index]["all frames"].append(numpy.average(lipids_op_ff_tailA[l_index][f_nb], lipids_op_ff_tailB[l_index][f_nb]))

	return
def detect_clusters_connectivity(dist, box_dim):
	
	#use networkx algorithm
	connected=(dist<args.cutoff_connect)
	network=nx.Graph(connected)
	groups=nx.connected_components(network)
	
	return groups
def detect_clusters_density(dist, box_dim):
	
	#run DBSCAN algorithm
	dbscan_output=DBSCAN(eps=args.dbscan_dist,metric='precomputed',min_samples=args.dbscan_nb).fit(dist)

	#build 'groups' structure i.e. a list whose element are all the clusters identified
	groups=[]
	for c_lab in numpy.unique(dbscan_output.labels_):
		tmp_pos=numpy.argwhere(dbscan_output.labels_==c_lab)
		if c_lab==-1:
			for p in tmp_pos:
				groups.append([p[0]])
		else:
			g=[]
			for p in tmp_pos:
				g.append(p[0])
			groups.append(g)

	return groups
def rolling_avg(loc_list):
	
	loc_arr=numpy.asarray(loc_list)
	shape=(loc_arr.shape[-1]-args.nb_smoothing+1,args.nb_smoothing)
	strides=(loc_arr.strides[-1],loc_arr.strides[-1])   	
	return numpy.average(numpy.lib.stride_tricks.as_strided(loc_arr, shape=shape, strides=strides), -1)
def update_colours_sizes():
	global colours_sizes
	colours_sizes = {}
	tmp_cmap = cm.get_cmap('jet')
	colours_sizes_value = tmp_cmap(numpy.linspace(0, 1, len(radial_sizes["all frames"])))
	for c_index in range(0, len(radial_sizes["all frames"])):
		c_size = radial_sizes["all frames"][c_index]
		colours_sizes[c_size] = colours_sizes_value[c_index]

	return
def sort_and_smooth_data():
	
	global time_sorted
	global z_lower_sorted
	global z_upper_sorted
	global z_ff_sorted	
	global time_smooth
	global z_lower_smooth
	global z_upper_smooth
	global z_ff_smooth
	time_sorted = []
	z_lower_sorted = []
	z_upper_sorted = []
	z_ff_sorted = []
	time_smooth = []
	z_lower_smooth = []
	z_upper_smooth = []
	z_ff_smooth = []

	#sort data
	#=========
	#time and bilayer coords
	#-----------------------
	for f_nb in sorted(time_stamp.keys()):
		time_sorted.append(time_stamp[f_nb])
		z_upper_sorted.append(z_upper[f_nb])
		z_lower_sorted.append(z_lower[f_nb])

	#thickness
	#---------
	if args.perturb == 1 or args.perturb == 3:
		for f_nb in sorted(time_stamp.keys()):
			for s in leaflet_species["both"] + ["all species"]:
				lipids_thick_nff["sorted"]["avg"][s].append(numpy.average(lipids_thick_nff["data"][s][f_nb]))
				lipids_thick_nff["sorted"]["std"][s].append(numpy.std(lipids_thick_nff["data"][s][f_nb]))

	#order parameter
	#---------------
	if args.perturb == 2 or args.perturb == 3:
		#nff lipids
		for f_nb in sorted(time_stamp.keys()):
			for l in ["lower","upper"]:
				for s in op_lipids_handled[l]:
					for tail in ["tailA","tailB","tails"]:
						lipids_op_nff["sorted"]["avg"][tail][l][s].append(numpy.average(lipids_op_nff["data"][tail][l][s][f_nb]))
						lipids_op_nff["sorted"]["std"][tail][l][s].append(numpy.std(lipids_op_nff["data"][tail][l][s][f_nb]))
					
		#ff lipids
		if args.selection_file_ff != "no":
			for l in range(0,lipids_ff_nb):
				for f_nb in sorted(time_stamp.keys()):		
					lipids_op_ff_tailA_sorted[l].append(lipids_op_ff_tailA[l][f_nb])
					lipids_op_ff_tailB_sorted[l].append(lipids_op_ff_tailB[l][f_nb])
					lipids_op_ff_tails_sorted[l].append(lipids_op_ff_tails[l][f_nb])
				lipids_op_ff_tailA[l]["all frames"] = numpy.average(lipids_op_ff_tailA[l]["all frames"])
				lipids_op_ff_tailB[l]["all frames"] = numpy.average(lipids_op_ff_tailB[l]["all frames"])
				lipids_op_ff_tails[l]["all frames"] = numpy.average(lipids_op_ff_tails[l]["all frames"])

	#smooth data
	#===========
	if args.nb_smoothing > 1:

		#time and bilayer coords
		#-----------------------
		time_smooth = rolling_avg(time_sorted)
		z_upper_smooth = rolling_avg(z_upper_sorted)
		z_lower_smooth = rolling_avg(z_lower_sorted)	

		#thickness
		#---------
		if args.perturb == 1 or args.perturb == 3:
			for s in leaflet_species["both"] + ["all species"]:
				lipids_thick_nff["smoothed"]["avg"][s] = rolling_avg(lipids_thick_nff["sorted"]["avg"][s])
				lipids_thick_nff["smoothed"]["std"][s] = rolling_avg(lipids_thick_nff["sorted"]["std"][s])
		
		#order parameter
		#---------------
		if args.perturb == 2 or args.perturb == 3:
			for l in ["lower","upper"]:
				for s in op_lipids_handled[l]:
					for tail in ["tailA","tailB","tails"]:
						lipids_op_nff["smoothed"]["avg"][tail][l][s] = rolling_avg(lipids_op_nff["sorted"]["avg"][tail][l][s])
						lipids_op_nff["smoothed"]["std"][tail][l][s] = rolling_avg(lipids_op_nff["sorted"]["avg"][tail][l][s])

			#ff lipids
			if args.selection_file_ff!= "no":
				for f_nb in sorted(time_stamp.keys()):		
					for l in range(0,lipids_ff_nb):
						z_ff_sorted[l].append(z_ff[l][f_nb])
						lipids_op_ff_tailA_sorted[l].append(lipids_op_ff_tailA[l][f_nb])
						lipids_op_ff_tailB_sorted[l].append(lipids_op_ff_tailB[l][f_nb])
						lipids_op_ff_tails_sorted[l].append(lipids_op_ff_tails[l][f_nb])
				for l in range(0,lipids_ff_nb):
					z_ff_smooth[l] = rolling_avg(z_ff_sorted[l])
					lipids_op_ff_tailA_smooth[l] = rolling_avg(lipids_op_ff_tailA_sorted[l])
					lipids_op_ff_tailB_smooth[l] = rolling_avg(lipids_op_ff_tailB_sorted[l])
					lipids_op_ff_tails_smooth[l] = rolling_avg(lipids_op_ff_tails_sorted[l])

	return

################################################################################################################################################
# FUNCTIONS: outputs
################################################################################################################################################

#thickness
def thick_xvg_write():
	
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/xvg/1_2_thickness_species.txt'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/xvg/1_2_thickness_species.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by thickness v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_thickness_species.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of bilayer thickness by lipid specie\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"thickness\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(leaflet_species["both"])) + "\n")
	for s_index in range(0,len(leaflet_species["both"])):
		s = leaflet_species["both"][s_index]
		output_xvg.write("@ s" + str(s_index) + " legend \"" + str(s) + " (avg)\"\n")
		output_txt.write("1_2_thickness_species.xvg," + str(s_index+1) + "," + str(s) + " (avg)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
	for s_index in range(0,len(leaflet_species["both"])):
		s=leaflet_species["both"][s_index]
		output_xvg.write("@ s" + str(len(leaflet_species["both"])+s_index) + " legend \"" + str(s) + " (std)\"\n")
		output_txt.write("1_2_thickness_species.xvg," + str(len(leaflet_species["both"])+s_index+1) + "," + str(s) + " (std)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
	output_txt.close()
	for f_nb in sorted(time_stamp.iterkeys()):
		results = str(time_stamp[f_nb])
		for s in leaflet_species["both"]:
			results += "	" + str(round(numpy.average(lipids_thick_nff["data"][s][f_nb]),2))
		for s in leaflet_species["both"]:
			results += "	" + str(round(numpy.std(lipids_thick_nff["data"][s][f_nb]),2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def thick_xvg_write_smoothed():
	
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/smoothed/xvg/1_4_thickness_species_smoothed.txt'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/smoothed/xvg/1_4_thickness_species_smoothed.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by thickness v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_4_thickness_species_smoothed.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of bilayer thickness by lipid specie\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"thickness\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(leaflet_species["both"])) + "\n")
	for s_index in range(0,len(leaflet_species["both"])):
		s = leaflet_species["both"][s_index]
		output_xvg.write("@ s" + str(s_index) + " legend \"" + str(s) + " (avg)\"\n")
		output_txt.write("1_4_thickness_species_smoothed.xvg," + str(s_index+1) + "," + str(s) + " (avg)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
	for s_index in range(0,len(leaflet_species["both"])):
		s=leaflet_species["both"][s_index]
		output_xvg.write("@ s" + str(len(leaflet_species["both"]) + s_index) + " legend \"" + str(s) + " (std)\"\n")
		output_txt.write("1_4_thickness_species_smoothed.xvg," + str(len(leaflet_species["both"]) + s_index + 1) + "," + str(s) + " (std)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
	output_txt.close()
	for t_index in range(0, len(time_smoothed)):
		results = str(time_smoothed[t_index])
		for s in leaflet_species["both"]:
			results += "	" + str(round(lipids_thick_nff["smoothed"]["avg"][s][t_index],2))
		for s in leaflet_species["both"]:
			results += "	" + str(round(lipids_thick_nff["smoothed"]["std"][s][t_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def thick_xvg_graph():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/png/1_1_thickness.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/1_1_thickness.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 4))
	fig.suptitle("Evolution of lipid bilayer thickness")
					
	#plot data:
	#----------
	ax1 = fig.add_subplot(111)
	p_upper={}
	for s in leaflet_species["both"]:
		p_upper[s]=plt.plot(time_sorted, lipids_thick_nff["sorted"]["avg"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		p_upper[str(s + "_err")]=plt.fill_between(time_sorted, numpy.asarray(lipids_thick_nff["sorted"]["avg"][s])-numpy.asarray(lipids_thick_nff["sorted"]["std"][s]), numpy.asarray(lipids_thick_nff["sorted"]["avg"][s])+numpy.asarray(lipids_thick_nff["sorted"]["std"][s]), color=colours_lipids[s], alpha=0.2)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('thickness ($\AA$)', fontsize="small")
	
	#save figure
	#-----------
	ax1.set_ylim(numpy.average(lipids_thick_nff["sorted"]["avg"]["all species"]) - numpy.average(lipids_thick_nff["sorted"]["std"]["all species"]) - 5, numpy.average(lipids_thick_nff["sorted"]["avg"]["all species"]) + numpy.average(lipids_thick_nff["sorted"]["std"]["all species"]) + 5)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def thick_xvg_graph_smoothed():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/smoothed/png/1_3_thickness.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/smoothed/1_3_thickness.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 4))
	fig.suptitle("Evolution of lipid bilayer thickness")
					
	#plot data:
	#----------
	ax1 = fig.add_subplot(111)
	p_upper={}
	for s in leaflet_species["both"]:
		p_upper[s]=plt.plot(time_smoothed, lipids_thick_nff["smoothed"]["avg"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		p_upper[str(s + "_err")]=plt.fill_between(time_smoothed, numpy.asarray(lipids_thick_nff["smoothed"]["avg"][s])-numpy.asarray(lipids_thick_nff["smoothed"]["std"][s]), numpy.asarray(lipids_thick_nff["smoothed"]["avg"][s])+numpy.asarray(lipids_thick_nff["smoothed"]["std"][s]), color=colours_lipids[s], alpha=0.2)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('thickness ($\AA$)', fontsize="small")
	
	#save figure
	#-----------
	ax1.set_ylim(numpy.average(lipids_thick_nff["smoothed"]["avg"]["all species"]) - numpy.average(lipids_thick_nff["smoothed"]["std"]["all species"]) - 5, numpy.average(lipids_thick_nff["smoothed"]["avg"]["all species"]) + numpy.average(lipids_thick_nff["smoothed"]["std"]["all species"]) + 5)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def thick_frame_write_stat(f_nb, f_time):

	#create file
	if f_nb == "all frames":
		filename_details = os.getcwd() + '/' + str(args.output_folder) + '/thickness/1_species/1_0_thickness.stat'
	else:
		filename_details = os.getcwd() + '/' + str(args.output_folder) + '/thickness/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_thickness_' + str(int(f_time)).zfill(5) + 'ns.stat'
	output_stat = open(filename_details, 'w')		
	output_stat.write("[thickness statistics - written by thickness v" + str(version_nb) + "]\n")
	output_stat.write("\n")
	
	#general info
	output_stat.write("1. membrane composition: \n")
	output_stat.write(membrane_comp["upper"][:-1] + "\n")
	output_stat.write(membrane_comp["lower"][:-1] + "\n")	
	if args.xtcfilename != "no":
		output_stat.write("\n")
		output_stat.write("2. nb frames processed:	" + str(nb_frames_processed) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if f_nb != "all frames":		
		output_stat.write("\n")
		output_stat.write("3. time: " + str(f_time) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")		
	output_stat.write("\n")
	
	#results data
	output_stat.write("Bilayer thickness:\n")
	output_stat.write("-----------------\n")
	output_stat.write("avg=" + str(round(numpy.average(lipids_thick_nff["data"]["all species"][f_nb]),2)) + "\n")
	output_stat.write("std=" + str(round(numpy.std(lipids_thick_nff["data"]["all species"][f_nb]),2)) + "\n")
	output_stat.write("max=" + str(round(numpy.max(lipids_thick_nff["data"]["all species"][f_nb]),2)) + "\n")
	output_stat.write("min=" + str(round(numpy.min(lipids_thick_nff["data"]["all species"][f_nb]),2)) + "\n")
	output_stat.write("\n")
	output_stat.write("Average bilayer thickness for each specie:\n")
	output_stat.write("------------------------------------------\n")
	for s in leaflet_species["both"]:
		output_stat.write(str(s) + "	" + str(round(numpy.average(lipids_thick_nff["data"][s][f_nb]),2)) + "	(" + str(round(numpy.std(lipids_thick_nff["data"][s][f_nb]),2)) + ")\n")
	output_stat.write("\n")
	output_stat.close()		
	
	return
def thick_frame_write_snapshot(f_nb, f_time):

	#store order parameter info in beta factor field
	for l in ["lower","upper"]:
		for s in leaflet_species[l]:
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				lipids_sele_nff[l][s][r_index].set_bfactor(lipids_thick_nff[l][s][r_index][f_nb])

	#case: gro file
	if args.xtcfilename=="no":
		all_atoms.write(os.getcwd() + '/' + str(args.output_folder) + '/thickness/2_snapshots/' + args.grofilename[:-4] + '_annotated_thickness', format="PDB")

	#case: xtc file
	else:
		tmp_name=os.getcwd() + "/" + str(args.output_folder) + '/thickness/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_thickness_' + str(int(f_time)).zfill(5) + 'ns.pdb'
		W=Writer(tmp_name, nb_atoms)
		W.write(all_atoms)
	
	return
def thick_frame_write_annotation(f_nb, f_time):
	
	#create file
	if args.xtcfilename == "no":
		filename_details=os.getcwd() + "/" + str(args.output_folder) + '/thickness/2_snapshots/' + args.grofilename[:-4] + '_annotated_thickness.txt'
	else:
		filename_details = os.getcwd() + "/" + str(args.output_folder) + '/thickness/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_thickness_' + str(int(f_time)).zfill(5) + 'ns.txt'
	output_stat = open(filename_details, 'w')		

	#create selection string
	tmp_sele_string = ""
	for l in ["lower","upper"]:
		for s in leaflet_species[l]:
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				tmp_sele_string += "." + lipids_sele_nff_VMD_string[l][s][r_index]
	output_stat.write(tmp_sele_string[1:] + "\n")

	#write min and max boundaries of thickness
	output_stat.write(str(round(numpy.min(lipids_thick_nff["data"]["all species"][f_nb]),2)) + ";" + str(round(numpy.max(lipids_thick_nff["data"]["all species"][f_nb]),2)) + "\n")
	
	#ouptut thickness for each lipid
	tmp_thick = "1"
	for l in ["lower","upper"]:
		for s in leaflet_species[l]:
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				tmp_thick += ";" + str(round(lipids_thick_nff[l][s][r_index][f_nb],2))
	output_stat.write(tmp_thick + "\n")			
	output_stat.close()

	return
def thick_xtc_write_annotation():
	
	#create file
	filename_details=os.getcwd() + '/' + str(args.output_folder) + '/thickness/3_VMD/' + args.xtcfilename[:-4] + '_annotated_thickness_dt' + str(args.frames_dt) + '.txt'
	output_stat = open(filename_details, 'w')		

	#create selection string
	tmp_sele_string = ""
	for l in ["lower","upper"]:
		for s in leaflet_species[l]:
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				tmp_sele_string += "." + lipids_sele_nff_VMD_string[l][s][r_index]
	output_stat.write(tmp_sele_string[1:] + "\n")

	#write min and max boundaries of thickness
	output_stat.write(str(round(numpy.min(lipids_thick_nff["data"]["all species"]["all frames"]),2)) + ";" + str(round(numpy.max(lipids_thick_nff["data"]["all species"]["all frames"]),2)) + "\n")
	
	#ouptut thickness for each lipid
	for f_nb in sorted(time_stamp.iterkeys()):
		tmp_thick = str(f_nb)
		for l in ["lower","upper"]:
			for s in leaflet_species[l]:
				for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
					tmp_thick += ";" + str(round(lipids_thick_nff[l][s][r_index][f_nb],2))
		output_stat.write(tmp_thick + "\n")
			
	output_stat.close()

	return

#order parameters
def op_xvg_ff_write():
	
	#flipflops: upper to lower
	if numpy.size(lipids_ff_u2l_index)>0:
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/xvg/4_3_order_param_ff_u2l.txt'
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/xvg/4_3_order_param_ff_u2l.xvg'
		output_txt  =  open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 4_3_order_param_ff_u2l.xvg.\n")
		output_xvg  =  open(filename_xvg, 'w')
		output_xvg.write("@ title \"Evolution of the tail order parameters of flipflopping lipids\"\n")
		output_xvg.write("@ xaxis  label \"time (ns)\"\n")
		output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(lipids_ff_nb*3) + "\n")
		for l_index in range(0,len(lipids_ff_u2l_index)):
			l = lipids_ff_u2l_index[l_index]
			output_xvg.write("@ s" + str(3*l_index) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A\"\n")
			output_xvg.write("@ s" + str(3*l_index+1) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B\"\n")
			output_xvg.write("@ s" + str(3*l_index+2) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both\"\n")
			output_txt.write("4_3_order_param_ff_u2l.xvg," + str((3*l_index)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A,auto\n")
			output_txt.write("4_3_order_param_ff_u2l.xvg," + str((3*l_index+1)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B,auto\n")
			output_txt.write("4_3_order_param_ff_u2l.xvg," + str((3*l_index+2)+1) +"," + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both,auto\n")
		output_txt.close()
		for f_nb in sorted(time_stamp.iterkeys()):
			results = str(time_stamp[f_nb])
			for l in lipids_ff_u2l_index:
				results += "	" + str(round(lipids_op_ff_tailA[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tailB[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tails[l][f_nb],2))
			output_xvg.write(results + "\n")
		output_xvg.close()
	
	#flipflops: lower to upper
	if numpy.size(lipids_ff_l2u_index)>0:
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/xvg/4_3_order_param_ff_l2u.txt'
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/xvg/4_3_order_param_ff_l2u.xvg'
		output_txt  =  open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 4_3_order_param_ff_l2u.xvg.\n")
		output_xvg  =  open(filename_xvg, 'w')
		output_xvg.write("@ title \"Evolution of the tail order parameters of flipflopping lipids\"\n")
		output_xvg.write("@ xaxis  label \"time (ns)\"\n")
		output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(lipids_ff_nb*3) + "\n")
		for l_index in range(0,len(lipids_ff_l2u_index)):
			l = lipids_ff_l2u_index[l_index]
			output_xvg.write("@ s" + str(3*l_index) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A\"\n")
			output_xvg.write("@ s" + str(3*l_index+1) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B\"\n")
			output_xvg.write("@ s" + str(3*l_index+2) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both\"\n")
			output_txt.write("4_3_order_param_ff_l2u.xvg," + str((3*l_index)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A,auto\n")
			output_txt.write("4_3_order_param_ff_l2u.xvg," + str((3*l_index+1)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B,auto\n")
			output_txt.write("4_3_order_param_ff_l2u.xvg," + str((3*l_index+2)+1) +"," + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both,auto\n")
		output_txt.close()
		for f_nb in sorted(time_stamp.iterkeys()):
			results = str(time_stamp[f_nb])
			for l in lipids_ff_l2u_index:
				results+= "	" + str(round(lipids_op_ff_tailA[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tailB[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tails[l][f_nb],2))
			output_xvg.write(results + "\n")
		output_xvg.close()

	return
def op_xvg_ff_write_smoothed():

	#flipflops: upper to lower
	if numpy.size(lipids_ff_u2l_index)>0:
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/xvg/4_5_order_param_ff_u2l_smoothed.txt'
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/xvg/4_5_order_param_ff_u2l_smoothed.xvg'
		output_txt  =  open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 4_5_order_param_ff_u2l_smoothed.xvg.\n")
		output_xvg  =  open(filename_xvg, 'w')
		output_xvg.write("@ title \"Evolution of the tail order parameters of flipflopping lipids\"\n")
		output_xvg.write("@ xaxis  label \"time (ns)\"\n")
		output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(lipids_ff_nb*3) + "\n")
		for l_index in range(0,len(lipids_ff_u2l_index)):
			l = lipids_ff_u2l_index[l_index]
			output_xvg.write("@ s" + str(3*l_index) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A\"\n")
			output_xvg.write("@ s" + str(3*l_index+1) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B\"\n")
			output_xvg.write("@ s" + str(3*l_index+2) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both\"\n")
			output_txt.write("4_3_order_param_ff_u2l_smoothed.xvg," + str((3*l_index)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A,auto\n")
			output_txt.write("4_3_order_param_ff_u2l_smoothed.xvg," + str((3*l_index+1)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B,auto\n")
			output_txt.write("4_3_order_param_ff_u2l_smoothed.xvg," + str((3*l_index+2)+1) +"," + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both,auto\n")
		output_txt.close()
		for t_index in range(0, numpy.size(time_smooth)):
			results = str(time_smooth[t_index])
			for l in lipids_ff_u2l_index:
				results+= "	" + str(round(lipids_op_ff_tailA_smooth[l][t_index],2)) + "	" + str(round(lipids_op_ff_tailB_smooth[l][t_index],2)) + "	" + str(round(lipids_op_ff_tails_smooth[l][t_index],2))
			output_xvg.write(results + "\n")
		output_xvg.close()
	
	#flipflops: lower to upper
	if numpy.size(lipids_ff_l2u_index)>0:
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/xvg/4_5_order_param_ff_l2u_smoothed.txt'
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/xvg/4_5_order_param_ff_l2u_smoothed.xvg'
		output_txt  =  open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 4_5_order_param_ff_l2u_smoothed.xvg.\n")
		output_xvg  =  open(filename_xvg, 'w')
		output_xvg.write("@ title \"Evolution of the tail order parameters of flipflopping lipids\"\n")
		output_xvg.write("@ xaxis  label \"time (ns)\"\n")
		output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(lipids_ff_nb*3) + "\n")
		for l_index in range(0,len(lipids_ff_l2u_index)):
			l = lipids_ff_l2u_index[l_index]
			output_xvg.write("@ s" + str(3*l_index) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A\"\n")
			output_xvg.write("@ s" + str(3*l_index+1) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B\"\n")
			output_xvg.write("@ s" + str(3*l_index+2) + " legend \"" + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both\"\n")
			output_txt.write("4_3_order_param_ff_l2u_smoothed.xvg," + str((3*l_index)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail A,auto\n")
			output_txt.write("4_3_order_param_ff_l2u_smoothed.xvg," + str((3*l_index+1)+1) + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " tail B,auto\n")
			output_txt.write("4_3_order_param_ff_l2u_smoothed.xvg," + str((3*l_index+2)+1) +"," + "," + str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]) + " both,auto\n")
		output_txt.close()
		for t_index in range(0, numpy.size(time_smooth)):
			results = str(time_smooth[t_index])
			for l in lipids_ff_l2u_index:
				results+= "	" + str(round(lipids_op_ff_tailA_smooth[l][t_index],2)) + "	" + str(round(lipids_op_ff_tailB_smooth[l][t_index],2)) + "	" + str(round(lipids_op_ff_tails_smooth[l][t_index],2))
			output_xvg.write(results + "\n")
		output_xvg.close()

	return
def op_xvg_ff_graph():

	#upper to lower flipflops
	#========================
	if numpy.size(lipids_ff_u2l_index)>0:

		#create filenames
		#----------------
		filename_png=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/png/4_2_order_param_ff_u2l.png'
		filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/4_2_order_param_ff_u2l.svg'
		
		#create figure
		#-------------
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("Flipflopping lipids: upper to lower")
	
		#create data
		#-----------
		tmp_time=[]
		tmp_z_ff={}
		tmp_z_upper=[]
		tmp_z_lower=[]
		tmp_lipids_op_ff_tails={}
		for l in lipids_ff_u2l_index:
			tmp_z_ff[l]=[]
			tmp_lipids_op_ff_tails[l]=[]
		for f_nb in sorted(time_stamp.iterkeys()):
			tmp_time.append(time_stamp[f_nb])
			tmp_z_upper.append(z_upper[f_nb])
			tmp_z_lower.append(z_lower[f_nb])
			for l in lipids_ff_u2l_index:
				tmp_z_ff[l].append(z_ff[l][f_nb])
				tmp_lipids_op_ff_tails[l].append(lipids_op_ff_tails[l][f_nb])
	
		#plot data: order paramter
		#-------------------------
		ax1 = fig.add_subplot(211)
		p_upper={}
		for l in lipids_ff_u2l_index:
			p_upper[l]=plt.plot(tmp_time, tmp_lipids_op_ff_tails[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
	
		#plot data: z coordinate
		#-----------------------
		ax2 = fig.add_subplot(212)
		p_lower={}
		p_lower["upper"]=plt.plot(tmp_time, tmp_z_upper, linestyle='dashed', color='k')
		p_lower["lower"]=plt.plot(tmp_time, tmp_z_lower, linestyle='dashed', color='k')
		for l in lipids_ff_u2l_index:
			p_lower[l]=plt.plot(tmp_time, tmp_z_ff[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax2.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('z coordinate', fontsize="small")
		
		#save figure
		#-----------
		ax1.set_ylim(-0.5, 1)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		ax2.set_ylim(-40, 40)
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=3))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	#lower to upper flipflops
	#========================
	if numpy.size(lipids_ff_l2u_index)>0:

		#create filenames
		#----------------
		filename_png=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/png/4_2_order_param_ff_l2u.png'
		filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/4_2_order_param_ff_l2u.svg'
		
		#create figure
		#-------------
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("Flipflopping lipids: upper to lower")
	
		#create data
		#-----------
		tmp_time=[]
		tmp_z_ff={}
		tmp_z_upper=[]
		tmp_z_lower=[]
		tmp_lipids_op_ff_tails={}
		for l in lipids_ff_l2u_index:
			tmp_z_ff[l]=[]
			tmp_lipids_op_ff_tails[l]=[]
		for f_nb in sorted(time_stamp.iterkeys()):
			tmp_time.append(time_stamp[f_nb])
			tmp_z_upper.append(z_upper[f_nb])
			tmp_z_lower.append(z_lower[f_nb])
			for l in lipids_ff_l2u_index:
				tmp_z_ff[l].append(z_ff[l][f_nb])
				tmp_lipids_op_ff_tails[l].append(lipids_op_ff_tails[l][f_nb])
	
		#plot data: order paramter
		#-------------------------
		ax1 = fig.add_subplot(211)
		p_upper={}
		for l in lipids_ff_l2u_index:
			p_upper[l]=plt.plot(tmp_time, tmp_lipids_op_ff_tails[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
	
		#plot data: z coordinate
		#-----------------------
		ax2 = fig.add_subplot(212)
		p_lower={}
		p_lower["upper"]=plt.plot(tmp_time, tmp_z_upper, linestyle='dashed', color='k')
		p_lower["lower"]=plt.plot(tmp_time, tmp_z_lower, linestyle='dashed', color='k')
		for l in lipids_ff_l2u_index:
			p_lower[l]=plt.plot(tmp_time, tmp_z_ff[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('z coordinate', fontsize="small")
		
		#save figure
		#-----------
		ax1.set_ylim(-0.5, 1)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		ax2.set_ylim(-40, 40)
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=3))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	return
def op_xvg_ff_graph_smoothed():
	
	#upper to lower flipflops
	#========================
	if numpy.size(lipids_ff_u2l_index)>0:

		#create filenames
		#----------------
		filename_png=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/png/4_4_order_param_ff_u2l_smoothed.png'
		filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/4_4_order_param_ff_u2l_smoothed.svg'
		
		#create figure
		#-------------
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("Flipflopping lipids: upper to lower")
	
		#create data
		#-----------
		tmp_time=[]
		tmp_z_ff={}
		tmp_z_upper=[]
		tmp_z_lower=[]
		tmp_lipids_op_ff_tails={}
		for l in lipids_ff_u2l_index:
			tmp_z_ff[l]=[]
			tmp_lipids_op_ff_tails[l]=[]
		for t_index in range(0, numpy.size(time_smooth)):
			tmp_time.append(time_smooth[t_index])
			tmp_z_upper.append(z_upper_smooth[t_index])
			tmp_z_lower.append(z_lower_smooth[t_index])
			for l in lipids_ff_u2l_index:
				tmp_z_ff[l].append(z_ff_smooth[l][t_index])
				tmp_lipids_op_ff_tails[l].append(lipids_op_ff_tails_smooth[l][t_index])
		
		#plot data: order paramter
		#-------------------------
		ax1 = fig.add_subplot(211)
		p_upper={}
		for l in lipids_ff_u2l_index:
			p_upper[l]=plt.plot(tmp_time, tmp_lipids_op_ff_tails[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
	
		#plot data: z coordinate
		#-----------------------
		ax2 = fig.add_subplot(212)
		p_lower={}
		p_lower["upper"]=plt.plot(tmp_time, tmp_z_upper, linestyle='dashed', color='k')
		p_lower["lower"]=plt.plot(tmp_time, tmp_z_lower, linestyle='dashed', color='k')
		for l in lipids_ff_u2l_index:
			p_lower[l]=plt.plot(tmp_time, tmp_z_ff[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax2.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('z coordinate', fontsize="small")
		
		#save figure
		#-----------
		ax1.set_ylim(-0.5, 1)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		ax2.set_ylim(-40, 40)
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=3))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	#lower to upper flipflops
	#========================
	if numpy.size(lipids_ff_l2u_index)>0:

		#create filenames
		#----------------
		filename_png=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/png/4_4_order_param_ff_l2u_smoothed.png'
		filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/4_ff/smoothed/4_4_order_param_ff_l2u_smoothed.svg'
		
		#create figure
		#-------------
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("Flipflopping lipids: lower to lower")
		
		#create data
		#-----------
		tmp_time=[]
		tmp_z_ff={}
		tmp_z_upper=[]
		tmp_z_lower=[]
		tmp_lipids_op_ff_tails={}
		for l in lipids_ff_l2u_index:
			tmp_z_ff[l]=[]
			tmp_lipids_op_ff_tails[l]=[]
		for t_index in range(0, numpy.size(time_smooth)):
			tmp_time.append(time_smooth[t_index])
			tmp_z_upper.append(z_upper_smooth[t_index])
			tmp_z_lower.append(z_lower_smooth[t_index])
			for l in lipids_ff_l2u_index:
				tmp_z_ff[l].append(z_ff_smooth[l][t_index])
				tmp_lipids_op_ff_tails[l].append(lipids_op_ff_tails_smooth[l][t_index])

		#plot data: order paramter
		#-------------------------
		ax1 = fig.add_subplot(211)
		p_upper={}
		for l in lipids_ff_l2u_index:
			p_upper[l]=plt.plot(tmp_time, tmp_lipids_op_ff_tails[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
	
		#plot data: z coordinate
		#-----------------------
		ax2 = fig.add_subplot(212)
		p_lower={}
		p_lower["upper"]=plt.plot(tmp_time, tmp_z_upper, linestyle='dashed', color='k')
		p_lower["lower"]=plt.plot(tmp_time, tmp_z_lower, linestyle='dashed', color='k')
		for l in lipids_ff_l2u_index:
			p_lower[l]=plt.plot(tmp_time, tmp_z_ff[l], label=str(lipids_ff_info[l][0]) + " " + str(lipids_ff_info[l][1]))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('time (ns)', fontsize="small")
		plt.ylabel('z coordinate', fontsize="small")
		
		#save figure
		#-----------
		ax1.set_ylim(-0.5, 1)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		ax2.set_ylim(-40, 40)
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=3))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	return
def op_xvg_nff_write():
	
	#lipids in upper leaflet
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/xvg/1_3_order_param_nff_upper.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/xvg/1_3_order_param_nff_upper.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_3_order_param_nff_upper.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of lipid tails order parameters in upper leaflet\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(2*len(op_lipids_handled["upper"])*3) + "\n")
	for s_index in range(0,len(op_lipids_handled["upper"])):
		output_xvg.write("@ s" + str(3*s_index) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail A (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+1) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail B (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+2) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " both (avg)\"\n")
		output_txt.write("1_3_order_param_nff_upper.xvg," + str((3*s_index)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail A (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper.xvg," + str((3*s_index+1)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail B (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper.xvg," + str((3*s_index+2)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " both (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
	for s_index in range(0,len(op_lipids_handled["upper"])):
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["upper"])+3*s_index) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail A (std)\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["upper"])+3*s_index+1) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail (std)B\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["upper"])+3*s_index+2) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " both (std)\"\n")
		output_txt.write("1_3_order_param_nff_upper.xvg," + str(3*len(op_lipids_handled["upper"])+(3*s_index)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail A (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper.xvg," + str(3*len(op_lipids_handled["upper"])+(3*s_index+1)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail B (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper.xvg," + str(3*len(op_lipids_handled["upper"])+(3*s_index+2)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " both (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
	output_txt.close()
	for f_nb in sorted(time_stamp.iterkeys()):
		results=str(time_stamp[f_nb])
		for s in op_lipids_handled["upper"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(numpy.average(lipids_op_nff["data"][tail]["upper"][s][f_nb]),2))
		for s in op_lipids_handled["upper"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(numpy.std(lipids_op_nff["data"][tail]["upper"][s][f_nb]),2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	#lipids in lower leaflet
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/xvg/1_3_order_param_nff_lower.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/xvg/1_3_order_param_nff_lower.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_3_order_param_nff_lower.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of lipid tails order parameters in lower leaflet\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(2*len(op_lipids_handled["lower"])*3) + "\n")
	for s_index in range(0,len(op_lipids_handled["lower"])):
		output_xvg.write("@ s" + str(3*s_index) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail A (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+1) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail B (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+2) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " both (avg)\"\n")
		output_txt.write("1_3_order_param_nff_lower.xvg," + str((3*s_index)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail A (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower.xvg," + str((3*s_index+1)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail B (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower.xvg," + str((3*s_index+2)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " both (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
	for s_index in range(0,len(op_lipids_handled["lower"])):
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["lower"])+3*s_index) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail A (std)\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["lower"])+3*s_index+1) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail (std)B\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["lower"])+3*s_index+2) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " both (std)\"\n")
		output_txt.write("1_3_order_param_nff_lower.xvg," + str(3*len(op_lipids_handled["lower"])+(3*s_index)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail A (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower.xvg," + str(3*len(op_lipids_handled["lower"])+(3*s_index+1)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail B (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower.xvg," + str(3*len(op_lipids_handled["lower"])+(3*s_index+2)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " both (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
	output_txt.close()
	for f_nb in sorted(time_stamp.iterkeys()):
		results=str(time_stamp[f_nb])
		for s in op_lipids_handled["lower"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(numpy.average(lipids_op_nff["data"][tail]["lower"][s][f_nb]),2))
		for s in op_lipids_handled["lower"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(numpy.std(lipids_op_nff["data"][tail]["lower"][s][f_nb]),2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def op_xvg_nff_write_smoothed():
	
	#lipids in upper leaflet
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/smoothed/xvg/1_5_order_param_nff_upper_smoothed.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/smoothed/xvg/1_5_order_param_nff_upper_smoothed.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_5_order_param_nff_upper_smoothed.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of lipid tails order parameters in upper leaflet\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(2*len(op_lipids_handled["upper"])*3) + "\n")
	for s_index in range(0,len(op_lipids_handled["upper"])):
		output_xvg.write("@ s" + str(3*s_index) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail A (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+1) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail B (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+2) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " both (avg)\"\n")
		output_txt.write("1_3_order_param_nff_upper_smoothed.xvg," + str((3*s_index)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail A (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper_smoothed.xvg," + str((3*s_index+1)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail B (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper_smoothed.xvg," + str((3*s_index+2)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " both (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
	for s_index in range(0,len(op_lipids_handled["upper"])):
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["upper"])+3*s_index) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail A (std)\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["upper"])+3*s_index+1) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " tail (std)B\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["upper"])+3*s_index+2) + " legend \"" + str(op_lipids_handled["upper"][s_index]) + " both (std)\"\n")
		output_txt.write("1_3_order_param_nff_upper_smoothed.xvg," + str(3*len(op_lipids_handled["upper"])+(3*s_index)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail A (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper_smoothed.xvg," + str(3*len(op_lipids_handled["upper"])+(3*s_index+1)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " tail B (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_upper_smoothed.xvg," + str(3*len(op_lipids_handled["upper"])+(3*s_index+2)+1) +"," + str(op_lipids_handled["upper"][s_index]) + " both (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["upper"][s_index])]) + "\n")
	output_txt.close()
	for t_index in range(0, numpy.size(time_smooth)):
		results=str(time_smooth[t_index])
		for s in op_lipids_handled["upper"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(lipids_op_nff["smoothed"]["avg"][tail]["upper"][s][t_index],2))
		for s in op_lipids_handled["upper"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(lipids_op_nff["smoothed"]["std"][tail]["upper"][s][t_index],2))
		output_xvg.write(results)
	output_xvg.close()

	#lipids in lower leaflet
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/smoothed/xvg/1_5_order_param_nff_lower_smoothed.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/smoothed/xvg/1_5_order_param_nff_lower_smoothed.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_5_order_param_nff_lower_smoothed.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Evolution of lipid tails order parameters in lower leaflet\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(2*len(op_lipids_handled["lower"])*3) + "\n")
	for s_index in range(0,len(op_lipids_handled["lower"])):
		output_xvg.write("@ s" + str(3*s_index) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail A (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+1) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail B (avg)\"\n")
		output_xvg.write("@ s" + str(3*s_index+2) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " both (avg)\"\n")
		output_txt.write("1_3_order_param_nff_lower_smoothed.xvg," + str((3*s_index)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail A (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower_smoothed.xvg," + str((3*s_index+1)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail B (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower_smoothed.xvg," + str((3*s_index+2)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " both (avg)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
	for s_index in range(0,len(op_lipids_handled["lower"])):
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["lower"])+3*s_index) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail A (std)\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["lower"])+3*s_index+1) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " tail (std)B\"\n")
		output_xvg.write("@ s" + str(3*len(op_lipids_handled["lower"])+3*s_index+2) + " legend \"" + str(op_lipids_handled["lower"][s_index]) + " both (std)\"\n")
		output_txt.write("1_3_order_param_nff_lower_smoothed.xvg," + str(3*len(op_lipids_handled["lower"])+(3*s_index)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail A (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower_smoothed.xvg," + str(3*len(op_lipids_handled["lower"])+(3*s_index+1)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " tail B (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
		output_txt.write("1_3_order_param_nff_lower_smoothed.xvg," + str(3*len(op_lipids_handled["lower"])+(3*s_index+2)+1) +"," + str(op_lipids_handled["lower"][s_index]) + " both (std)," + mcolors.rgb2hex(colours_lipids[str(op_lipids_handled["lower"][s_index])]) + "\n")
	output_txt.close()
	for t_index in range(0, numpy.size(time_smooth)):
		results=str(time_smooth[t_index])
		for s in op_lipids_handled["lower"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(lipids_op_nff["smoothed"]["avg"][tail]["lower"][s][t_index],2))
		for s in op_lipids_handled["lower"]:
			for tail in ["tailA", "tailB", "tails"]:
				results += "	" + str(round(lipids_op_nff["smoothed"]["std"][tail]["lower"][s][t_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def op_xvg_nff_graph():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/png/1_2_order_param_nff.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/1_2_order_param_nff.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of lipid tails order parameter")
					
	#plot data: upper leafet
	#-----------------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	for s in op_lipids_handled["upper"]:
		p_upper[s] = plt.plot(time_sorted, lipids_op_nff["sorted"]["avg"]["tails"]["upper"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		p_upper[str(s + "_err")] = plt.fill_between(time_sorted, numpy.asarray(lipids_op_nff["sorted"]["avg"]["tails"]["upper"][s]) - numpy.asarray(lipids_op_nff["sorted"]["std"]["tails"]["upper"][s]), numpy.asarray(lipids_op_nff["sorted"]["avg"]["tails"]["upper"][s]) + numpy.asarray(lipids_op_nff["sorted"]["std"]["tails"]["upper"][s]), color=colours_lipids[s], alpha=0.2)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("upper leaflet", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('order parameter', fontsize="small")
	
	#plot data: lower leafet
	#-----------------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	for s in op_lipids_handled["lower"]:
		p_lower[s] = plt.plot(time_sorted, lipids_op_nff["sorted"]["avg"]["tails"]["lower"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		p_lower[str(s + "_err")] = plt.fill_between(time_sorted, numpy.asarray(lipids_op_nff["sorted"]["avg"]["tails"]["lower"][s]) - numpy.asarray(lipids_op_nff["sorted"]["std"]["tails"]["lower"][s]), numpy.asarray(lipids_op_nff["sorted"]["avg"]["tails"]["lower"][s]) + numpy.asarray(lipids_op_nff["sorted"]["std"]["tails"]["lower"][s]), color=colours_lipids[s], alpha=0.2)
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("lower leaflet", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('order parameter', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(-0.5, 1)
	ax2.set_ylim(-0.5, 1)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=7))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def op_xvg_nff_graph_smoothed():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/smoothed/png/1_4_order_param_nff_smoothed.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/smoothed/1_4_order_param_nff_smoothed.svg'
	
	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of lipid tails order parameter")
					
	#plot data: upper leafet
	#-----------------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	for s in op_lipids_handled["upper"]:
		p_upper[s] = plt.plot(time_smoothed, lipids_op_nff["smoothed"]["avg"]["tails"]["upper"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		p_upper[str(s + "_err")] = plt.fill_between(time_smoothed, numpy.asarray(lipids_op_nff["smoothed"]["avg"]["tails"]["upper"][s]) - numpy.asarray(lipids_op_nff["smoothed"]["std"]["tails"]["upper"][s]), numpy.asarray(lipids_op_nff["smoothed"]["avg"]["tails"]["upper"][s]) + numpy.asarray(lipids_op_nff["smoothed"]["std"]["tails"]["upper"][s]), color=colours_lipids[s], alpha=0.2)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("upper leaflet", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('order parameter', fontsize="small")
	
	#plot data: lower leafet
	#-----------------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	for s in op_lipids_handled["lower"]:
		p_lower[s] = plt.plot(time_smoothed, lipids_op_nff["smoothed"]["avg"]["tails"]["lower"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		p_lower[str(s + "_err")] = plt.fill_between(time_smoothed, numpy.asarray(lipids_op_nff["smoothed"]["avg"]["tails"]["lower"][s]) - numpy.asarray(lipids_op_nff["smoothed"]["std"]["tails"]["lower"][s]), numpy.asarray(lipids_op_nff["smoothed"]["avg"]["tails"]["lower"][s]) + numpy.asarray(lipids_op_nff["smoothed"]["std"]["tails"]["lower"][s]), color=colours_lipids[s], alpha=0.2)
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("lower leaflet", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('order parameter', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(-0.5, 1)
	ax2.set_ylim(-0.5, 1)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=7))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def op_frame_write_stat(f_nb, f_time):

	#nff lipids
	#==========
	#create file
	if f_nb == "all frames":
		filename_details = os.getcwd() + '/' + str(args.output_folder) + '/order_param/1_nff/1_1_order_param_nff.stat'
	else:
		filename_details = os.getcwd() + '/' + str(args.output_folder) + '/order_param/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_orderparam_' + str(int(f_time)).zfill(5) + 'ns_nff.stat'
	output_stat = open(filename_details, 'w')		
	output_stat.write("[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("1. membrane composition\n")
	output_stat.write(membrane_comp["upper"] + "\n")
	output_stat.write(membrane_comp["lower"] + "\n")
	tmp_string=str(op_lipids_handled["both"][0])
	for s in op_lipids_handled["both"][1:]:
		tmp_string+=", " + str(s)
	output_stat.write("\n")
	output_stat.write("2. lipid species processed: " + str(tmp_string) + "\n")
	if args.xtcfilename != "no":
		output_stat.write("\n")
		output_stat.write("3. nb frames processed:	" + str(nb_frames_processed) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
	if f_nb != "all frames":	
		output_stat.write("\n")
		output_stat.write("4. time: " + str(f_time) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")
	output_stat.write("\n")
	output_stat.write("lipid orientation with bilayer normal\n")
	output_stat.write(" P2=1    : parallel\n")
	output_stat.write(" P2=0    : random\n")
	output_stat.write(" P2=-0.5 : orthogonal\n")
		
	#lipids in upper leaflet
	output_stat.write("\n")
	output_stat.write("upper leaflet\n")
	output_stat.write("=============\n")
	output_stat.write("avg	nb	tail A	tail B	 both\n")
	output_stat.write("-------------------------------------\n")
	for s in op_lipids_handled["upper"]:
		tmp_output = str(s) + "	" + str(leaflet_sele["upper"][s].numberOfResidues())
		for tail in ["tailA", "tailB", "tails"]:
			tmp_output += "	" + str(round(numpy.average(lipids_op_nff["data"][tail]["upper"][s][f_nb]),2))
		output_stat.write(tmp_output + "\n")
	output_stat.write("\n")
	output_stat.write("std	nb	tail A	tail B	 both\n")
	output_stat.write("-------------------------------------\n")
	for s in op_lipids_handled["upper"]:
		tmp_output = str(s) + "	" + str(leaflet_sele["upper"][s].numberOfResidues())
		for tail in ["tailA", "tailB", "tails"]:
			tmp_output += "	" + str(round(numpy.std(lipids_op_nff["data"][tail]["upper"][s][f_nb]),2))
		output_stat.write(tmp_output + "\n")

	#lipids in lower leaflet
	output_stat.write("\n")
	output_stat.write("lower leaflet\n")
	output_stat.write("=============\n")
	output_stat.write("avg	nb	tail A	tail B	 both\n")
	output_stat.write("-------------------------------------\n")
	for s in op_lipids_handled["lower"]:
		tmp_output = str(s) + "	" + str(leaflet_sele["lower"][s].numberOfResidues())
		for tail in ["tailA", "tailB", "tails"]:
			tmp_output += "	" + str(round(numpy.average(lipids_op_nff["data"][tail]["lower"][s][f_nb]),2))
		output_stat.write(tmp_output + "\n")
	output_stat.write("\n")
	output_stat.write("std	nb	tail A	tail B	 both\n")
	output_stat.write("-------------------------------------\n")
	for s in op_lipids_handled["lower"]:
		tmp_output = str(s) + "	" + str(leaflet_sele["lower"][s].numberOfResidues())
		for tail in ["tailA", "tailB", "tails"]:
			tmp_output += "	" + str(round(numpy.std(lipids_op_nff["data"][tail]["lower"][s][f_nb]),2))
		output_stat.write(tmp_output + "\n")
	output_stat.close()

	#ff lipids
	#=========
	if args.selection_file_ff!="no":
		filename_details=os.getcwd() + '/' + str(args.output_folder) + '/order_param/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_orderparam_' + str(int(f_time)).zfill(5) + 'ns_ff.stat'
		output_stat = open(filename_details, 'w')		
		output_stat.write("[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_stat.write("\n")
	
		#general info
		output_stat.write("1. membrane composition\n")
		output_stat.write(membrane_comp["upper"] + "\n")
		output_stat.write(membrane_comp["lower"] + "\n")
		tmp_string=str(op_lipids_handled["both"][0])
		for s in op_lipids_handled["both"][1:]:
			tmp_string+=", " + str(s)
		output_stat.write("\n")
		output_stat.write("2. lipid species processed: " + str(tmp_string) + "\n")
		output_stat.write("\n")
		output_stat.write("3. nb frames processed:	" + str(nb_frames_processed) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		output_stat.write("\n")
		output_stat.write("4. time: " + str(f_time) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")
		output_stat.write("\n")
		output_stat.write("lipid orientation with bilayer normal\n")
		output_stat.write(" P2=1    : parallel\n")
		output_stat.write(" P2=0    : random\n")
		output_stat.write(" P2=-0.5 : orthogonal\n")
	
		#upper to lower
		if numpy.size(lipids_ff_u2l_index)>0:
			output_stat.write("\n")
			output_stat.write("upper to lower\n")
			output_stat.write("==============\n")
			output_stat.write("specie	resid	tail A	tail B  both\n")
			output_stat.write("-------------------------------------\n")
			for l in lipids_ff_u2l_index:
				output_stat.write(str(lipids_ff_info[l][0]) + "	" + str(lipids_ff_info[l][1]) + "	" + str(round(lipids_op_ff_tailA[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tailB[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tails[l][f_nb],2)) + "\n")
		
		#lower to upper
		if numpy.size(lipids_ff_l2u_index)>0:
			output_stat.write("\n")
			output_stat.write("lower to upper\n")
			output_stat.write("==============\n")
			output_stat.write("specie	resid	tail A	tail B  both\n")
			output_stat.write("-------------------------------------\n")
			for l in lipids_ff_l2u_index:
				output_stat.write(str(lipids_ff_info[l][0]) + "	" + str(lipids_ff_info[l][1]) + "	" + str(round(lipids_op_ff_tailA[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tailB[l][f_nb],2)) + "	" + str(round(lipids_op_ff_tails[l][f_nb],2)) + "\n")
		output_stat.close()
	
	return
def op_frame_write_snapshot(f_nb, f_time):

	#store order parameter info in beta factor field: nff lipids
	for l in ["lower","upper"]:
		for s in op_lipids_handled[l]:
			for r_index in range(0,leaflet_sele_atoms[l][s].numberOfResidues()):
				lipids_sele_nff[l][s][r_index].set_bfactor(lipids_op_nff[l][s][r_index][f_nb])
	
	#store order parameter info in beta factor field: ff lipids
	if args.selection_file_ff != "no":
		for l in range(0,lipids_ff_nb):
			lipids_sele_ff[l].set_bfactor(lipids_op_ff_tails[l][f_nb])

	#case: gro file
	if args.xtcfilename == "no":
		all_atoms.write(os.getcwd() + '/' + str(args.output_folder) + '/order_param/2_snapshots/' + args.grofilename[:-4] + '_annotated_orderparam', format="PDB")

	#case: xtc file
	else:
		tmp_name = os.getcwd() + "/" + str(args.output_folder) + '/order_param/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_orderparam_' + str(int(f_time)).zfill(5) + 'ns.pdb'
		W=Writer(tmp_name, nb_atoms)
		W.write(all_atoms)
	
	return
def op_frame_write_annotation(f_nb, f_time):
	
	#create file
	if args.xtcfilename == "no":
		filename_details = os.getcwd() + "/" + str(args.output_folder) + '/order_param/2_snapshots/' + args.grofilename[:-4] + '_annotated_orderparam.txt'
	else:
		filename_details = os.getcwd() + "/" + str(args.output_folder) + '/order_param/2_snapshots/' + args.xtcfilename[:-4] + '_annotated_orderparam_' + str(int(f_time)).zfill(5) + 'ns.txt'
	output_stat = open(filename_details, 'w')		

	#output selection strings: nff lipids
	tmp_sele_string = ""
	for l in ["lower","upper"]:
		for s in op_lipids_handled[l]:
			for r_index in range(0,leaflet_sele[l][s].numberOfResidues()):
				tmp_sele_string +="." + lipids_sele_nff_VMD_string[l][s][r_index]
	
	#output selection strings: ff lipids
	if args.selection_file_ff != "no":
		for l in range(0,lipids_ff_nb):
			tmp_sele_string += "." + lipids_sele_ff_VMD_string[l]
	output_stat.write(tmp_sele_string[1:] + "\n")
	
	#ouptut order param for each lipid for current frame: nff lipids
	tmp_ops = "1"
	for l in ["lower","upper"]:
		for s in op_lipids_handled[l]:
			for r_index in lipids_op_nff[l][s]:
				tmp_ops += ";" + str(round(lipids_op_nff[l][s][r_index][f_nb],2))
	
	#ouptut order param for each lipid for current frame: ff lipids
	if args.selection_file_ff != "no":
		for l in range(0,lipids_ff_nb):
			tmp_ops += ";" + str(round(lipids_op_ff_tails[l][f_nb],2))
	output_stat.write(tmp_ops + "\n")
	output_stat.close()

	return
def op_xtc_write_annotation():
	
	#create file
	filename_details = os.getcwd() + '/' + str(args.output_folder) + '/order_param/3_VMD/' + args.xtcfilename[:-4] + '_annotated_orderparam_dt' + str(args.frames_dt) + '.txt'
	output_stat = open(filename_details, 'w')		

	#output selection strings
	#------------------------
	#nff lipids
	tmp_sele_string=""
	for l in ["lower","upper"]:
		for s in op_lipids_handled[l]:
			for r_index in range(0,leaflet_sele_atoms[l][s].numberOfResidues()):
				tmp_sele_string+="." + lipids_sele_nff_VMD_string[l][s][r_index]
	#ff lipids
	if args.selection_file_ff!="no":
		for l in range(0,lipids_ff_nb):
			tmp_sele_string+="." + lipids_sele_ff_VMD_string[l]
	output_stat.write(tmp_sele_string[1:] + "\n")
	
	#ouptut order param for each lipid
	#---------------------------------
	for f_nb in sorted(time_stamp.iterkeys()):
		tmp_ops = str(f_nb)
		#nff lipids
		for l in ["lower","upper"]:
			for s in op_lipids_handled[l]:
				for r_index in lipids_op_nff[l][s]:
					tmp_ops += ";" + str(round(lipids_op_nff[l][s][r_index][f_nb],2))
		#ff lipids
		if args.selection_file_ff!="no":
			for l in range(0,lipids_ff_nb):
				tmp_ops += ";" + str(round(lipids_op_ff_tails[l][f_nb],2))
		output_stat.write(tmp_ops + "\n")
	output_stat.close()

	return

#annotations
def xtc_write_snapshots():

	#NB: - this will always output the first and final frame snapshots
	#    - it will also output intermediate frames according to the -w option	

	loc_nb_frames_processed=0
	for ts in U.trajectory:
		#case: frames before specified time boundaries
		if ts.time/float(1000)<args.t_start:
			progress = '\r -skipping frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)

		#case: frames within specified time boundaries
		elif ts.time/float(1000) > args.t_start and ts.time/float(1000) < args.t_end:
			progress = '\r -writing snapshots...   frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			if ((ts.frame-1) % args.frames_dt) == 0:
				if ((loc_nb_frames_processed) % args.frames_write_dt) == 0 or loc_nb_frames_processed == nb_frames_processed-1:
					#density
					#-------
					if args.radial:
						radial_density_frame_xvg_write(ts.frame, ts.time/float(1000))
						radial_density_frame_xvg_graph(ts.frame, ts.time/float(1000))
					#thickness
					#---------
					if args.perturb == 1 or args.perturb == 3:
						#statistics
						thick_frame_write_stat(ts.frame, ts.time/float(1000))
						#vmd annotations
						thick_frame_write_snapshot(ts.frame, ts.time/float(1000))
						thick_frame_write_annotation(ts.frame, ts.time/float(1000))
						#radial plots
						if args.radial:
							radial_thick_frame_xvg_write(ts.frame, ts.time/float(1000))
							radial_thick_frame_xvg_graph(ts.frame, ts.time/float(1000))
					#order parameter
					#---------------
					if args.perturb == 2 or args.perturb == 3:
						#statistics
						op_frame_write_stat(ts.frame, ts.time/float(1000))
						#vmd annotations
						op_frame_write_snapshot(ts.frame, ts.time/float(1000))
						op_frame_write_annotation(ts.frame, ts.time/float(1000))
						#radial plots
						if args.radial:
							radial_op_frame_xvg_write(ts.frame, ts.time/float(1000))
							radial_op_frame_xvg_graph(ts.frame, ts.time/float(1000))
				loc_nb_frames_processed += 1
		
		#case: frames after specified time boundaries
		elif ts.time/float(1000)>args.t_end:
			break
	print ''

	return

################################################################################################################################################
# FUNCTIONS: radial plots
################################################################################################################################################

#density #rename snapshots files
def radial_density_frame_xvg_write(f_nb, f_time):	
	global radial_step
	
	#species: summary of sizes influence for each specie
	#===================================================
	for s in leaflet_species["both"]:
		tmp_leaflets = []
		for l in ["lower","upper"]:
			if s in leaflet_species[l]:
				tmp_leaflets.append(l)
		if f_nb == "all frames":
			tmp_filename = 'radial_density_species_' + str(s)
		else:
			tmp_filename = 'radial_density_species_' + str(int(f_time)).zfill(5) + 'ns_' + str(s)
		if f_nb == "all frames":
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/species/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/species/xvg/' + str(tmp_filename) + '.xvg'
		else:
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/species/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/species/xvg/' + str(tmp_filename) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_txt.write("@[lipid density statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_filename) + ".xvg.\n")
		output_xvg = open(filename_xvg, 'w')
		output_xvg.write("@ title \"radial evolution of lipids density\n")
		output_xvg.write("@ xaxis  label \"distance from cluster center of geometry (Angstrom)\"\n")
		output_xvg.write("@ yaxis  label \"lipids density (%)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(tmp_leaflets)*len(radial_sizes[f_nb])) + "\n")
		for leaflet_index in range(0,len(tmp_leaflets)):
			#nb
			for c_index in range(0,len(radial_sizes[f_nb])):
				c_size = radial_sizes[f_nb][c_index]
				output_xvg.write("@ s" + str(leaflet_index*len(radial_sizes[f_nb]) + c_index) + " legend \"" + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (nb)\"\n")
				output_txt.write(str(tmp_filename) + ".xvg," + str(leaflet_index*len(radial_sizes[f_nb]) + c_index + 1) + "," + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (nb)," + mcolors.rgb2hex(colours_sizes[c_size]) + "\n")
			#%
			for c_index in range(0,len(radial_sizes[f_nb])):
				c_size = radial_sizes[f_nb][c_index]
				output_xvg.write("@ s" + str(leaflet_index*len(radial_sizes[f_nb]) + len(radial_sizes[f_nb]) + c_index) + " legend \"" + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (%)\"\n")
				output_txt.write(str(tmp_filename) + ".xvg," + str(leaflet_index*len(radial_sizes[f_nb]) + len(radial_sizes[f_nb]) + c_index + 1) + "," + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (%)," + mcolors.rgb2hex(colours_sizes[c_size]) + "\n")
		output_txt.close()
		for n in range(0,args.radial_nb_bins):
			results = str(n*radial_step)
			for leaflet_index in range(0,len(tmp_leaflets)):
				#nb
				for c_index in range(0,len(radial_sizes[f_nb])):
					c_size = radial_sizes[f_nb][c_index]
					if f_nb in radial_density[tmp_leaflets[leaflet_index]][s][c_size]["nb"].keys():
						results += "	" + str(radial_density[tmp_leaflets[leaflet_index]][s][c_size]["nb"][f_nb][n])
					else:
						results += "	0" 
				#%
				for c_index in range(0,len(radial_sizes[f_nb])):
					c_size = radial_sizes[f_nb][c_index]
					if f_nb in radial_density[tmp_leaflets[leaflet_index]][s][c_size]["nb"].keys():
						results += "	" + str(radial_density[tmp_leaflets[leaflet_index]][s][c_size]["pc"][f_nb][n])
					else:
						results += "	0" 
			output_xvg.write(results + "\n")
		output_xvg.close()

	#sizes: each species
	#===================
	for c_size in radial_sizes[f_nb] + ["all sizes"]:		
		if f_nb == "all frames":
			if c_size == "all sizes":
				tmp_filename = 'radial_density_sizes_all'
			else:
				tmp_filename = 'radial_density_sizes_' + str(c_size)
		else:
			if c_size == "all sizes":
				tmp_filename = 'radial_density_' + str(int(f_time)).zfill(5) + 'ns_sizes_all'
			else:
				tmp_filename = 'radial_density_' + str(int(f_time)).zfill(5) + 'ns_sizes_' + str(c_size)
		if f_nb == "all frames":
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/sizes/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/sizes/xvg/' + str(tmp_filename) + '.xvg'
		else:
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/sizes/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/sizes/xvg/' + str(tmp_filename) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_filename) +".xvg.\n")
		output_xvg = open(filename_xvg, 'w')
		output_xvg.write("@ title \"radial evolution of lipids density\n")
		output_xvg.write("@ xaxis  label \"distance from cluster center of geometry (Angstrom)\"\n")
		output_xvg.write("@ yaxis  label \"lipids density (%)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(leaflet_species["both"])) + "\n")
		#captions: lower leaflet
		#-----------------------
		#nb
		for s_index in range(0,len(leaflet_species["lower"])):
			s = leaflet_species["lower"][s_index]
			output_xvg.write("@ s" + str(s_index) + " legend \" lower" + str(s) + " (nb)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(s_index + 1) + ",lower" + str(s) + " (nb)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#%
		for s_index in range(0,len(leaflet_species["lower"])):
			s = leaflet_species["lower"][s_index]
			output_xvg.write("@ s" + str(len(leaflet_species["lower"]) + s_index) + " legend \" lower" + str(s) + " (%)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(len(leaflet_species["lower"]) + s_index + 1) + ",lower" + str(s) + " (%)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#captions: upper leaflet
		#-----------------------
		#nb
		for s_index in range(0,len(leaflet_species["upper"])):
			s = leaflet_species["upper"][s_index]
			output_xvg.write("@ s" + str(2*len(leaflet_species["lower"]) + s_index) + " legend \" upper" + str(s) + " (nb)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(len(leaflet_species["lower"]) + s_index + 1) + ",upper" + str(s) + " (nb)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#%
		for s_index in range(0,len(leaflet_species["upper"])):
			s = leaflet_species["upper"][s_index]
			output_xvg.write("@ s" + str(2*len(leaflet_species["lower"]) + len(leaflet_species["upper"]) + s_index) + " legend \" upper" + str(s) + " (nb)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(2*len(leaflet_species["lower"]) + len(leaflet_species["upper"]) + s_index + 1) + ",upper" + str(s) + " (nb)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		output_txt.close()
		#data
		#----
		for n in range(0,args.radial_nb_bins):
			results = str(n*radial_step)
			#data: lower leaflet
			for s_index in range(0,len(leaflet_species["lower"])):
				s = leaflet_species["lower"][s_index]
				if f_nb in radial_density["lower"][s][c_size]["nb"].keys():
					results += "	" + str(radial_density["lower"][s][c_size]["nb"][f_nb][n])
				else:
					results += "	0"
			for s_index in range(0,len(leaflet_species["lower"])):
				s = leaflet_species["lower"][s_index]
				if f_nb in radial_density["lower"][s][c_size]["nb"].keys():
					results += "	" + str(radial_density["lower"][s][c_size]["pc"][f_nb][n])
				else:
					results += "	0"
			
			#data: upper leaflet
			for s_index in range(0,len(leaflet_species["upper"])):
				s = leaflet_species["upper"][s_index]
				if f_nb in radial_density["upper"][s][c_size]["nb"].keys():
					results += "	" + str(radial_density["upper"][s][c_size]["nb"][f_nb][n])
				else:
					results += "	0"
			for s_index in range(0,len(leaflet_species["upper"])):
				s = leaflet_species["upper"][s_index]
				if f_nb in radial_density["upper"][s][c_size]["nb"].keys():				
					results += "	" + str(radial_density["upper"][s][c_size]["pc"][f_nb][n])
				else:
					results += "	0"
			output_xvg.write(results + "\n")
		output_xvg.close()	
		
	return
def radial_density_frame_xvg_graph(f_nb, f_time):
	
	global radial_step
	
	#species: summary of sizes influence for each specie
	#===================================================
	for s in leaflet_species["both"]:
		#create filenames
		if f_nb == "all frames":
			filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/species/png/radial_density_species' + str(s) + '.png'
			filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/species/radial_density_species' + str(s) + '.svg'
		else:
			filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/species/png/radial_density_species' + str(int(f_time)).zfill(5) + 'ns_' + str(s) + '.png'
			filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/species/radial_density_species' + str(int(f_time)).zfill(5) + 'ns_' + str(s) + '.svg'
	
		#create figure
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("radial evolution of lipids density")
		
		#create data
		loc_radial_bins=[]
		for n in range(0,args.radial_nb_bins):
			loc_radial_bins.append(n*radial_step)
		tmp_radial={}
		for l in ["lower","upper"]:
			tmp_radial[l] = {}
			for c_size in radial_sizes[f_nb]:
				tmp_radial[l][c_size]=numpy.zeros(args.radial_nb_bins)
				if s in leaflet_species[l]:
					if f_nb in radial_density[l][s][c_size]["pc"].keys():
						for n in range(0,args.radial_nb_bins):
							tmp_radial[l][c_size][n] = radial_density[l][s][c_size]["pc"][f_nb][n]
					else:
						for n in range(0,args.radial_nb_bins):
							tmp_radial[l][c_size][n] = numpy.nan

		#plot data: upper leafet
		ax1 = fig.add_subplot(211)
		p_upper={}
		if s in leaflet_species["upper"]:
			for c_size in radial_sizes[f_nb]:
				p_upper[c_size]=plt.plot(loc_radial_bins, tmp_radial["upper"][c_size], color=colours_sizes[c_size], linewidth=3.0, label=str(c_size))
			fontP.set_size("small")
			ax1.legend(prop=fontP)
		plt.title("upper leaflet", fontsize="small")
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('lipids density', fontsize="small")
		
		#plot data: lower leafet
		ax2 = fig.add_subplot(212)
		p_lower={}
		if s in leaflet_species["lower"]:
			for c_size in radial_sizes[f_nb]:
				p_lower[c_size]=plt.plot(loc_radial_bins, tmp_radial["lower"][c_size], color=colours_sizes[c_size], linewidth=3.0, label=str(c_size))
			fontP.set_size("small")
			ax2.legend(prop=fontP)
		plt.title("lower leaflet", fontsize="small")
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('lipids density', fontsize="small")
	
		#save figure
		ax1.set_xlim(0, args.radial_radius)
		ax1.set_ylim(0, 100)
		ax2.set_xlim(0, args.radial_radius)
		ax2.set_ylim(0, 100)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=10))
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=10))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	#sizes: each species
	#===================
	for c_size in radial_sizes[f_nb] + ["all sizes"]:
		#create filenames
		if f_nb == "all frames":
			if c_size == "all sizes":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/sizes/png/radial_density_sizes_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/sizes/radial_density_sizes_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/sizes/png/radial_density_sizes_' + str(c_size) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/sizes/radial_density_sizes_' + str(c_size) + '.svg'
		else:
			if c_size == "all sizes":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/sizes/png/radial_density_' + str(int(f_time)).zfill(5) + 'ns_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/sizes/radial_density_' + str(int(f_time)).zfill(5) + 'ns_all.svg'			
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/sizes/png/radial_density_' + str(int(f_time)).zfill(5) + 'ns_' + str(c_size) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/density/snapshots/sizes/radial_density_' + str(int(f_time)).zfill(5) + 'ns_' + str(c_size) + '.svg'

		#create figure
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("radial evolution of lipids density")
		
		#create data
		loc_radial_bins=[]
		for n in range(0,args.radial_nb_bins):
			loc_radial_bins.append(n*radial_step)
		tmp_radial = {}
		for l in ["lower","upper"]:
			tmp_radial[l]={}
			for s in leaflet_species[l]:
				tmp_radial[l][s] = numpy.zeros(args.radial_nb_bins)
				if f_nb in radial_density[l][s][c_size]["pc"].keys():
					for n in range(0,args.radial_nb_bins):
						tmp_radial[l][s][n] = radial_density[l][s][c_size]["pc"][f_nb][n]
				else:
					for n in range(0,args.radial_nb_bins):
						tmp_radial[l][s][n] = numpy.nan					

		#plot data: upper leafet
		ax1 = fig.add_subplot(211)
		p_upper={}
		for s in leaflet_species["upper"]:
			p_upper[s] = plt.plot(loc_radial_bins, tmp_radial["upper"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.title("upper leaflet", fontsize="small")
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('lipids density', fontsize="small")
		
		#plot data: lower leafet
		ax2 = fig.add_subplot(212)
		p_lower={}
		for s in leaflet_species["lower"]:
			p_lower[s] = plt.plot(loc_radial_bins, tmp_radial["lower"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
		fontP.set_size("small")
		ax2.legend(prop=fontP)
		plt.title("lower leaflet", fontsize="small")
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('lipids density', fontsize="small")
	
		#save figure
		ax1.set_xlim(0, args.radial_radius)
		ax1.set_ylim(0, 100)
		ax2.set_xlim(0, args.radial_radius)
		ax2.set_ylim(0, 100)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=10))
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=10))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	return

#thickness
def radial_thick_frame_xvg_write(f_nb, f_time):	
	
	global radial_step
	
	#species: summary of sizes influence for each specie
	#===================================================
	for s in leaflet_species["both"] + ["all species"]:
		#create filename
		if f_nb == "all frames":
			if s == "all species":
				tmp_filename = 'radial_thickness_species_all'
			else:
				tmp_filename = 'radial_thickness_species_' + str(s)
		else:
			if s == "all species":
				tmp_filename = 'radial_thickness_' + str(int(f_time)).zfill(5) + 'ns_species_all'
			else:
				tmp_filename = 'radial_thickness_' + str(int(f_time)).zfill(5) + 'ns_species_' + str(s)
		if f_nb == "all frames":
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/species/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/species/xvg/' + str(tmp_filename) + '.xvg'
		else:
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/species/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/species/xvg/' + str(tmp_filename) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_txt.write("@[bilayer thickness statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_filename) + ".xvg.\n")
		output_xvg = open(filename_xvg, 'w')
		output_xvg.write("@ title \"radial evolution of bilayer thickness\n")
		output_xvg.write("@ xaxis  label \"distance from cluster center of geometry (Angstrom)\"\n")
		output_xvg.write("@ yaxis  label \"bilayer thickness (Angstrom)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(radial_sizes[f_nb])) + "\n")
		#average values
		for c_index in range(0,len(radial_sizes[f_nb])):
			c_size = radial_sizes[f_nb][c_index]
			output_xvg.write("@ s" + str(c_index) + " legend \"" + str(c_size) + " (avg)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(c_index + 1) + "," + str(c_size) + " (avg)," + mcolors.rgb2hex(colours_sizes[c_size]) + "\n")
		#std values
		for c_index in range(0,len(radial_sizes[f_nb])):
			c_size = radial_sizes[f_nb][c_index]
			output_xvg.write("@ s" + str(len(radial_sizes[f_nb]) + len(radial_sizes[f_nb]) + c_index) + " legend \"" + str(c_size) + " (std)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(len(radial_sizes[f_nb]) + len(radial_sizes[f_nb]) + c_index + 1) + "," + str(c_size) + " (std)," + mcolors.rgb2hex(colours_sizes[c_size]) + "\n")
		output_txt.close()
		for n in range(0,args.radial_nb_bins):
			results = str(n*radial_step)
			#average values
			for c_index in range(0,len(radial_sizes[f_nb])):
				c_size = radial_sizes[f_nb][c_index]
				if f_nb in radial_thick[s]["avg"][c_size].keys():
					results += "	" + str(radial_thick[s]["avg"][c_size][f_nb][n])
				else:
					results += "	0"
			#std values
			for c_index in range(0,len(radial_sizes[f_nb])):
				c_size = radial_sizes[f_nb][c_index]
				if f_nb in radial_thick[s]["std"][c_size].keys():
					results += "	" + str(radial_thick[s]["std"][c_size][f_nb][n])
				else:
					results += "	0"
			output_xvg.write(results + "\n")
		output_xvg.close()

	#sizes: each species
	#===================
	for c_size in radial_sizes[f_nb] + ["all sizes"]:
		if f_nb == "all frames":
			if c_size == "all sizes":
				tmp_filename = 'radial_thickness_sizes_all'
			else:
				tmp_filename = 'radial_thickness_sizes_' + str(c_size)
		else:
			if c_size == "all sizes":
				tmp_filename = 'radial_thickness_' + str(int(f_time)).zfill(5) + 'ns_sizes_all'
			else:
				tmp_filename = 'radial_thickness_' + str(int(f_time)).zfill(5) + 'ns_sizes_' + str(c_size)
		if f_nb == "all frames":
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/xvg/' + str(tmp_filename) + '.xvg'
		else:
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/sizes/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/sizes/xvg/' + str(tmp_filename) + '.txt'
		filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/xvg/radial_thickness_' + str(c_size) + '.txt'
		filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/xvg/radial_thickness_' + str(c_size) + '.xvg'
		output_txt = open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_filename) + ".xvg.\n")
		output_xvg = open(filename_xvg, 'w')
		output_xvg.write("@ title \"radial evolution of bilayer thickness\n")
		output_xvg.write("@ xaxis  label \"distance from cluster center of geometry (Angstrom)\"\n")
		output_xvg.write("@ yaxis  label \"bilayer thickness (Angstrom)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(len(leaflet_species["both"])) + "\n")
		#avg values
		for s_index in range(0,len(leaflet_species["both"])):
			s = leaflet_species["lower"][s_index]
			output_xvg.write("@ s" + str(s_index) + " legend \"" + str(s) + " (avg)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(s_index + 1) + "," + str(s) + " (avg)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#std values
		for s_index in range(0,len(leaflet_species["both"])):
			s = leaflet_species["lower"][s_index]
			output_xvg.write("@ s" + str(len(leaflet_species["both"]) + s_index) + " legend \"" + str(s) + " (std)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(len(leaflet_species["both"]) + s_index + 1) + "," + str(s) + " (std)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		output_txt.close()
		#data
		#----
		for n in range(0,args.radial_nb_bins):
			results = str(n*radial_step)
			#avg values
			for s_index in range(0,len(leaflet_species["both"])):
				s = leaflet_species["both"][s_index]
				if f_nb in radial_thick[s]["avg"][c_size].keys():
					results += "	" + str(radial_thick[s]["avg"][c_size][f_nb][n])
				else:
					results += "	0"
			#std values
			for s_index in range(0,len(leaflet_species["lower"])):
				s = leaflet_species["both"][s_index]
				if f_nb in radial_thick[s]["std"][c_size].keys():
					results += "	" + str(radial_thick[s]["std"][c_size][f_nb][n])
				else:
					results += "	0"
			output_xvg.write(results + "\n")
		output_xvg.close()	
	
	return
def radial_thick_frame_xvg_graph(f_nb, f_time):

	global radial_step
	
	#species: summary of sizes influence for each specie
	#===================================================
	for s in leaflet_species["both"] + ["all species"]:
		#create filenames
		if f_nb == "all frames":
			if s == "all species":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/species/png/radial_thickness_species_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/species/radial_thickness_species_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/species/png/radial_thickness_species_' + str(s) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/species/radial_thickness_species_' + str(s) + '.svg'
		else:
			if s == "all species":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/species/png/radial_thickness_species_' + str(int(f_time)).zfill(5) + 'ns_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/species/radial_thickness_species_' + str(int(f_time)).zfill(5) + 'ns_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/species/png/radial_thickness_species_' + str(int(f_time)).zfill(5) + 'ns_' + str(s) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/species/radial_thickness_species_' + str(int(f_time)).zfill(5) + 'ns_' + str(s) + '.svg'
			
		#create figure
		fig=plt.figure(figsize=(8, 5))
		fig.suptitle("radial evolution of bilayer thickness")
		
		#create data
		loc_radial_bins=[]
		for n in range(0,args.radial_nb_bins):
			loc_radial_bins.append(n*radial_step)
		tmp_thick_avg={}
		tmp_thick_std={}
		for c_size in radial_sizes[f_nb] + ["all sizes"]:
			tmp_thick_avg[c_size] = numpy.zeros(args.radial_nb_bins)
			tmp_thick_std[c_size] = numpy.zeros(args.radial_nb_bins)
			if f_nb in radial_thick[s]["avg"][c_size].keys():
				for n in range(0,args.radial_nb_bins):
					tmp_thick_avg[c_size][n] = radial_thick[s]["avg"][c_size][f_nb][n]
					tmp_thick_std[c_size][n] = radial_thick[s]["std"][c_size][f_nb][n]
			else:
				for n in range(0,args.radial_nb_bins):
					tmp_thick_avg[c_size][n] = numpy.nan
					tmp_thick_std[c_size][n] = numpy.nan

		#plot data
		ax1 = fig.add_subplot(111)
		p_upper={}
		for c_size in radial_sizes[f_nb]:
			p_upper[c_size]=plt.plot(loc_radial_bins, tmp_thick_avg[c_size], color=colours_sizes[c_size], linewidth=3.0, label=str(c_size))
			p_upper[str(c_size) + "_err"]=plt.fill_between(loc_radial_bins, tmp_thick_avg[c_size]-tmp_thick_std[c_size], tmp_thick_avg[c_size]+tmp_thick_std[c_size], color=colours_sizes[c_size], alpha=0.2)
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('bilayer thickness', fontsize="small")
		
		#save figure
		ax1.set_xlim(0, args.radial_radius)
		ax1.set_ylim(numpy.min(lipids_thick_nff["data"]["all species"]["all frames"]), numpy.max(lipids_thick_nff["data"]["all species"]["all frames"]))
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	#sizes: each species
	#===================
	for c_size in radial_sizes[f_nb] + ["all sizes"]:
		#create filenames
		if f_nb == "all frames":
			if c_size == "all sizes":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/png/radial_thickness_sizes_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/radial_thickness_sizes_all.svg'		
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/png/radial_thickness_sizes_' + str(c_size) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/sizes/radial_thickness_sizes_' + str(c_size) + '.svg'
		else:
			if c_size == "all sizes":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/sizes/png/radial_thickness_sizes_' + str(int(f_time)).zfill(5) + 'ns_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/sizes/radial_thickness_sizes_' + str(int(f_time)).zfill(5) + 'ns_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/sizes/png/radial_thickness_sizes_' + str(int(f_time)).zfill(5) + 'ns_' + str(c_size) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/thickness/snapshots/sizes/radial_thickness_sizes_' + str(int(f_time)).zfill(5) + 'ns_' + str(c_size) + '.svg'

		#create figure
		fig=plt.figure(figsize=(8, 5))
		fig.suptitle("radial evolution of bilayer thickness")
		
		#create data
		loc_radial_bins=[]
		for n in range(0,args.radial_nb_bins):
			loc_radial_bins.append(n*radial_step)
		tmp_thick_avg={}
		tmp_thick_std={}
		for s in leaflet_species["both"]:
			tmp_thick_avg[s]=numpy.zeros(args.radial_nb_bins)
			tmp_thick_std[s]=numpy.zeros(args.radial_nb_bins)			
			if f_nb in radial_thick[s]["avg"][c_size].keys():
				for n in range(0,args.radial_nb_bins):
					tmp_thick_avg[s][n] = radial_thick[s]["avg"][c_size][f_nb][n]
					tmp_thick_std[s][n] = radial_thick[s]["std"][c_size][f_nb][n]
			else:
				for n in range(0,args.radial_nb_bins):
					tmp_thick_avg[s][n] = numpy.nan
					tmp_thick_std[s][n] = numpy.nan
						
		#plot data: upper leafet
		ax1 = fig.add_subplot(111)
		p_upper={}
		for s in leaflet_species["both"]:
			p_upper[s]=plt.plot(loc_radial_bins, tmp_thick_avg[s], color=colours_lipids[s], linewidth=3.0, label=str(s))
			p_upper[str(s + "_err")]=plt.fill_between(loc_radial_bins, tmp_thick_avg[s]-tmp_thick_std[s], tmp_thick_avg[s]+tmp_thick_std[s], color=colours_lipids[s], alpha=0.2)
		fontP.set_size("small")
		ax1.legend(prop=fontP)
		plt.xlabel('distance from cluster center of geometry($\AA$)', fontsize="small")
		plt.ylabel('bilayer thickness ($\AA$)', fontsize="small")
			
		#save figure
		ax1.set_xlim(0, args.radial_radius)		
		ax1.set_ylim(numpy.min(lipids_thick_nff["data"]["all species"]["all frames"]), numpy.max(lipids_thick_nff["data"]["all species"]["all frames"]))
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	return

#order parameters
def radial_op_frame_xvg_write(f_nb, f_time):	
	
	global radial_step
	
	#species: summary of sizes influence for each specie
	#===================================================
	for s in op_lipids_handled["both"] + ["all species"]:
		#find out in which leaflets the specie is present
		tmp_leaflets = []
		if s == "all species":
			for l in ["lower","upper"]:
				if len(op_lipids_handled[l]) > 0:
					tmp_leaflets.append(l)
		else:
			for l in ["lower","upper"]:
				if s in op_lipids_handled[l]:
					tmp_leaflets.append(l)
		#create filename
		if f_nb == "all frames":
			if s == "all species":
				tmp_filename = 'radial_order_param_species_all'
			else:
				tmp_filename = 'radial_order_param_species_' + str(s)
		else:
			if s == "all species":
				tmp_filename = 'radial_order_param_' + str(int(f_time)).zfill(5) + 'ns_species_all'
			else:
				tmp_filename = 'radial_order_param_' + str(int(f_time)).zfill(5) + 'ns_species_' + str(s)
		if f_nb == "all frames":
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/species/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/species/xvg/' + str(tmp_filename) + '.xvg'
		else:
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/species/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/species/xvg/' + str(tmp_filename) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_filename) + ".xvg.\n")
		output_xvg = open(filename_xvg, 'w')
		output_xvg.write("@ title \"radial evolution of lipid order parameters\n")
		output_xvg.write("@ xaxis  label \"distance from cluster center of geometry (Angstrom)\"\n")
		output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(2*len(tmp_leaflets)*len(radial_sizes[f_nb])) + "\n")
		for leaflet_index in range(0,len(tmp_leaflets)):
			#average values
			for c_index in range(0,len(radial_sizes[f_nb])):
				c_size = radial_sizes[f_nb][c_index]
				output_xvg.write("@ s" + str(leaflet_index*2*len(radial_sizes[f_nb]) + c_index) + " legend \"" + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (avg)\"\n")
				output_txt.write(str(tmp_filename) + ".xvg," + str(leaflet_index*2*len(radial_sizes[f_nb]) + c_index + 1) + "," + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (avg)," + mcolors.rgb2hex(colours_sizes[c_size]) + "\n")
			#std values
			for c_index in range(0,len(radial_sizes[f_nb])):
				c_size = radial_sizes[f_nb][c_index]
				output_xvg.write("@ s" + str(leaflet_index*2*len(radial_sizes[f_nb]) + len(radial_sizes[f_nb]) + c_index) + " legend \"" + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (std)\"\n")
				output_txt.write(str(tmp_filename) + ".xvg," + str(leaflet_index*2*len(radial_sizes[f_nb]) + len(radial_sizes[f_nb]) + c_index + 1) + "," + str(tmp_leaflets[leaflet_index]) + " " + str(c_size) + " (std)," + mcolors.rgb2hex(colours_sizes[c_size]) + "\n")
		output_txt.close()
		for n in range(0,args.radial_nb_bins):
			results = str(n*radial_step)
			for leaflet_index in range(0,len(tmp_leaflets)):
				#average values
				for c_index in range(0,len(radial_sizes[f_nb])):
					c_size = radial_sizes[f_nb][c_index]
					if f_nb in radial_op[tmp_leaflets[leaflet_index]][s]["avg"][c_size].keys():
						results += "	" + str(radial_op[tmp_leaflets[leaflet_index]][s]["avg"][c_size][f_nb][n])
					else:
						results += "	0"
				#std values
				for c_index in range(0,len(radial_sizes[f_nb])):
					c_size = radial_sizes[f_nb][c_index]
					if f_nb in radial_op[tmp_leaflets[leaflet_index]][s]["std"][c_size].keys():
						results += "	" + str(radial_op[tmp_leaflets[leaflet_index]][s]["std"][c_size][f_nb][n])
					else:
						results += "	0"
			output_xvg.write(results + "\n")
		output_xvg.close()

	#sizes: each species
	#===================
	for c_size in radial_sizes[f_nb] + ["all sizes"]:
		if f_nb == "all frames":
			if c_size == "all sizes":
				tmp_filename = 'radial_order_param_sizes_all'
			else:
				tmp_filename = 'radial_order_param_sizes_' + str(c_size)
		else:
			if c_size == "all sizes":
				tmp_filename = 'radial_order_param_' + str(int(f_time)).zfill(5) + 'ns_sizes_all'
			else:
				tmp_filename = 'radial_order_param_' + str(int(f_time)).zfill(5) + 'ns_sizes_' + str(c_size)
		if f_nb == "all frames":
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/xvg/' + str(tmp_filename) + '.xvg'
		else:
			filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/sizes/xvg/' + str(tmp_filename) + '.txt'
			filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/sizes/xvg/' + str(tmp_filename) + '.txt'
		filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/xvg/radial_order_param_' + str(c_size) + '.txt'
		filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/xvg/radial_order_param_' + str(c_size) + '.xvg'
		output_txt = open(filename_txt, 'w')
		output_txt.write("@[lipid tail order parameters statistics - written by bilayer_perturbations v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_filename) + ".xvg.\n")
		output_xvg = open(filename_xvg, 'w')
		output_xvg.write("@ title \"radial evolution of lipid order parameters\n")
		output_xvg.write("@ xaxis  label \"distance from cluster center of geometry (Angstrom)\"\n")
		output_xvg.write("@ yaxis  label \"order parameter P2\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(len(op_lipids_handled["both"])) + "\n")
		#captions: lower leaflet
		#-----------------------
		#avg values
		for s_index in range(0,len(op_lipids_handled["lower"])):
			s = op_lipids_handled["lower"][s_index]
			output_xvg.write("@ s" + str(s_index) + " legend \" lower" + str(s) + " (avg)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(s_index + 1) + ",lower" + str(s) + " (avg)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#std values
		for s_index in range(0,len(op_lipids_handled["lower"])):
			s = op_lipids_handled["lower"][s_index]
			output_xvg.write("@ s" + str(len(op_lipids_handled["lower"]) + s_index) + " legend \" lower" + str(s) + " (std)\"\n")
			output_txt.write(str(tmp_filename) + ".xvg," + str(len(op_lipids_handled["lower"]) + s_index + 1) + ",lower" + str(s) + " (std)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#captions: upper leaflet
		#-----------------------
		#avg values
		for s_index in range(0,len(op_lipids_handled["upper"])):
			s = op_lipids_handled["upper"][s_index]
			output_xvg.write("@ s" + str(2*len(op_lipids_handled["lower"]) + s_index) + " legend \" upper" + str(s) + " (avg)\"\n")
			output_txt.write("radial_order_param_" + str(c_size) + ".xvg," + str(2*len(op_lipids_handled["lower"]) + s_index + 1) + ",upper" + str(s) + " (avg)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		#std values
		for s_index in range(0,len(op_lipids_handled["upper"])):
			s = op_lipids_handled["upper"][s_index]
			output_xvg.write("@ s" + str(2*len(op_lipids_handled["upper"]) + s_index) + " legend \" upper" + str(s) + " (std)\"\n")
			output_txt.write("radial_order_param_" + str(c_size) + ".xvg," + str(2*len(op_lipids_handled["upper"]) + s_index + 1) + ",upper" + str(s) + " (std)," + mcolors.rgb2hex(colours_lipids[s]) + "\n")
		output_txt.close()
		#data
		#----
		for n in range(0,args.radial_nb_bins):
			results = str(n*radial_step)
			#data: lower leaflet
			#avg values
			for s_index in range(0,len(op_lipids_handled["lower"])):
				s = op_lipids_handled["lower"][s_index]
				if f_nb in radial_op["lower"][s]["avg"][c_size].keys():
					results += "	" + str(radial_op["lower"][s]["avg"][c_size][f_nb][n])
				else:
					results += "	0"
			#std values
			for s_index in range(0,len(op_lipids_handled["lower"])):
				s = op_lipids_handled["lower"][s_index]
				if f_nb in radial_op["lower"][s]["std"][c_size].keys():
					results += "	" + str(radial_op["lower"][s]["std"][c_size][f_nb][n])
				else:
					results += "	0"
			#data: upper leaflet
			#avg values
			for s_index in range(0,len(op_lipids_handled["upper"])):
				s = op_lipids_handled["upper"][s_index]
				if f_nb in radial_op["upper"][s]["avg"][c_size].keys():
					results += "	" + str(radial_op["upper"][s]["avg"][c_size][f_nb][n])
				else:
					results += "	0"
			#std values
			for s_index in range(0,len(op_lipids_handled["upper"])):
				s = op_lipids_handled["upper"][s_index]
				if f_nb in radial_op["upper"][s]["std"][c_size].keys():
					results += "	" + str(radial_op["upper"][s]["std"][c_size][f_nb][n])
				else:
					results += "	0"
			output_xvg.write(results + "\n")
		output_xvg.close()	
	
	return
def radial_op_frame_xvg_graph(f_nb, f_time):
	
	global radial_step
	
	#species: summary of sizes influence for each specie
	#===================================================
	for s in op_lipids_handled["both"] + ["all species"]:
		#create filenames
		if f_nb == "all frames":
			if s == "all species":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/species/png/radial_order_param_species_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/species/radial_order_param_species_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/species/png/radial_order_param_species_' + str(s) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/species/radial_order_param_species_' + str(s) + '.svg'
		else:
			if s == "all species":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/species/png/radial_order_param_species_' + str(int(f_time)).zfill(5) + 'ns_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/species/radial_order_param_species_' + str(int(f_time)).zfill(5) + 'ns_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/species/png/radial_order_param_species_' + str(int(f_time)).zfill(5) + 'ns_' + str(s) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/species/radial_order_param_species_' + str(int(f_time)).zfill(5) + 'ns_' + str(s) + '.svg'
			
		#create figure
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("radial evolution of lipid tails order parameter")
		
		#create data
		loc_radial_bins=[]
		for n in range(0,args.radial_nb_bins):
			loc_radial_bins.append(n*radial_step)
		tmp_op_avg={}
		tmp_op_std={}
		for l in ["lower","upper"]:
			tmp_op_avg[l]={}
			tmp_op_std[l]={}
			for c_size in radial_sizes[f_nb]:
				tmp_op_avg[l][c_size] = numpy.zeros(args.radial_nb_bins)
				tmp_op_std[l][c_size] = numpy.zeros(args.radial_nb_bins)
				if s in op_lipids_handled[l] or (s == "all species" and len(op_lipids_handled[l]) > 0):
					if f_nb in radial_op[l][s]["avg"][c_size].keys():
						for n in range(0,args.radial_nb_bins):
							tmp_op_avg[l][c_size][n] = radial_op[l][s]["avg"][c_size][f_nb][n]
							tmp_op_std[l][c_size][n] = radial_op[l][s]["std"][c_size][f_nb][n]
					else:
						for n in range(0,args.radial_nb_bins):
							tmp_op_avg[l][c_size][n] = numpy.nan
							tmp_op_std[l][c_size][n] = numpy.nan

		#plot data: upper leafet
		ax1 = fig.add_subplot(211)
		p_upper={}
		if s in op_lipids_handled["upper"] or ( s == "all species" and len(op_lipids_handled["upper"]) > 0):
			for c_size in radial_sizes[f_nb]:
				p_upper[c_size]=plt.plot(loc_radial_bins, tmp_op_avg["upper"][c_size], color=colours_sizes[c_size], linewidth=3.0, label=str(c_size))
				p_upper[str(c_size) + "_err"]=plt.fill_between(loc_radial_bins, tmp_op_avg["upper"][c_size]-tmp_op_std["upper"][c_size], tmp_op_avg["upper"][c_size]+tmp_op_std["upper"][c_size], color=colours_sizes[c_size], alpha=0.2)
			fontP.set_size("small")
			ax1.legend(prop=fontP)
		plt.title("upper leaflet", fontsize="small")
		plt.xlabel('distance from cluster ($\AA$)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
		
		#plot data: lower leafet
		ax2 = fig.add_subplot(212)
		p_lower={}
		if s in op_lipids_handled["lower"] or ( s == "all species" and len(op_lipids_handled["lower"]) > 0):
			for c_size in radial_sizes[f_nb]:
				p_lower[c_size]=plt.plot(loc_radial_bins, tmp_op_avg["lower"][c_size], color=colours_sizes[c_size], linewidth=3.0, label=str(c_size))
				p_lower[str(c_size) + "_err"]=plt.fill_between(loc_radial_bins, tmp_op_avg["lower"][c_size]-tmp_op_std["lower"][c_size], tmp_op_avg["lower"][c_size]+tmp_op_std["lower"][c_size], color=colours_sizes[c_size], alpha=0.2)
			fontP.set_size("small")
			ax2.legend(prop=fontP)
		plt.title("lower leaflet", fontsize="small")
		plt.xlabel('distance from cluster ($\AA$)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
	
		#save figure
		ax1.set_xlim(0, args.radial_radius)
		ax1.set_ylim(-0.5, 1)
		ax2.set_xlim(0, args.radial_radius)
		ax2.set_ylim(-0.5, 1)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=7))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	#sizes: each species
	#===================
	for c_size in radial_sizes[f_nb] + ["all sizes"]:
		#create filenames
		if f_nb == "all frames":
			if c_size == "all sizes":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/png/radial_order_param_sizes_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/radial_order_param_sizes_all.svg'		
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/png/radial_order_param_sizes_' + str(c_size) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/sizes/radial_order_param_sizes_' + str(c_size) + '.svg'
		else:
			if c_size == "all sizes":
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/sizes/png/radial_order_param_sizes_' + str(int(f_time)).zfill(5) + 'ns_all.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/sizes/radial_order_param_sizes_' + str(int(f_time)).zfill(5) + 'ns_all.svg'
			else:
				filename_png=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/sizes/png/radial_order_param_sizes_' + str(int(f_time)).zfill(5) + 'ns_' + str(c_size) + '.png'
				filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/radial/order_param/snapshots/sizes/radial_order_param_sizes_' + str(int(f_time)).zfill(5) + 'ns_' + str(c_size) + '.svg'

		#create figure
		fig=plt.figure(figsize=(8, 6.2))
		fig.suptitle("radial evolution of lipid tails order parameter")
		
		#create data
		loc_radial_bins=[]
		for n in range(0,args.radial_nb_bins):
			loc_radial_bins.append(n*radial_step)
		tmp_op_avg={}
		tmp_op_std={}
		for l in ["lower","upper"]:
			tmp_op_avg[l]={}
			tmp_op_std[l]={}
			for s in op_lipids_handled[l]:
				tmp_op_avg[l][s] = numpy.zeros(args.radial_nb_bins)
				tmp_op_std[l][s] = numpy.zeros(args.radial_nb_bins)
				if f_nb in radial_op[l][s]["avg"][c_size].keys():
					for n in range(0,args.radial_nb_bins):
						tmp_op_avg[l][s][n] = radial_op[l][s]["avg"][c_size][f_nb][n]
						tmp_op_std[l][s][n] = radial_op[l][s]["std"][c_size][f_nb][n]
				else:
					for n in range(0,args.radial_nb_bins):
						tmp_op_avg[l][s][n] = numpy.nan
						tmp_op_std[l][s][n] = numpy.nan

		#plot data: upper leafet
		ax1 = fig.add_subplot(211)
		p_upper={}
		for s in op_lipids_handled["upper"]:
			p_upper[s]=plt.plot(loc_radial_bins, tmp_op_avg["upper"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
			p_upper[str(s + "_err")]=plt.fill_between(loc_radial_bins, tmp_op_avg["upper"][s]-tmp_op_std["upper"][s], tmp_op_avg["upper"][s]+tmp_op_std["upper"][s], color=colours_lipids[s], alpha=0.2)
		if len(op_lipids_handled["upper"]) > 0:
			fontP.set_size("small")
			ax1.legend(prop=fontP)
		plt.title("upper leaflet", fontsize="small")
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
		
		#plot data: lower leafet
		ax2 = fig.add_subplot(212)
		p_lower={}
		for s in op_lipids_handled["lower"]:
			p_lower[s]=plt.plot(loc_radial_bins, tmp_op_avg["lower"][s], color=colours_lipids[s], linewidth=3.0, label=str(s))
			p_lower[str(s + "_err")]=plt.fill_between(loc_radial_bins, tmp_op_avg["lower"][s]-tmp_op_std["lower"][s], tmp_op_avg["lower"][s]+tmp_op_std["lower"][s], color=colours_lipids[s], alpha=0.2)
		if len(op_lipids_handled["lower"]) > 0:
			fontP.set_size("small")
			ax2.legend(prop=fontP)
		plt.title("lower leaflet", fontsize="small")
		plt.xlabel('distance from cluster center of geometry ($\AA$)', fontsize="small")
		plt.ylabel('order parameter', fontsize="small")
	
		#save figure
		ax1.set_xlim(0, args.radial_radius)
		ax1.set_ylim(-0.5, 1)
		ax2.set_xlim(0, args.radial_radius)
		ax2.set_ylim(-0.5, 1)
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=7))
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=args.radial_nb_bins))
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=7))
		plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
		plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
		plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()
	
	return

################################################################################################################################################
# ALGORITHM
################################################################################################################################################

#=======================================================================
#process inputs
#=======================================================================
set_forcefield_lipids_info()
load_MDA_universe()
if args.selection_file_ff != "no":
	identify_ff()
if args.selection_file_prot != "no":
	identify_proteins()
identify_leaflets()
identify_species()
initialise_colours_and_groups()

#=======================================================================
# initialise data structures
#=======================================================================
print "\nInitialising data structures..."
data_struct_time_coords()
#case: thickness
if args.perturb == 1 or args.perturb == 3:
	data_struct_thick()
#case: order parameter
if args.perturb == 2 or args.perturb == 3:
	data_struct_op_nff()
	if args.selection_file_ff != "no":
		data_struct_op_ff()
#case: radial
if args.radial:
	data_struct_radial()

#=======================================================================
# generate data
#=======================================================================
print "\nCalculating bilayer perturbations..."
#case: gro file
if args.xtcfilename=="no":
	time_stamp[1]=0
	if args.perturb == 1 or args.perturb == 3:
		calculate_thickness(1)
	if args.perturb == 2 or args.perturb == 3:
		calculate_order_parameters(1)
	if args.radial:
		calculate_radial(1)
#case: xtc file
else:
	for ts in U.trajectory:
		if ts.time/float(1000) < args.t_start:
			progress = '\r -skipping frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)

		elif ts.time/float(1000) > args.t_start and ts.time/float(1000) < args.t_end:
			progress = '\r -processing frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '                      '  
			sys.stdout.flush()
			sys.stdout.write(progress)
			if ((ts.frame-1) % args.frames_dt) == 0:
				nb_frames_processed += 1
				time_stamp[ts.frame] = ts.time/float(1000)
				get_z_coords(ts.frame)
				if args.perturb == 1 or args.perturb == 3:
					calculate_thickness(ts.frame)
				if args.perturb == 2 or args.perturb == 3:
					calculate_order_parameters(ts.frame)
				if args.radial:
					calculate_radial(ts.frame)
				nb_frames_processed += 1
		
		elif ts.time/float(1000)>args.t_end:
			break
	print ''

#=======================================================================
# process data
#=======================================================================
print "\nProcessing results..."
if args.radial:
	calculate_radial_data()
	update_colours_sizes()
		
#=======================================================================
# produce outputs
#=======================================================================
print "\nWriting outputs..."
#case: gro file
#--------------
if args.xtcfilename == "no":

	#writing statistics
	if args.perturb != 0:
		print " -writing statistics..."
		if args.perturb == 1 or args.perturb == 3:
			thick_frame_write_stat("all frames", 0)
		if args.perturb == 2 or args.perturb == 3:
			op_frame_write_stat("all frames", 0)

	#write annotation files for VMD
	print " -writing VMD annotation files..."
	if args.perturb == 1 or args.perturb == 3:
		thick_frame_write_snapshot(1, 0)
		thick_frame_write_annotation(1, 0)
	if args.perturb == 2 or args.perturb == 3:
		op_frame_write_snapshot(1, 0)
		op_frame_write_annotation(1, 0)

	#radials plots
	if args.radial:
		if len(radial_sizes["all frames"]) == 0:
			print "Warning: no TM cluster dected."
		else:
			print " -writing radials perturbations..."
			radial_density_frame_xvg_write("all frames", 0)
			radial_density_frame_xvg_graph("all frames", 0)
			if args.perturb == 1 or args.perturb == 3:
				radial_thick_frame_xvg_write("all frames", 0)
				radial_thick_frame_xvg_graph("all frames", 0)
			if args.perturb == 2 or args.perturb == 3:
				radial_op_frame_xvg_write("all frames", 0)
				radial_op_frame_xvg_graph("all frames", 0)
			if len(radial_sizes["all frames"]) == 1:
				print 
				print "Warning: a single TM cluster size (", str(radial_sizes["all frames"][0]), ") was detected. Check the cluster detection options (see bilayer_perturbations -h)."

#case: xtc file
#--------------
else:
	sort_and_smooth_data()
	#writing statistics
	if args.perturb != 0:
		print " -writing statistics..."
		if args.perturb == 1 or args.perturb == 3:
			thick_frame_write_stat("all frames", 0)
		if args.perturb == 2 or args.perturb == 3:
			op_frame_write_stat("all frames", 0)
	
	#output snapshots
	xtc_write_snapshots()
	
	#write annotation files for VMD
	print " -writing VMD xtc annotation files..."
	if args.perturb == 1 or args.perturb == 3:
		thick_xtc_write_annotation()
	if args.perturb == 2 or args.perturb == 3:
		op_xtc_write_annotation()
	
	#write xvg and graphs
	print " -writing xvg and graphs..."
	if args.perturb == 1 or args.perturb == 3:
		thick_xvg_write()
		thick_xvg_graph()	
		if args.nb_smoothing > 1:
			thick_xvg_write_smoothed()
			thick_xvg_graph_smoothed()		
	if args.perturb == 2 or args.perturb == 3:
		op_xvg_nff_write()
		op_xvg_nff_graph()	
		if args.nb_smoothing > 1:
			op_xvg_nff_write_smoothed()
			op_xvg_nff_graph_smoothed()		
		if args.selection_file_ff != "no":
			op_xvg_ff_write()
			op_xvg_ff_graph()	
			if args.nb_smoothing > 1:
				op_xvg_ff_write_smoothed()
				op_xvg_ff_graph_smoothed()	
	
	#radials plots
	if args.radial:
		if len(radial_sizes["all frames"]) == 0:
			print "Warning: no TM cluster dected."
		else:
			print " -writing radials perturbations..."
			radial_density_frame_xvg_write("all frames", 0)
			radial_density_frame_xvg_graph("all frames", 0)
			if args.perturb == 1 or args.perturb == 3:
				radial_thick_frame_xvg_write("all frames", 0)
				radial_thick_frame_xvg_graph("all frames", 0)
			if args.perturb == 2 or args.perturb == 3:
				radial_op_frame_xvg_write("all frames", 0)
				radial_op_frame_xvg_graph("all frames", 0)
			if len(radial_sizes["all frames"]) == 1:
				print 
				print "Warning: a single TM cluster size (", str(radial_sizes["all frames"][0]), ") was detected. Check the cluster detection options (see bilayer_perturbations -h)."
					
#=======================================================================
# exit
#=======================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
