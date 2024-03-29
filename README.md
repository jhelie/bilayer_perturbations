bilayer_perturbations
=====================

Python utility to calculate bilayer perturbations (lipids distribution, bilayer thickness, lipids tails orders parameters)

Requirements
------------
The following Python modules are required:
- MDAnalysis
- numpy
- scipy
- matplotlib
- networkX (if option --radial is used and option --algorithm set to 'min' or 'cog')
- sklearn (if option --radial is used and option --algorithm set to 'density')

Usage
-----
python bilayer_perturbations.py --help

Features
--------
- calculate evoution of bilayer thickness and lipids tail order parameters
- calculate radial influence of TM clusters on lipids distribution, bilayer thickness, lipids tails order parameters
- statistics and graphs broken by species / leaflet
- automatic detection of protein and TM clusters 
- output files for visualisation in VMD
- customisable lipids selection and definition
- different leaflet identifications method (1 suitable for very large systems)
- calculation of TM clusters center of geometry take pbc into account
- on the fly outputting

