#! /usr/bin/python
#
# Simple parser to extract contacts
# Chain ids or TER not required
#
__author__ = "gelpi"
__date__ = "$29-ago-2017 16:14:26$"

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from ForceField import VdwParamset
from ResLib import  ResiduesDataLib
from matplotlib import pyplot as plt



import pylab
import sys
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np

COVLNK = 2.0
HBLNK  = 3.5

all_polars = [
    'N', 'ND1', 'ND2', 'NE',  'NE1', 'NE2', 'NH1', 'NH2', 'NZ',
    'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG',  'OG1', 'OH',
    'S', 'SD',  'SG'
]
backbone_polars =  ['N','O']
waternames = ['WAT','HOH']

residues = ['N', 'ND1', 'ND2', 'NE',  'NE1', 'NE2', 'NH1', 'NH2', 'NZ',
    'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG',  'OG1', 'OH',
    'S', 'SD',  'SG', 'N','O', 'WAT','HOH']

def main():
	parser = argparse.ArgumentParser(
	prog='polarContacts',
	description='Polar contacts detector')

	parser.add_argument(
	'--backonly',
	action='store_true',
	dest='backonly',
	help='Restrict to backbone')

	parser.add_argument(
	'--nowats',
	action='store_true',
	dest='nowats',
	help='Exclude water molecules')
    
	parser.add_argument(
	'--diel',
	type= float,
	action='store',
	dest='diel',
	default = 1.0,
	help='Relative dielectric constant')
    
	parser.add_argument(
	'--vdw',
	action='store',
	dest='vdwprm',
	help='VDW Paramters file')
    
	parser.add_argument(
	'--rlib',
	action='store',
	dest='reslib',
	help='AminoAcid library')

	parser.add_argument('pdb_path')

	args = parser.parse_args()

	print ("Settings")
	print ("--------")
	for k,v in vars(args).items():
		print ('{:10}:'.format(k),v)

	backonly = args.backonly
	nowats =args.nowats
	pdb_path = args.pdb_path
	vdwprm = args.vdwprm
	reslib = args.reslib
	diel = args.diel
    
# Load VDW parameters
	vdwParams = VdwParamset(vdwprm)
	print ("{} atom types loaded".format(vdwParams.ntypes))

# Load AA Library
	aaLib = ResiduesDataLib(reslib)
	print ("{} amino acid atoms loaded".format(aaLib.nres))
    
	if not pdb_path:
		parser.print_help()
		sys.exit(2)

	parser = PDBParser(PERMISSIVE=1)

	try:
		st = parser.get_structure('st', pdb_path)
	except OSError:
		print ("#ERROR: loading PDB")
		sys.exit(2)

# Checking for models
	if len(st) > 1:
		print ("#WARNING: Several Models found, using only first")

# Using Model 0 any way
	st = st[0]

# Making a list of polar atoms
	polats = []
	if backonly:
		selected_atoms = backbone_polars
	else:
		selected_atoms = all_polars

	for at in st.get_atoms():
		if at.id in selected_atoms:
			polats.append(at)
#Searching for contacts under HNLNK on diferent residues
	nbsearch = NeighborSearch(polats)
	hblist = []
	for at1, at2 in nbsearch.search_all(HBLNK):
		if at1.get_parent() == at2.get_parent():
			continue
 #Discard covalents and neighbours
		if (at1-at2) < COVLNK:
			continue
		if abs(at2.get_parent().id[1] - at1.get_parent().id[1]) == 1:
			continue
# remove waters
		if nowats:
			if at1.get_parent().get_resname() in waternames \
                	or at2.get_parent().get_resname() in waternames:
				continue
                
   #     atom1 = Atom(at1,1,aaLib,vdwParams)
   #     atom2 = Atom(at2,1,aaLib,vdwParams)        
		if at1.get_serial_number() < at2.get_serial_number():
			hblist.append([at1, at2])
		else:
			hblist.append([at2, at1])
       
	print ()
	print ("Polar contacts")
	print ('{:13} {:13} {:6} '.format(
            'Atom1','Atom2','Dist (A)'))
    
	for hb in sorted (hblist,key=lambda i: i[0].get_serial_number()):
		r1 = hb[0].get_parent()
		r2 = hb[1].get_parent()
		print ('{:14} {:14} {:6.3f} '.format(
		r1.get_resname()+' '+str(r1.id[1])+hb[0].id,
		r2.get_resname()+' '+str(r2.id[1])+hb[1].id,
		hb[0] - hb[1]
		)
	)
		print ()
		print ("Residue interactions")
    

# Making list or residue pairs to avoid repeated pairs
	respairs = []
	for hb in hblist:
		r1 = hb[0].get_parent()
		r2 = hb[1].get_parent()
		if [r1,r2] not in respairs:
			respairs.append([r1,r2])
            
	l = []
    
	for rpair in sorted(respairs, key=lambda i: i[0].id[1]):            
		eint=0.
		evdw=0.
		
		for at1 in rpair[0].get_atoms():
			resid1 = rpair[0].get_resname()
			atid1 = at1.id
			atparam1 = aaLib.getParams(resid1,atid1)
			vdwprm1 = vdwParams.atTypes[atparam1.atType]
			
			for at2 in rpair[1].get_atoms():
				resid2 = rpair[1].get_resname()
				atid2 = at2.id
				atparam2 = aaLib.getParams(resid2,atid2)
				vdwprm2 = vdwParams.atTypes[atparam2.atType]
				eint = eint + 332.16 * atparam1.charg * atparam2.charg/diel/(at1-at2)
				eps = math.sqrt(vdwprm1.eps*vdwprm2.eps)
				sig = math.sqrt(vdwprm1.sig*vdwprm2.sig)
				evdw = evdw + 4 * eps *( (sig/(at1-at2))**12-(sig/(at1-at2))**6)
				
			print (resid1,rpair[0].id[1],resid2,rpair[1].id[1],eint,evdw, eint+evdw)            
			l.append([resid1,rpair[0].id[1],resid2,rpair[1].id[1],eint,evdw, eint+evdw])
			
	#here we have the code for finding the most stable contacts and plotting each energy component 
	#(global, electrostatic and vdw) with respect to the residue number involved in these contacts		
	print("Five most stable contacts")
	stable = []
	for index, element in enumerate(sorted(l, key=lambda i: i[6])):
		if index < 5:
			stable.append(element)
			print(element)

	n_groups = 5
	eint = (-96.879048334060371, -89.262401988309293,-63.650369322307412,-51.488465980661772,-50.345308360049728)
	vdw = (-1.0577685827505385, -1.5931226867662258,1.5038605656892994,-2.9601475910966375,-2.0108017155800604)
	etot = (-97.936816916810912, -90.855524675075515,-62.146508756618111,-54.448613571758408,-52.356110075629786)
	fig, ax = plt.subplots()
	index = np.arange(n_groups)
	bar_width = 0.15
	opacity = 0.5
	inf1 = plt.bar(index, eint, bar_width,alpha=opacity,color='b',label='electrostatic energies')
	inf2 = plt.bar(index+bar_width, vdw, bar_width,alpha=opacity,color='g',label='van der waals energies')
	inf3 = plt.bar(index+bar_width, etot, bar_width,alpha=opacity,color='r',label='total energies')
	plt.title('Energies for pair of residues')
	plt.xlabel('Contacts')
	plt.ylabel('Energy')
	plt.xticks(index + bar_width, ('LYS-ASP','LYS-GLU','GLU-LYS','GLU-ARG','ASP-ARG'))
	plt.legend()
	plt.tight_layout()
	plt.show()
		
	
	both_main_x = []
	both_main_y = []
	both_side_x = []
	both_side_y = []
	main_side_x = []
	main_side_y = []
	side_main_x = []
	side_main_y = []
	for hb in sorted (hblist, key=lambda i: i[0].get_serial_number()):
		if hb[0].id in backbone_polars:
			where0 = 'main'
		else:
			where0 = 'side'
		if hb[1].id in backbone_polars:
			where1 = 'main'
		else:
			where1 = 'side'
		label = where0+':'+where1
		if label[0] == label[5] and label[0] == 'm':
			value = 1
			both_main_x.append(hb[0].get_parent().id[1])
			both_main_y.append(hb[1].get_parent().id[1])
		elif label[0] == label[5] and label[0] == 's':
			value = 2
			both_side_x.append(hb[0].get_parent().id[1])
			both_side_y.append(hb[1].get_parent().id[1])
		elif label[0] != label[5] and label[0] == 'm':
			value = 3
			main_side_x.append(hb[0].get_parent().id[1])
			main_side_y.append(hb[1].get_parent().id[1])
		elif label[0] != label[5] and label[0] == 's':
			value = 4
			side_main_x.append(hb[0].get_parent().id[1])
			side_main_y.append(hb[1].get_parent().id[1])
		linking = [label,value,hb[0].id,hb[1].id,hb[0]-hb[1]]
		print ('{:14}{:14}{:14}{:14}{:6.3f}'.format(label,value,hb[0].id,hb[1].id,hb[0]-hb[1]))
	
	plt.figure(figsize=(10, 8))
	plt.scatter(both_main_x, both_main_y, c='red', label='both_main')
	plt.scatter(both_side_x, both_side_y, c='green', label='both_side')
	plt.scatter(main_side_x, main_side_y, c='blue', label='main_side')
	plt.scatter(side_main_x, side_main_y, c='yellow', label='side_main')
	plt.title('Interaction')
	plt.xlabel('Residue1 number')
	plt.ylabel('Residue2 number')
	plt.legend(loc='upper right')
	plt.show()
	
	#The surface residues are the ones with an Area<5 (http://cib.cf.ocha.ac.jp/bitool/ASA/display.php?id=1513152459.2996)
	surface_res=[['ILE',3],['VAL',5],['ILE',23],['VAL',26],['ILE',30],['GLN',41],['LEU',43],['LEU',56],['ILE',61],['LEU',67],['LEU',69]]
	
	
	for rpair in sorted(respairs, key=lambda i: i[0].id[1]):            
		eint=0.
		
		for atom1 in rpair[0].get_atoms():
			resname1 = rpair[0].get_resname()
			atid1 = at1.id
			atparam1 = aaLib.getParams(resid1,atid1)
			
			for atom2 in rpair[1].get_atoms():
				resname2 = rpair[1].get_resname()
				atid2 = at2.id
				atparam2 = aaLib.getParams(resid2,atid2)
				
				for values in surface_res:
					for values2 in surface_res:
						if resname1==values[0] and rpair[0].id[1]==values[1]:
							if resname2==values2[0] and rpair[1].id[1]==values2[1]:
								eint = eint + 80 * atparam1.charg * atparam2.charg/diel/(atom1-atom2)
		if eint!=0:
			print (resid1,rpair[0].id[1],resid2,rpair[1].id[1],eint,evdw, eint+evdw)

	
	
'''
A) Evaluate the electrostatic and vdw interaction energies between the residues involved in the
polar contacts. Proper evaluation of interaction energies requires to add hydrogen atoms to the
downloaded PDB files (See note below).

Obtain a list of the five most stable contacts. Indicate in a table the residues involved and the
global, electrostatic and vdw energies of these five contacts.

Make a representation of each energy component (global, electrostatic and vdw) with respect
to the residue number involved in these contacts.

residue1, pairnumber1, residue2, pairnumber2, electrostatic energies, van der waals energies, total energies
['LYS', 27, 'ASP', 52, -96.879048334060371, -1.0577685827505385, -97.936816916810912]
['LYS', 11, 'GLU', 34, -89.262401988309293, -1.5931226867662258, -90.855524675075515]
['GLU', 16, 'LYS', 29, -63.650369322307412, 1.5038605656892994, -62.146508756618111]
['GLU', 51, 'ARG', 54, -51.488465980661772, -2.9601475910966375, -54.448613571758408]
['ASP', 39, 'ARG', 72, -50.345308360049728, -2.0108017155800604, -52.356110075629786]
'''


# Making list or residue pairs to avoid repeated pairs
#    for hb in hblist:
#        r1 = Residue(hb[0].at.get_parent(), 1, aaLib, vdwParams)
#        r2 = Residue(hb[1].at.get_parent(), 1, aaLib, vdwParams)
#        if [r1,r2] not in respairs:
#            respairs.append([r1,r2])
# 
#    for rpair in sorted(respairs, key=lambda i: i[0].resNum()):            
#        eint = rpair[0].elecInt(rpair[1],diel)
#        evdw = rpair[0].vdwInt(rpair[1])
#        print (
#            '{:10} {:10} {: 8.4f} {: 8.4f} {: 8.4f}'.format(
#                rpair[0].resid(), 
#                rpair[1].resid(),
#                eint,
#                evdw,
#                eint+evdw)
#        )
        
if __name__ == "__main__":
	main()
