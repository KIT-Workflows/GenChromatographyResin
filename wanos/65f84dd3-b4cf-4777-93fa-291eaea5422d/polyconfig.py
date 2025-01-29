import sys
import os

#{{ wano["num_Monomer"] }} {{ wano["Monomer type for the polymeric backbone"] }}, {{ wano["num_Crosslinker"] }} {{ wano["Crosslinker type for the polymer network"] }} , and {{ wano["num_Ligand"] }} {{ wano["Ligand type"] }} , with spacer size {{ wano["Spacer units"] }} are mixed

nMonomer, tMonomer, nXLinker, tXLinker, nLigand, tLigand, nSpacer=sys.argv[1:]
if tLigand=="Tetramethylammonium":
    tLigand="Q"
elif tLigand=="Sulfonic_acid":
    tLigand="S"

with open("poly.config", "w") as fw:
    fw.write("{} {}, {} {}, and {} {}, with spacer size {} are mixed".format(nMonomer, tMonomer, nXLinker, tXLinker, nLigand, tLigand, nSpacer))
    
"""
	python3 stitchFF.py {{ wano["Monomer type for the polymeric backbone"] }} {{ wano["Crosslinker type for the polymer network"] }} {{ wano["Ligand type"] }}
	cp {{ wano["Monomer type for the polymeric backbone"] }}2.lmp Monomer.lmp
	cp FF.{{ wano["Monomer type for the polymeric backbone"] }}2 FF.Monomer
	cp {{ wano["Crosslinker type for the polymer network"] }}half.lmp Crosslinkerhalf.lmp
	cp FF.{{ wano["Crosslinker type for the polymer network"] }}half FF.Crosslinkerhalf
	cp mol.{{ wano["Crosslinker type for the polymer network"] }} mol.Crosslinker
	cp {{ wano["Ligand type"] }}.lmp Ligand.lmp
	cp FF.{{ wano["Ligand type"] }} FF.Ligand
"""
os.system("python3 stitchFF.py {} {} {}".format(tMonomer, tXLinker, tLigand))
os.system("cp {}2.lmp Monomer.lmp".format(tMonomer))
os.system("cp {}half.lmp Crosslinkerhalf.lmp".format(tXLinker))
os.system("cp {}.lmp Ligand.lmp".format(tLigand))

os.system("cp FF.{}2 FF.Monomer".format(tMonomer))
os.system("cp FF.{}half FF.Crosslinkerhalf".format(tXLinker))
os.system("cp FF.{} FF.Ligand".format(tLigand))

os.system("cp mol.{} mol.Crosslinker".format(tXLinker))
