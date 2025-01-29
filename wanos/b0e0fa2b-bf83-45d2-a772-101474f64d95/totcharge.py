#! /usr/bin/python3

####################################
# Script to insert water molecules #
# into the zeolite structure.      #
####################################

#import argparse
import LammpsData as LD

########################
## Start script here! ##
########################

## Parse command line arguments
## ----------------------------

#parser = argparse.ArgumentParser(description=
#'''Insert cation/anion pairs into MOF structure'''
#)

#parser.add_argument("input", type=str, help="target file to re-ID")
#args = parser.parse_args()

def totcharge(inFile):

    ld = LD.LammpsData()
    #ld.readFile(args.input)
    ld.readFile(inFile)

    totch=[]
    totweight=[]
    molmax = max([ a[0] for a in ld.atoms ])
    for i in range(molmax):
	    totchtemp=0.0
	    weighttemp=0.0
	    for atom in ld.atoms:
		    if atom[0]==i+1:		
			    totchtemp+=atom[2]
			    weighttemp+=ld.masses[atom[1]-1]
			    print(weighttemp)
	    totch.append(totchtemp)
	    totweight.append(weighttemp)
    tot=sum(totch)
    totw=sum(totweight)
    print("====================================")
    #print("CHARGES and WEIGHTS in "+args.input)
    print("------------------------------------")
    for item in range(len(totch)):
	    print("molecule:\t"+str(item+1)+"\t"+str(totch[item])+"\t"+str(totweight[item]))
    print("------------------------------------")
    print("total charge:\t"+str(tot))
    print("total weight:\t"+str(totw))
    print("====================================")
    # Write out new data file

    return tot
