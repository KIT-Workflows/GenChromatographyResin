import LammpsDataDict as LD

ld=LD.LammpsData()
ld.readFile("fREF_polymer")
#ld.masses=ld.masses[:-2] #disregard initiators
del ld.masses[26]
del ld.masses[27]
ld.writeFile("fRef_pos", reorder=True)
