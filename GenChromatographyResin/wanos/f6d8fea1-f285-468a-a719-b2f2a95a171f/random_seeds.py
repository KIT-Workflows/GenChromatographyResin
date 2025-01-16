import random
import os
from datetime import datetime
import sys

rnd=sys.argv[1]
if rnd=="True":
    rnd=True
else:
    rnd=False

print(sys.argv)
boxx, boxy, boxz, numMonomer, numCrosslinker, numLigand, numPorogen=sys.argv[2:]

if rnd:
    random.seed(datetime.now().timestamp())
    rngs=random.sample(range(1, 10000), 5)
else:
    defaultseeds=[1234, 2345, 3456, 4567, 5678]
    rngs=defaultseeds
    
os.system("lmp_new_colvars -in in.createbox -var dimX {boxx} -var dimY {boxy} -var dimZ {boxz} -var numHEMA {numMonomer} -var numEGDMA {numCrosslinker} -var numTRP {numLigand} -var numPoro {numPorogen} -var rng1 {rng1} -var rng2 {rng2} -var rng3 {rng3} -var rng4 {rng4} -var rng5 {rng5} ".format(boxx=boxx, boxy=boxy, boxz=boxz, numMonomer=numMonomer, numCrosslinker=numCrosslinker, numLigand=numLigand, numPorogen=numPorogen, rng1=rngs[0], rng2=rngs[1], rng3=rngs[2], rng4=rngs[3], rng5=rngs[4]))
