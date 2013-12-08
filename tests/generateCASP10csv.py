#!/usr/bin/env python

import csv
import os
import sys

from prody import *

caspdir = sys.argv[1]
confProDy(verbosity='none')

with open(sys.argv[2], 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='"')
    writer.writerow(['Target File', 'Models'])

    for root, dirs, files in os.walk(caspdir):
        if root == caspdir:
            targets = dirs
        else:
            target = os.path.basename(root)
            targetfile = os.path.join(root, target + '.pdb')

            assert os.path.isfile(targetfile), ('Target %s cannot be found' %
                                                targetfile)

            targetPDB = parsePDB(targetfile)
            if targetPDB.numResidues() != targetPDB.ca.numAtoms():
                   print "Target %s is invalid" % targetfile
                   print "Numresidues: %d vs. numCA:%d" % \
                   (targetPDB.numResidues(), targetPDB.ca.numAtoms())
                   continue

            files.remove(os.path.basename(targetfile))
            files.sort()
            files = [os.path.join(root, f) for f in files]
            filesFiltered = []
            for f in files:
                fPDB = parsePDB(f)
                if not fPDB:
                    print 'Invalid PDB %s, skipping...' % f
                    continue

                if fPDB.numResidues() == fPDB.ca.numAtoms():
                    filesFiltered.append(f)
                else:
                    print 'Invalied number of CA in %s, skipping...' % f

            writer.writerow([targetfile] + filesFiltered)
