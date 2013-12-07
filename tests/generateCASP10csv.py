#!/usr/bin/env python

import csv
import os
import sys

caspdir = sys.argv[1]

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
            files.remove(os.path.basename(targetfile))
            files.sort()
            writer.writerow([targetfile] + [os.path.join(root, f) for f in files])
