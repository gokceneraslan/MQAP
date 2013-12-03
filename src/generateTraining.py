## mqap
# Copyright (C) 2013 Gokcen Eraslan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import subprocess
import tempfile
import os
import shutil
import csv

import numpy as np
import pandas as pd
import prody

import netsurfp

def generateFeatures(models, target=None, output=None, distance=None):
    if target:
        targetPDB = prody.parsePDB(target)
        assert distance, "Distance is not valid"

    modelFilenames = models.split(',')

    if target:
        dataframe = pd.DataFrame(columns=['STRIDEarea', 'STRIDEss', 'DSSPacc',
                                          'DSSPss', 'ClassLabel'])
    else:
        dataframe = pd.DataFrame(columns=['STRIDEarea', 'STRIDEss', 'DSSPacc',
                                          'DSSPss'])

    for i, modelFilename in enumerate(modelFilenames):
        datadict = {}
        tempdir = tempfile.mkdtemp()
        modelPDB = prody.parsePDB(modelFilename)

        #if model has no chainID, let's assign one. That makes STRIDE parser
        #happy
        if np.unique(modelPDB.getChids()) == ' ':
            modelPDB.all.setChids('A')
            modelFilename = os.path.join(tempdir,
                                         os.path.basename(modelFilename))

            modelFilename = prody.writePDB(modelFilename,
                                           modelPDB,
                                           autoext=False)

        if target:
            match = prody.matchAlign(modelPDB, targetPDB, tarsel='calpha')
            mapmodel = match[1]
            maptarget = match[2]

        #run STRIDE
        prody.parseSTRIDE(prody.execSTRIDE(modelFilename, outputdir=tempdir),
                                           modelPDB)

        if target:
            datadict['STRIDEarea'] = \
                pd.Series(modelPDB.ca.getData('stride_area')[mapmodel.getResindices()],
                          index=maptarget.getResindices())

            datadict['STRIDEss'] = pd.Series(mapmodel.getSecstrs(),
                                             index=maptarget.getResindices())
        else:
            datadict['STRIDEarea'] = pd.Series(modelPDB.ca.getData('stride_area'))
            datadict['STRIDEss'] = pd.Series(modelPDB.ca.getSecstrs())

        #run DSSP
        prody.parseDSSP(prody.execDSSP(modelFilename, outputdir=tempdir), modelPDB)

        if target:
            datadict['DSSPacc'] = \
                pd.Series(modelPDB.ca.getData('dssp_acc')[mapmodel.getResindices()],
                          index=maptarget.getResindices())

            ss = pd.Series(mapmodel.getSecstrs(), index=maptarget.getResindices())
            ss[ss == ''] = '-' #empty strings cause trouble in csv load/save
            datadict['DSSPss'] = ss
        else:
            datadict['DSSPacc'] = pd.Series(modelPDB.ca.getData('dssp_acc'))
            ss = pd.Series(modelPDB.ca.getSecstrs())
            ss[ss == ''] = '-'
            datadict['DSSPss'] = ss


        #run NetSurfP
        netsurfp.parseNetSurfP(netsurfp.execNetSurfP(modelFilename,
                                                     outputdir=tempdir),
                                                     modelPDB)

        if target:
            datadict['NetSurfP_exp'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_exposure')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['NetSurfP_asa'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_asa')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['NetSurfP_rsa'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_rsa')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['NetSurfP_alpha'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_alphascore')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['NetSurfP_beta'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_betascore')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['NetSurfP_coil'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_coilscore')[mapmodel.getResindices()],
                          index=maptarget.getResindices())

        else:
            datadict['NetSurfP_exp'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_exposure'))
            datadict['NetSurfP_asa'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_asa'))
            datadict['NetSurfP_rsa'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_rsa'))
            datadict['NetSurfP_alpha'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_alphascore'))
            datadict['NetSurfP_beta'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_betascore'))
            datadict['NetSurfP_coil'] = \
                pd.Series(modelPDB.ca.getData('netsurfp_coilscore'))


        #remove temporary directory
        shutil.rmtree(tempdir, ignore_errors=True)

        if target:
            datadict['ClassLabel'] = pd.Series((np.abs(
                prody.calcDistance(maptarget.copy(), mapmodel.copy())) < distance).astype(int),
                index=maptarget.getResindices())

        dataframe = pd.concat([dataframe, pd.DataFrame(datadict)])

    if output:
        dataframe.to_csv(output, index=False, quoting=csv.QUOTE_NONNUMERIC)
        #dataframe.to_csv(output)

    print(dataframe)
    return dataframe


# calculate transformation matrix using DeepAlign tool
def calculateDeepAlignTransformation(fixed, moving):
    tempdir = tempfile.mkdtemp()
    tempalignfile = os.path.join(tempdir, 'deepalign')
    p = subprocess.Popen(['DeepAlign', '-t', moving,
                          '-q', fixed, '-u', '1',
                          '-o', tempalignfile], universal_newlines=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    #read translation/rotation matrix from .score file generated by deepalign
    with open(tempalignfile + '.score') as scorefile:
        matrix = scorefile.readlines()[-3:]

    shutil.rmtree(tempdir, ignore_errors=True)
    matrix = np.array([x.strip().split() for x in matrix]).astype('d')
    translation = matrix[:,-1]
    rotation = np.transpose(matrix[:, 0:3])
    return ((rotation, translation))


