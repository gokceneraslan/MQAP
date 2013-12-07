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
import psipred

def copyDataFromTarget(target, model, labels=('psipred', 'netsurfp')):
    data = {}
    for label in target.getDataLabels():
        if not label.startswith(tuple(labels)):
            continue
        data[label] = np.zeros(model.numAtoms(), dtype=target.getData(label).dtype)

    for tchain in target.iterChains():

        #try to find matching chain first
        if tchain.getChid() in [x.getChid() for x in model.iterChains()]:
            mchain = tchain.getChid()
        else:
            mchain = prody.matchChains(model, target[tchain.getChid()])[0][0]
            mchain = mchain.copy().getChids()[0]

        for tres in tchain.ca.copy().iterResidues():
            for label in data:
                try:
                    indices = model[mchain][tres.getResnum()].getIndices()
                    data[label][indices] = tres.ca.getData(label)
                except AttributeError: #target may have more residues than
                                       #the model
                    pass

    for label in data:
        if data[label].dtype.char == 'S':
            data[label][data[label] == ''] = '-'
        model.setData(label, data[label])

def generateFeatures(modelFilename):
    tempdir = tempfile.mkdtemp()

    datadict = {}
    modelPDB = prody.parsePDB(modelFilename)

    if not modelPDB:
        raise Exception('Model file %s cannot be parsed' % modelFilename)

    #if model has no chainID, let's assign one. That makes STRIDE parser
    #happy
    if np.unique(modelPDB.getChids()) == ' ':
        modelPDB.all.setChids('A')
        modelFilename = os.path.join(tempdir,
                                     os.path.basename(modelFilename))

        modelFilename = prody.writePDB(modelFilename,
                                       modelPDB,
                                       autoext=False)

    #run STRIDE
    prody.parseSTRIDE(prody.execSTRIDE(modelFilename, outputdir=tempdir),
                                       modelPDB)

    datadict['STRIDEarea'] = pd.Series(modelPDB.ca.getData('stride_area'))
    ss = pd.Series(modelPDB.ca.getSecstrs())
    ss[ss == ''] = '-' #empty strings cause trouble in csv load/save
    datadict['STRIDEss'] = ss

    #run DSSP
    prody.parseDSSP(prody.execDSSP(modelFilename, outputdir=tempdir), modelPDB)

    datadict['DSSPacc'] = pd.Series(modelPDB.ca.getData('dssp_acc'))
    ss = pd.Series(modelPDB.ca.getSecstrs())
    ss[ss == ''] = '-'
    datadict['DSSPss'] = ss

    #run NetSurfP
    netsurfp.parseNetSurfP(netsurfp.execNetSurfP(modelFilename,
                                                 outputdir=tempdir),
                                                 modelPDB)
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

    #run PSIPRED
    psipred.parsePSIPRED(psipred.execPSIPRED(modelFilename,
                                             outputdir=tempdir),
                                                 modelPDB)
    datadict['PSIPRED_ss'] = \
        pd.Series(modelPDB.ca.getData('psipred_ss'))
    datadict['PSIPRED_coilscore'] = \
        pd.Series(modelPDB.ca.getData('psipred_coilscore'))
    datadict['PSIPRED_helixscore'] = \
        pd.Series(modelPDB.ca.getData('psipred_helixscore'))
    datadict['PSIPRED_strandscore'] = \
        pd.Series(modelPDB.ca.getData('psipred_strandscore'))


    dataframe = pd.DataFrame(datadict)

    #remove temporary directory
    shutil.rmtree(tempdir, ignore_errors=True)

    return dataframe


# TODO: so many code duplication in generateFeatures & generateTrainingSet
# refoactor to reuse common code
def generateTrainingSet(inputdict, distance, output=None):

    devnull = open(os.devnull, 'w')
    subprocess.check_call('dssp --version', shell=True, stdout=devnull,
                          stderr=devnull)
    subprocess.check_call('stride -h', shell=True, stdout=devnull,
                          stderr=devnull)
    subprocess.check_call('netsurfp -h', shell=True, stdout=devnull,
                          stderr=devnull)
    subprocess.check_call('runpsipred', shell=True, stdout=devnull,
                          stderr=devnull)
    devnull.close()

    dataframe = pd.DataFrame()

    for target, models in inputdict.items():

        targetPDB = prody.parsePDB(target)
        assert distance, "Distance is not valid"

        tempdir = tempfile.mkdtemp()

        # we don't want to run NetSurfP and PSIPRED over and over again for all
        # model structures. we compute them for target structure and just reuse
        # on model structures
        netsurfp.parseNetSurfP(netsurfp.execNetSurfP(target, outputdir=tempdir),
                                                     targetPDB)

        psipred.parsePSIPRED(psipred.execPSIPRED(target, outputdir=tempdir),
                                                 targetPDB)

        for i, modelFilename in enumerate(models):
            datadict = {}
            modelPDB = prody.parsePDB(modelFilename)

            if not modelPDB:
                print('Model file %s cannot be parsed, skipping...' %
                      modelFilename)
                continue

            #if model has no chainID, let's assign one. That makes STRIDE parser
            #happy
            if np.unique(modelPDB.getChids()) == ' ':
                modelPDB.all.setChids('A')
                modelFilename = os.path.join(tempdir,
                                             os.path.basename(modelFilename))

                modelFilename = prody.writePDB(modelFilename,
                                               modelPDB,
                                               autoext=False)

            #superimpose model onto target structure
            match = prody.matchAlign(modelPDB, targetPDB, tarsel='calpha')
            mapmodel = match[1]
            maptarget = match[2]

            #and copy NetSurfP and PSIPRED data from target to model
            copyDataFromTarget(targetPDB, modelPDB)

            #run STRIDE
            prody.parseSTRIDE(prody.execSTRIDE(modelFilename, outputdir=tempdir),
                                               modelPDB)

            datadict['STRIDEarea'] = \
                pd.Series(modelPDB.ca.getData('stride_area')[mapmodel.getResindices()],
                          index=maptarget.getResindices())

            ss = pd.Series(mapmodel.getSecstrs(), index=maptarget.getResindices())
            ss[ss == ''] = '-' #empty strings cause trouble in csv load/save
            datadict['STRIDEss'] = ss

            #run DSSP
            prody.parseDSSP(prody.execDSSP(modelFilename, outputdir=tempdir), modelPDB)

            datadict['DSSPacc'] = \
                pd.Series(modelPDB.ca.getData('dssp_acc')[mapmodel.getResindices()],
                          index=maptarget.getResindices())

            ss = pd.Series(mapmodel.getSecstrs(), index=maptarget.getResindices())
            ss[ss == ''] = '-' #empty strings cause trouble in csv load/save
            datadict['DSSPss'] = ss


            #save NetSurfP data
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

            #save PSIPRED data
            datadict['PSIPRED_ss'] = \
                pd.Series(modelPDB.ca.getData('psipred_ss')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['PSIPRED_coilscore'] = \
                pd.Series(modelPDB.ca.getData('psipred_coilscore')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['PSIPRED_helixscore'] = \
                pd.Series(modelPDB.ca.getData('psipred_helixscore')[mapmodel.getResindices()],
                          index=maptarget.getResindices())
            datadict['PSIPRED_strandscore'] = \
                pd.Series(modelPDB.ca.getData('psipred_strandscore')[mapmodel.getResindices()],
                          index=maptarget.getResindices())

            #Compute class labels based on the distance argument
            datadict['ClassLabel'] = pd.Series((np.abs(
                prody.calcDistance(maptarget.copy(), mapmodel.copy())) < distance).astype(int),
                index=maptarget.getResindices())

            dataframe = pd.concat([dataframe, pd.DataFrame(datadict)])

        #remove temporary directory
        shutil.rmtree(tempdir, ignore_errors=True)

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

