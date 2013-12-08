## MQAP
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

import os
import csv

from sklearn.cross_validation import KFold

from generateTraining import generateTrainingSet
from train import train
from predict import predict

def crossvalidate(input, fold, output, distance, ntrees, trainingdir,
                  randomforestdir, predictiondir):
    cases = []

    with open(input) as inputfile:
        reader = csv.reader(inputfile, delimiter=',', quotechar='"')
        for row in reader:
            cases.append((row[0], row[1:]))

    iteration = 1
    kfold = KFold(len(cases), n_folds=fold, shuffle=True)
    for train_indices, test_indices in kfold:
        generateTrainingInput = {}

        for train_index in train_indices:
            generateTrainingInput[cases[train_index][0]] = \
                                  cases[train_index][1]

        trainingset = generateTrainingSet(generateTrainingInput, distance,
                                          output=trainingdir)
        rfoutput = None
        if randomforestdir:
            rfoutput = os.path.join(randomforestdir, 'cv-fold%d-rf.pickle' %
                                          iteration)

        rf = train(trainingset, ntrees, rfoutput)
        testmodels = []
        for test_index in test_indices:
            testmodels += cases[test_index][1]

        predictions = predict(rf, testmodels, saveOutput=predictiondir,
                              outputdir=predictiondir)
        iteration += 1
