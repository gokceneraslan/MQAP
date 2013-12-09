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

import joblib
import csv
import os
import pandas as pd

from encoders import MQAPmapper
from generateTraining import generateFeatures, generateTrainingSet

#needed for Joblib to reconstruct RandomForest object from the Joblib file
from sklearn.ensemble import RandomForestClassifier

def predict(rf, models, outputs=None, saveOutput=False, outputdir=None,
            templateFile=None):

    data = []
    assert isinstance(models, list), 'models parameter must be a list'

    assert outputs is None or isinstance(outputs, list), \
         ('outputs parameter must be a list or None')

    if outputs:
        assert len(models) == len(outputs), ('Number of models and number of'
                                            ' output file names do not match.')

    assert isinstance(rf, (str, RandomForestClassifier)), 'RF has invalid type'

    if isinstance(rf, str) and os.path.isfile(rf):
        # load trained RandomForest file in Joblib format
        rf = joblib.load(rf)

    # if templateFile argument is used, then netsurfp/psipred features will be
    # copied from this model to all model structures
    if templateFile:
        modelDFList = generateTrainingSet({templateFile: models}, 5,
                                          combineOutput=False)
        assert len(modelDFList) == len(models), "Some model PDBs cannot be parsed."
        modelDFList = [x.drop('ClassLabel', 1) for x in modelDFList]

    for i, modelFilename in enumerate(models):

        #cross-validation supplies the target model, we can use it to make
        #things faster by running psipred/netsurfp only once
        if not templateFile:
            modelDF = generateFeatures(modelFilename)
        else:
            modelDF = modelDFList[i]

        prediction = rf.predict(MQAPmapper.transform(modelDF))
        modelDF['ClassLabel'] = pd.Series(prediction)

        if saveOutput:
            if outputs: out = outputs[i]
            else:
                filename, ext = os.path.splitext(modelFilename)
                if ext != 'csv': out = filename + '.csv'
                else: out = filename + 'PREDICTION.csv'

            if outputdir:
                out = os.path.join(outputdir, os.path.basename(out))

            modelDF.to_csv(out, index=False, quoting=csv.QUOTE_NONNUMERIC)

        data.append(modelDF)

    return data
