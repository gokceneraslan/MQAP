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

import cPickle as pickle
import csv
import os
import pandas as pd

from encoders import MQAPmapper
from generateTraining import generateFeatures

#needed for pickle to reconstruct RandomForest object from the pickle file
from sklearn.ensemble import RandomForestClassifier

def predict(rfFilename, models, outputs=None, saveOutput=False):
    data = []
    assert isinstance(models, list), 'models parameter must be a list'
    assert isinstance(outputs, (list, None)), ('outputs parameter must be '
                                               'a list or None')

    if isinstance(outputs, list):
        assert len(models) == len(outputs), ('Number of models and number of'
                                            ' output file names do not match.')

    # load trained RandomForest file in pickle format
    with open(rfFilename) as rfFile:
        rf = pickle.load(rfFile)

    for i, modelFilename in enumerate(models):
        newdata = generateFeatures(modelFilename)
        prediction = rf.predict(MQAPmapper.transform(newdata))
        newdata['ClassLabel'] = pd.Series(prediction)

        if saveOutput:
            if outputs: out = outputs[i]
            else:
                filename, ext = os.path.splitext(modelFilename)
                if ext != 'csv': out = filename + '.csv'
                else: out = filename + 'PREDICTION.csv'

            newdata.to_csv(out, index=False, quoting=csv.QUOTE_NONNUMERIC)

        data.append(newdata)

    return data
