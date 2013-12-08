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

import cPickle as pickle
import csv

import pandas as pd
from sklearn.ensemble import RandomForestClassifier

from encoders import MQAPmapper

def train(training, ntrees, outputFile=None):
    if isinstance(training, pd.DataFrame):
        data = training
    else:
        data = pd.read_csv(training, quoting=csv.QUOTE_NONNUMERIC)

    columns = list(data.columns.values)
    assert len(data) > 0, 'Training set must have at least one instance!'

    labels = data['ClassLabel'].values

    columns.remove('ClassLabel')
    data = data[columns]

    rf = RandomForestClassifier(n_estimators=ntrees, compute_importances=True,
                                verbose=1, n_jobs=-1)

    rf.fit(MQAPmapper.transform(data), labels)

    if outputFile:
        with open(outputFile, 'wb') as output:
            pickle.dump(rf, output, pickle.HIGHEST_PROTOCOL)

    return rf
