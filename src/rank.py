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

import csv
import os

import pandas as pd

def rank(predictions, outputFilename):
    assert isinstance(predictions, list), 'Predictions argument must be a list'
    assert len(predictions) > 0, 'There must be at least one prediction files'

    scores = [] #list of tuples: (prediction, score)

    for predictionFile in predictions:
        prediction = pd.read_csv(predictionFile, quoting=csv.QUOTE_NONNUMERIC)
        predname = os.path.splitext(os.path.basename(predictionFile))[0]
        score = prediction.ClassLabel.values.mean()
        if predname in (x[0] for x in scores):
            predname += ":" + str(len([x for x in scores
                                         if x[0].startswith(predname+':')])+1)
        scores.append((predname, score))

    scores.sort(key=lambda x: x[1], reverse=True)
    df = pd.DataFrame([x[0] for x in scores], columns=['Predictions'])
    df['Scores'] = [x[1] for x in scores]

    df.to_csv(outputFilename, index=False, quoting=csv.QUOTE_NONNUMERIC)

