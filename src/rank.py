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

    for prediction in predictions:
        if isinstance(prediction, str):
            label = os.path.splitext(os.path.basename(prediction))[0]
            prediction = pd.read_csv(prediction, quoting=csv.QUOTE_NONNUMERIC)
        else:
            assert isinstance(prediction, pd.DataFrame), ("Prediction has "
                                                          "invalid type.")
            label = prediction.label

        score = prediction.ClassLabel.values.mean()
        if label in (x[0] for x in scores):
            label += ":" + str(len([x for x in scores
                                         if x[0].startswith(label+':')])+1)
        scores.append((label, score))

    df = pd.DataFrame([x[0] for x in scores], columns=['Predictions'])
    df['Scores'] = [x[1] for x in scores]
    df.sort(columns=['Scores', 'Predictions'], ascending=[0, 1], inplace=True)

    df.to_csv(outputFilename, index=False, quoting=csv.QUOTE_NONNUMERIC)

