#!/usr/bin/env python

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

import argparse
import subprocess
import csv
import os

import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--models', help='Model structures in '
                              'PDB format.', nargs='+', required=True)
parser.add_argument('-o', '--output', help='Output file in CSV format',
                              required=True)

args = parser.parse_args()

p = subprocess.Popen(['dfire2'] + args.models, universal_newlines=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)

out, err = p.communicate()
data = [(os.path.basename(x.split()[0]), float(x.split()[1])) for x in out.split('\n') if x]

df = pd.DataFrame([x[0] for x in data], columns=['Predictions'])
df['Scores'] = [x[1] for x in data]
df.sort(columns=['Scores'], inplace=True, ascending=True)

df.to_csv(args.output, index=False, quoting=csv.QUOTE_NONNUMERIC)
