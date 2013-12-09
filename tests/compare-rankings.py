#!/usr/bin/env python

## Compare rankings
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

import pandas as pd
import numpy as np
import scipy.stats

parser = argparse.ArgumentParser()

parser.add_argument('rankings', help='Ranking files in CSV format',
                               nargs=2, metavar='RANKING')

args = parser.parse_args()

file1 = pd.read_csv(args.rankings[0])
file2 = pd.read_csv(args.rankings[1])

rank1 = file1['Predictions'].values
rank2 = file2['Predictions'].values

encoded_rank1 = np.array(range(len(rank1)))
encoded_rank2 = np.array([np.where(rank1 == x)[0][0] for x in rank2])

#print encoded_rank1
#print '============'
#print encoded_rank2

sp, _ = scipy.stats.spearmanr(encoded_rank1, encoded_rank2)
kd, _ = scipy.stats.kendalltau(encoded_rank1, encoded_rank2)

print "Spearman rank-order correlation coefficient: %f" % sp
print "Kendall's tau rank-order correlation coefficient: %f" % kd
