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

from sklearn.preprocessing import LabelEncoder
from sklearn_pandas import DataFrameMapper

DSSPssEncoder = LabelEncoder()
DSSPssEncoder.fit(['H', 'B', 'E', 'G', 'I', 'T', 'S', '-'])

STRIDEssEncoder = LabelEncoder()
STRIDEssEncoder.fit(['H', 'B', 'E', 'G', 'I', 'T', 'S', 'C', 'b'])

MQAPmapper = DataFrameMapper([
                              ('DSSPacc', None),
                              ('DSSPss', DSSPssEncoder),
                              ('STRIDEarea', None),
                              ('STRIDEss', STRIDEssEncoder)])

