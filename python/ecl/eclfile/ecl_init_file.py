#  Copyright (C) 2015  Equinor ASA, Norway.
#
#  The file 'ecl_init_file.py' is part of ERT - Ensemble based Reservoir Tool.
#
#  ERT is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ERT is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
#  for more details.

from ecl import EclFileEnum, EclFileFlagEnum
from ecl.eclfile import Ecl3DFile, EclFile

ECL_FILE_DEFAULT = EclFileFlagEnum.ECL_FILE_DEFAULT


class EclInitFile(Ecl3DFile):
    def __init__(self, grid, filename, flags=ECL_FILE_DEFAULT):
        file_type, report_step, fmt_file = EclFile.get_file_type(filename)
        if file_type == EclFileEnum.ECL_INIT_FILE:
            super(EclInitFile, self).__init__(grid, filename, flags)
        else:
            err = 'The input filename "%s" does not correspond to an init file.'
            err += '  Please follow the Eclipse naming conventions.'
            raise ValueError(err % filename)
