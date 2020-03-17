#    README file for further information
#    Copyright (C) 2020  Federica Angiolari

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import coll
import sm


def main(mod, n_angol, program, scr_dir, num_t, num_p, n, path_program):
    num_atom = num_t + num_p
    if program == 'GAMESS':
        for i in range(n):
            print("Dynamic: {}".format(i))
            os.system("{} geometry_{}.inp > output_gamess".format(path_program, i))
        os.system("mv {}geometry*trj .".format(scr_dir))
    else:
        for i in range(n):
            os.system("mv geometry_{}.gen geometry.gen".format(i))
            print("Dynamic: {}".format(i))
            os.system("{} > output_dftb".format(path_program))
            os.system("mv NVE.xyz NVE_{}.xyz".format(i))
            os.system("mv md.out md_{}.out".format(i))
    if mod == 'coll':
        coll.main(n, num_atom, num_p, program, n_angol)
    else:
        sm.main(num_t, n)
