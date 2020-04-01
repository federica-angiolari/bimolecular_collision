#    Create initial geometries and input for dynamic
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

import numpy as np
import math
import system  # Hide me if you want just the geometries and input for dynamic :), remember to hide also line 472
import sys

# Constants
avogadro = 6.022e23  # mol^-1
kb = 1.380649e-23  # J/K
ms_to_bohrps = 0.018897  # m/s -> bohr/ps 
kcal_to_joule = 4184  # kcal -> J
ev_to_joule = 1.602176e-19  # eV -> J 
mass = {"I": 126.9044719e-3, "F": 18.9984032e-3, "C": 12.000e-3, "N": 14.003074e-3, "H": 1.007825032e-3, "O": 15.994914619e-3}  # Kg/mol
atomic_number = {"H": 1.0, "O": 8.0, "C": 6.0}
# Value for 3ob-3-1 > input_dftb
hubbard_der = {"Br":-0.0573, "C":-0.1492, "Ca":-0.0340, "Cl":-0.0697, "F":-0.1623, "H":-0.1857, "I":-0.0433, "K":-0.0339, "Mg":-0.02, "N":-0.1535, "Na":-0.0454, "O":-0.1575, "P":-0.14, "S":-0.11, "Zn":-0.03}
ang_m = {"Br":"d", "C":"p", "Ca":"p", "Cl":"d", "F":"p", "H":"s", "I":"d", "K":"p", "Mg":"p", "N":"p", "Na":"p", "O":"p", "P":"d", "S":"d", "Zn":"d"}
spin_const = {"H":-0.072, "C":-0.023, "N":-0.026, "O": -0.028 }



def mass_array(label):
    m = []
    if len(label) == 1:
        m.append(mass[str(label)])
    else:
        for i in range(len(label)):
            m.append(mass[label[i]])
    return m


def inertia_matrix(coord, label):

    i_m = list([0]*3)
    for i in range(3):
        i_m[i] = [0]*3

    for i in range(label.size):
        i_m[0][0] += mass[label[i]]*((coord[1][i]**2) + (coord[2][i]**2))
        i_m[1][1] += mass[label[i]]*((coord[0][i]**2) + (coord[2][i]**2))
        i_m[2][2] += mass[label[i]]*((coord[0][i]**2) + (coord[1][i]**2))
        i_m[0][1] += mass[label[i]]*(-coord[0][i]*coord[1][i])
        i_m[0][2] += mass[label[i]]*(-coord[0][i]*coord[2][i])
        i_m[1][2] += mass[label[i]]*(-coord[1][i]*coord[2][i])
    i_m[1][0] = i_m[0][1]
    i_m[2][0] = i_m[0][2]
    i_m[2][1] = i_m[1][2]
    return i_m


def rotation_single_geom_in_inertia_axis(m_inertia, coord, num):

    [eigenvalues, eigenvectors] = np.linalg.eigh(np.matrix(m_inertia))
    transpose_matrix = np.transpose(np.asarray(eigenvectors))
    coord_new = empty_matrix(len(coord[0]), len(coord))
    for i in range(num):
        coord_new[0][i] = transpose_matrix[0][0] * coord[0][i] + transpose_matrix[0][1] * coord[1][i] + transpose_matrix[0][2] * coord[2][i]
        coord_new[1][i] = transpose_matrix[1][0] * coord[0][i] + transpose_matrix[1][1] * coord[1][i] + transpose_matrix[1][2] * coord[2][i]
        coord_new[2][i] = transpose_matrix[2][0] * coord[0][i] + transpose_matrix[2][1] * coord[1][i] + transpose_matrix[2][2] * coord[2][i]
    return coord_new


def empty_matrix(rows, col):
    matrix = [float(0)]*col
    for i in range(col):
        matrix[i] = [float(0)]*rows
    return np.array(matrix)


def shift_cm(matrix, label):
    m = mass_array(label)
    tmpx, tmpy, tmpz = [], [], []
    if len(label) == 1:
        matrix_new = [[0], [0], [0]]
    else:
        for i in range(len(label)):
            tmpx.append(matrix[0][i]*mass[label[i]])
            tmpy.append(matrix[1][i]*mass[label[i]])
            tmpz.append(matrix[2][i]*mass[label[i]])
        cmx = sum(tmpx)/sum(m)
        cmy = sum(tmpy)/sum(m)
        cmz = sum(tmpz)/sum(m)
        matrix_new = empty_matrix(len(matrix[0]), len(matrix))
        for i in range(len(label)):
            matrix_new[0][i] = matrix[0][i] - cmx
            matrix_new[1][i] = matrix[1][i] - cmy
            matrix_new[2][i] = matrix[2][i] - cmz
    return matrix_new


# Velocities
def get_velocities(d, program, bullet_label, target_label, temperature_target, temperature_bullett, energy_input):
    if bullet_label.size== 1:
        bullet_label = np.str(bullet_label)
    if program == 'GAMESS':
        energy_input = energy_input * kcal_to_joule
    else:
        energy_input = energy_input * ev_to_joule * avogadro 
    target_mass = mass_array(target_label)
    total_target_mass = sum(target_mass)
    bullett_mass = mass_array(bullet_label)
    total_bullet_mass = sum(bullett_mass)
    red_mass = (total_bullet_mass * total_target_mass)/(total_bullet_mass + total_target_mass)

    v_total_coll_bullett = empty_matrix(len(bullet_label), 3)
    v_total_vib_bullett = empty_matrix(len(bullet_label), 3)
    v_total_coll_target = empty_matrix(len(target_label), 3)
    v_total_vib_target = empty_matrix(len(target_label), 3)

    if d < 0:
        for i in range(0, len(bullet_label)):
            v_total_coll_bullett[0][i] = np.sqrt(2*energy_input/red_mass)
    else:
        for i in range(0, len(bullet_label)):
            v_total_coll_bullett[0][i] = -np.sqrt(2*energy_input/red_mass)
    
    time_need_for_coll = float((abs(d*1e-10))/(abs(v_total_coll_bullett[0][0])))  # second 

    if len(np.str(bullet_label)) > 1:
        for i in range(0, len(bullet_label)):
            v_total_vib_bullett[0][i] = np.random.normal(0, np.sqrt(kb*temperature_bullett/(mass[bullet_label[i]]/avogadro)))
            v_total_vib_bullett[1][i] = np.random.normal(0, np.sqrt(kb*temperature_bullett/(mass[bullet_label[i]]/avogadro)))
            v_total_vib_bullett[2][i] = np.random.normal(0, np.sqrt(kb*temperature_bullett/(mass[bullet_label[i]]/avogadro)))
        v_total_vib_bullett_cm = shift_cm(v_total_vib_bullett, bullet_label)
    else:
        v_total_vib_bullett_cm = empty_matrix(len(bullet_label), 3)
    
    for i in range(0, len(target_label)):
        v_total_vib_target[0][i] = np.random.normal(0, np.sqrt(kb*temperature_target/(mass[target_label[i]]/avogadro)))
        v_total_vib_target[1][i] = np.random.normal(0, np.sqrt(kb*temperature_target/(mass[target_label[i]]/avogadro)))
        v_total_vib_target[2][i] = np.random.normal(0, np.sqrt(kb*temperature_target/(mass[target_label[i]]/avogadro)))
    v_total_vib_target_cm = shift_cm(v_total_vib_target, target_label)
    
    # sum vel
    v_total_fin = [[], [], []] 
    for i in range(0, len(bullet_label)):
        v_total_fin[0].append(v_total_coll_bullett[0][i] + v_total_vib_bullett_cm[0][i])
        v_total_fin[1].append(v_total_coll_bullett[1][i] + v_total_vib_bullett_cm[1][i])
        v_total_fin[2].append(v_total_coll_bullett[2][i] + v_total_vib_bullett_cm[2][i])
    for i in range(0, len(target_label)):
        v_total_fin[0].append(v_total_coll_target[0][i] + v_total_vib_target_cm[0][i])
        v_total_fin[1].append(v_total_coll_target[1][i] + v_total_vib_target_cm[1][i])
        v_total_fin[2].append(v_total_coll_target[2][i] + v_total_vib_target_cm[2][i])
    
    label_total = []
    for i in range(len(bullet_label)):
        label_total.append(bullet_label[i])
    for i in range(0, len(target_label)):
        label_total.append(target_label[i])
    v_total_fin_cm = shift_cm(v_total_fin, label_total)

    return v_total_fin_cm, time_need_for_coll


#def angles(num):
#   theta = [[0,0,0],[]..]   # [[angle_x, angle_y, angle_z],[...],[...]..]
#   return theta


def angles(num):
    theta = empty_matrix(3, num)
    k = 0
    while k < num:
        j = 0
        while j < 2:
            theta[k][j] = np.random.uniform(-math.pi, math.pi)
            j += 1
            theta[k][j] = np.random.uniform(0, math.pi)
            j += 1
            theta[k][j] = np.random.uniform(-math.pi, math.pi)
        k += 1
    return theta


def rotation(coord, theta):
    num_atomi = len(coord[0])
    coord_1 = empty_matrix(num_atomi, 3)
    coord_2 = empty_matrix(num_atomi, 3)
    coord_3 = empty_matrix(num_atomi, 3)
    for i in range(num_atomi):
        coord_1[1][i] = coord[1][i] * math.cos(theta[0]) - coord[2][i] * math.sin(theta[0])
        coord_1[2][i] = coord[1][i] * math.sin(theta[0]) + coord[2][i] * math.cos(theta[0])
        coord_1[0][i] = coord[0][i]
    for i in range(num_atomi):
        coord_2[2][i] = coord_1[2][i] * math.cos(theta[1]) - coord_1[0][i] * math.sin(theta[1])
        coord_2[0][i] = coord_1[2][i] * math.sin(theta[1]) + coord_1[0][i] * math.cos(theta[1])
        coord_2[1][i] = coord_1[1][i]
    for i in range(num_atomi):
        coord_3[0][i] = coord_2[0][i] * math.cos(theta[2]) - coord_2[1][i] * math.sin(theta[2])
        coord_3[1][i] = coord_2[0][i] * math.sin(theta[2]) + coord_2[1][i] * math.cos(theta[2])
        coord_3[2][i] = coord_2[2][i]
    return coord_3


def geometries_output(icharge, step, count_traj, label_tot, b, distance, coord_targ, coord_bullet, theta, program_used, scf_type, timestep, mult, v):
    coord_final = empty_matrix(len(label_tot), 3)
    counting_b = 0
    while counting_b < len(b):
        
        for i in range(0, len(coord_bullet[0])):
            coord_final[0][i] = coord_bullet[0][i] + distance
            coord_final[1][i] = coord_bullet[1][i] + b[counting_b]
            coord_final[2][i] = coord_bullet[2][i]
        k = len(coord_bullet[0])
        for i in range(0, len(coord_targ[0])):
            coord_final[0][k] = coord_targ[0][i]
            coord_final[1][k] = coord_targ[1][i]
            coord_final[2][k] = coord_targ[2][i]
            k += 1
        coord_final_cm = shift_cm(coord_final, label_tot)
        print_xyz(count_traj, coord_final_cm, label_tot)
        if program_used == 'DFTB+':
            print_gen(count_traj, coord_final_cm, label_tot)
        elif program_used == 'GAMESS':
            print_gamess_input(icharge, step, count_traj, distance, scf_type, timestep, coord_final_cm, label_tot, mult, v)
        print_riepilogo(count_traj, b[counting_b], theta)
        count_traj += 1
        counting_b += 1
    return count_traj


def print_gamess_input(icharge, step_total, n, dist, scf, dt, coord, label, molt, vel_ms):
    str1 = " $BASIS GBASIS=N311 NGAUSS=6 NDFUNC=1 NPFUNC=1 DIFFSP=.TRUE. DIFFS=.TRUE. $END"
    str2 = " $CONTRL MAXIT=200 SCFTYP={} RUNTYP=MD DFTTYP=B3LYP ICHARG={} MULT={} $END\n $MD READ=.TRUE. MBT=.FALSE. MBR=.FALSE. TTOTAL=0.000000E+00\n MDINT= VVERLET     DT={}  NVTNH= 0  NSTEPS=   {}\n RSTEMP=.F. DTEMP=   100.00 LEVERY=   10000\n RSRAND=.F. NRAND=  1000 NVTOFF=  0 JEVERY=    10\n PROD=.F.   KEVERY=     1 DELR=   0.020".format(scf, icharge, molt, dt, step_total)

    v_bohr = empty_matrix(len(coord[0]), 3)
    for i in range(len(coord[0])):
        v_bohr[0][i] = vel_ms[0][i] * ms_to_bohrps
        v_bohr[1][i] = vel_ms[1][i] * ms_to_bohrps
        v_bohr[2][i] = vel_ms[2][i] * ms_to_bohrps
    with open("geometry_{}.inp".format(n), "w") as f:
        f.write("{}\n".format(str1))
        f.write("{}\n".format(str2))
        f.write(" TVELQM(1)=     ! QM ATOM TRANS. VELOCITIES (BOHR/PS) !\n")
        for i in range(len(coord[0])):
            f.write(" {} {} {}\n".format(v_bohr[0][i], v_bohr[1][i], v_bohr[2][i]))
        f.write(" $END\n")
        f.write(" $DATA\n Title\n C1")
        f.write("\n")
        for i in range(len(coord[0])):
            f.write(" {} {} {} {} {}\n".format(label[i], atomic_number[label[i]], coord[0][i], coord[1][i], coord[2][i]))
        f.write(" $END")


def print_riepilogo(n, b, angle):
    with open("summ_traj", "a") as f:
        f.write("Traj= {} ,b= {} ,angle= {} {} {}\n".format(n, round(b, 2), round(angle[0], 2), round(angle[1], 2), round(angle[2], 2)))


def print_xyz(n, coord, label):
    with open("geometry_{}.xyz".format(n), "w") as f:
        f.write("{}\n\n".format(label.size))
        for i in range(0, label.size):
            f.write("{}\t{}\t{}\t{}\n".format(label[i], coord[0][i], coord[1][i], coord[2][i]))


def print_gen(n, coord, label):
    sel = list(dict.fromkeys(label))
    col1 = []
    for i in range(label.size):
        col1.append(i+1)
    col2 = []
    for i in range(label.size):
        for j in range(len(sel)):
            if label[i] == sel[j]:
                col2.append(j+1)
    with open("geometry_{}.gen".format(n), "w") as f:
        f.write(" {} C\n".format(label.size))
        for i in range(len(sel)):
            f.write(" {} ".format(sel[i]))
        f.write("\n")
        for i in range(len(col1)):
            f.write(" {}\t{}\t{}\t{}\t{}\n".format(col1[i], col2[i], coord[0][i], coord[1][i], coord[2][i]))


def input_dftb(dt, path, label, step, temperature_conv, spin_spaiati, charge):
    sel_atoms = list(dict.fromkeys(label))
    if "Ca" in sel_atoms:
        atom_ca = ""
    else:
        atom_ca = "#"
    if "Cl" in sel_atoms:
        atom_cl = ""
    else:
        atom_cl = "#"
    if "O" in sel_atoms:
        atom_o = ""
    else:
        atom_o = "#"
    if "N" in sel_atoms:
        atom_n = ""
    else:
        atom_n = "#"
    if "C" in sel_atoms:
        atom_c = ""
    else:
        atom_c = "#"
    if "H" in sel_atoms:
        atom_h = ""
    else:
        atom_h = "#"
    if "Br" in sel_atoms:
        atom_br = ""
    else:
        atom_br = '#'
    f = open("dftb_in.hsd", "w+")
    f.write("Geometry = GenFormat {\n <<< \"geometry.gen\"\n}\n")
    f.write("Driver = VelocityVerlet {\n MovedAtoms=1:-1\n OutputPrefix=\"NVE\"\n Steps=%d\n Timestep[fs]=%f\n KeepStationary=yes\n MDRestartFrequency=10\n Thermostat=None {}\n Velocities [m/s] ={\n  <<<\"VELOC.DAT\"\n}}" % (step, dt))
    f.write("Hamiltonian = DFTB {\n  SCC = Yes\n  SCCTolerance=1.0e-7\n  Charge = %d\n" % charge)
    if spin_spaiati != 0:
        f.write("SpinPolarisation= Colinear{\n UnpairedElectrons = %d\n}\nSpinConstants ={\n %s O=%f\n %sH =%f\n%sN=%f\n%sC=%f\n}" %(spin_spaiati, atom_o, spin_const["O"], atom_h, spin_const["H"], atom_n, spin_const["N"], atom_c, spin_const["C"]))
    f.write("\n  SlaterKosterFiles = Type2FileNames {\n    Prefix = \"%s\"\n    Separator = \"-\"\n    Suffix = \".skf\" \n}" % path)
    f.write("MaxAngularMomentum{")
    f.write("\n%s O=\"%s\"\n%s H=\"%s\"\n%s C=\"%s\"\n%s N=\"%s\"\n%s Cl=\"%s\"\n%s Br = \"%s\"\n}" % (atom_o, ang_m["O"], atom_h, ang_m["H"], atom_c, ang_m["C"], atom_n, ang_m["N"], atom_cl, ang_m["Cl"], atom_br, ang_m["Br"]))
    f.write("Filling=Fermi{\n Temperature[Kelvin]=%d\n}" % temperature_conv)
    f.write("Dispersion = DftD3 {\nDamping = BeckeJohnson {}\n}\nThirdOrderFull = Yes\n")
    f.write("HubbardDerivs {\n%sO= %f\n%sH = %f\n %s N = %f\n %s C = %f\n %s Br = %f\n %s Cl = %f\n %sCa = %f\n}\n" % (atom_o, hubbard_der["O"], atom_h, hubbard_der["H"], atom_n, hubbard_der["N"],  atom_c, hubbard_der["C"], atom_br, hubbard_der["Br"], atom_cl, hubbard_der["Cl"], atom_ca, hubbard_der["Ca"]))
    f.write("HCorrection=Damping{\nExponent=4.00}\n}\nOptions{}\nAnalysis{}\nParserOptions{\n ParserVersion=5\n}")
    f.close()


def check_correct_insert():
        print("ERROR")
        print("Correct use:\npython[vers>=3.0] start.py program mod\nprogram = DFTB+ / GAMESS\nmod = coll / ms\nExample:\npython3.6 start.py DFTB+ coll\n")
        sys.exit()

        
def main():

    # Input
    if len(sys.argv) != 3:
        check_correct_insert()
    if str(sys.argv[1])!='DFTB+' and str(sys.argv[1])!='GAMESS':
        check_correct_insert()
    if str(sys.argv[2])!='coll' and str(sys.argv[2])!='ms':
        check_correct_insert()
    program = str(sys.argv[1])  # GAMESS/DFTB+
    mod = str(sys.argv[2])  # coll/ms

    # Get parameters
    if program == 'GAMESS':
        parameters_list = np.loadtxt("start_input_gamess", dtype=tuple, usecols=[1])
        mod = 'coll'
    elif program == 'DFTB+':
        parameters_list = np.loadtxt("start_input_dftbplus", dtype=tuple, usecols=[1])
    else:
        parameters_list = 0
        print("error")
    name_target = str(parameters_list[0])  # Name file target
    name_bullett = str(parameters_list[1])  # Name file bullet
    d = float(parameters_list[2])  # Distance [Angstrom]
    limit = float(parameters_list[3])  # B_max [Angstrom]
    delta_b = float(parameters_list[4])  # Delta_b[Angstrom]
    numero_angoli = int(parameters_list[5])  # Angle's number (for target rotation)
    temperature_t = float(parameters_list[6])  # Roto-vib temperature for target
    temperature_p = float(parameters_list[7])  # " " " for bullet
    energy_coll = float(parameters_list[8])  # Collisional energy (in CM frame) [Kcal/mol]
    dt = float(parameters_list[9])  # Timestep for dynamic, in fs for dftbplus and s in gamess
    step_to_add = int(parameters_list[10])
    path_program = str(parameters_list[11])
    if program == "DFTB+":  # Parameters for dynamic (external programs)
        path = str(parameters_list[12])
        charge = int(parameters_list[13])
        temp_conv = float(parameters_list[14])
        spin_spaiati = int(parameters_list[15])
        r_oruhf = 0
        mult = 0
        scr_dir = 0
    else:
        dt = float(dt)
        r_oruhf = str(parameters_list[12])
        mult = int(parameters_list[13])
        icharge = int(parameters_list[14])
        scr_dir = str(parameters_list[15])
        path = 0
        charge = 0
        temp_conv = 0
        spin_spaiati = 0

    # Variables
    # Label atoms
    atoms_b = np.loadtxt(name_bullett, dtype=np.str, skiprows=2, usecols=[0])
    atoms_t = np.loadtxt(name_target, dtype=np.str, skiprows=2, usecols=[0])
    atoms_tot = []
    if atoms_b.size > 1:
        for i in range(atoms_b.size):
            atoms_tot.append(atoms_b[i])
    else:
        atoms_tot.append(str(atoms_b))
    for i in range(atoms_t.size):
        atoms_tot.append(atoms_t[i])
    atoms_tot = np.array(atoms_tot)

    # Initial coordinate
    coord_target = (empty_matrix(atoms_t.size, 3))
    coord_bull = (empty_matrix(atoms_b.size, 3))
    count = 1
    for i in range(0, int(3)):
        coord_target[i] = np.loadtxt(name_target, dtype=float, skiprows=2, usecols=[count])  # Matrix_coord = [[x1, x2..xn], [y1,y2..yn], [z1,z2..zn]]
        coord_bull[i] = np.loadtxt(name_bullett, dtype=float, skiprows=2, usecols=[count])
        count += 1

    # VELOCITIES
    velocities_m_s, total_time_coll = get_velocities(d, program, atoms_b, atoms_t, temperature_p, temperature_t, energy_coll)
    if program == 'GAMESS':   # dt = second
        step = int(total_time_coll/dt) + step_to_add
    elif program == 'DFTB+':  # dt = fs > total_time_coll*1e15
        step = int(total_time_coll/(dt*1e-15)) + step_to_add
    # GEOMETRY
    # Shift in CM system for each molecule
    if atoms_b.size > 1:
        coord_bull_cm = shift_cm(coord_bull, atoms_b)
    else:
        coord_bull_cm = [[0], [0], [0]]
    coord_target_cm = shift_cm(coord_target, atoms_t)
    # Rotation in inertia axis (valid only if atoms.size > 1)
    if atoms_b.size > 1:
        matrix_inertia_b = inertia_matrix(coord_bull_cm, atoms_b)
        coord_rotate_b = rotation_single_geom_in_inertia_axis(matrix_inertia_b, coord_bull_cm, atoms_b.size)
    else:
        coord_rotate_b = [[0], [0], [0]]
    matrix_inertia_t = inertia_matrix(coord_target_cm, atoms_t)
    coord_rotate_t = rotation_single_geom_in_inertia_axis(matrix_inertia_t, coord_target_cm, atoms_t.size)
    b = []
    if mod == 'ms':
        r_temp = -limit
        while r_temp < limit + delta_b:
            b.append(r_temp)
            r_temp += delta_b
        traj_tot = geometries_output(icharge, step,0, atoms_tot, b, d, coord_rotate_t, coord_rotate_b, [0, 0, 0], program, r_oruhf, dt, mult, velocities_m_s)
    elif mod == 'coll':
        r_temp = 0
        while r_temp <= limit:
            if round(r_temp, 3) != round(0, 3):  # For the head-on collision, the cross section is null. So it is negectable
                b.append(r_temp)
            r_temp += delta_b
        theta_tot_angles = angles(numero_angoli)
        theta_tot_angles_degrees = empty_matrix(3, numero_angoli)
        for i in range(0, numero_angoli):
            for j in range(0, 3):
                theta_tot_angles_degrees[i][j] = math.degrees(theta_tot_angles[i][j])
        count_angle = 0
        counting_traj = 0 
        while count_angle < numero_angoli:
            coord_fin_targ = rotation(coord_rotate_t, theta_tot_angles[count_angle])
            counting_traj = geometries_output(icharge, step, counting_traj, atoms_tot, b, d, coord_fin_targ, coord_rotate_b, theta_tot_angles_degrees[count_angle], program, r_oruhf, dt, mult, velocities_m_s)
            count_angle += 1
        traj_tot = len(b)*numero_angoli
    # Input dftb+
    if program == 'DFTB+':
        input_dftb(dt, path, atoms_tot, step, temp_conv, spin_spaiati, charge)
        with open("VELOC.DAT", "w") as f:
            for i in range(len(atoms_tot)):
                f.write("{} {} {}\n".format(velocities_m_s[0][i], velocities_m_s[1][i], velocities_m_s[2][i]))

    system.main(mod, numero_angoli, program,scr_dir, coord_target[0].size, coord_bull[0].size,traj_tot, path_program)


main()
