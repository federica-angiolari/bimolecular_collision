#    This script analyzes all trajectories in your working dir.  
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


# Constant
conversion_bohrps_ms = 52.9177
conversion_kcalmol_Jmol = 4.184
conversion_eV_Jmol = 96.48533
conversion_Jmol_eVmol = 6.24150934e18
conversion_Jmol_kcalmol = 0.2390057
avogadro = 6.02214e23
conversion_angstromps_ms = 1e2
cut_off = 2.0  # Angstrom
mass = {"H": 1.007825032e-3, "C": 12.000e-3, "N": 14.003e-3, "O": 15.9949146e-3, "F": 18.998e-3}  # Kg/mol


def cm_finder(matrix0, label):

    m_array, tmpx, tmpy, tmpz = [], [], [], []
    for i in range(0, len(label)):
        tmpx.append(mass[label[i]] * matrix0[0][i])
        tmpy.append(mass[label[i]] * matrix0[1][i])
        tmpz.append(mass[label[i]] * matrix0[2][i])
        m_array.append(mass[label[i]])
    cmx = sum(tmpx)/sum(m_array)
    cmy = sum(tmpy)/sum(m_array)
    cmz = sum(tmpz)/sum(m_array)
    return cmx, cmy, cmz


def f7(seq):

    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


### Start: Find products 
def first_cycle(atom_choose, coord, index_check):

    index_bond, frammento = [], []
    index_bond.append(atom_choose)
    number_atoms = len(coord[0])
    for i in range(0, number_atoms):
        if i != atom_choose:
            temporary_val = get_distance(coord[0][0], coord[0][i], coord[1][0], coord[1][i], coord[2][0], coord[2][i])
            if temporary_val <= 2.0 and temporary_val != 0.0:
                index_bond.append(i)
    continue_index = second_cycle(coord, index_bond, index_check)
    for l in range(0, len(continue_index)):
        index_bond.append(continue_index[l])
    index_bond = f7(index_bond)
    for i in range(0, len(index_bond)):
        coord[0][index_bond[i]] = 1e20
        coord[1][index_bond[i]] = 1e20
        coord[2][index_bond[i]] = 1e20
    return coord, index_bond


def second_cycle(coord, index, index_check):

    temporary_distance = []
    for j in range(len(index)):
        for i in range(len(coord[0])):
            if i not in index_check and i not in index and coord[0][i]!=1e20:

                temporary = get_distance(coord[0][index[j]], coord[0][i], coord[1][index[j]], coord[1][i], coord[2][index[j]], coord[2][i])
                if temporary != 0.00:
                    temporary_distance.append(temporary)
                if temporary <= cut_off and temporary != 0.00:
                    index.append(i)
    if all(i > cut_off for i in temporary_distance):
        return index
    else:
        return second_cycle(coord, index, index_check)


def number_diff(coordinate_to_check):
    number = 0
    index_to_take = []
    for i in range(0, len(coordinate_to_check)):
        if coordinate_to_check[i] != 1e20:
            index_to_take.append(i)
            number += 1
    return number, index_to_take


def coordinate_check(coord_origin):
    coord_matrix = [[], [], []]
    for i in range(0, 3):
        for k in range(0, len(coord_origin[0])):
            coord_matrix[i].append(coord_origin[i][k])
    i = 0
    index = []
    index_to_check = []
    coord_matrix, index_to_app = first_cycle(i, coord_matrix, index_to_check)
    for i in range(0, len(index_to_app)):
        index_to_check.append(index_to_app[i])
    index.append(index_to_app)
    i += 1
    n, ind = number_diff(coord_matrix[0])
    if n == 1:
        index.append(ind)
        for i in range(0, len(coord_matrix[0])):
            coord_matrix[0][i] = 1e20 
        return "two", index
    if any(i != 1e20 for i in coord_matrix[0]):
        coord_matrix, index_to_app = first_cycle(i, coord_matrix, index_to_check)
        index.append(index_to_app)
        for i in range(0, len(index_to_app)):
            index_to_check.append(index_to_app[i])
        return "two", index
    else:
        return "one", index
    i += 1
    if any(i != 1e20 for i in coord_matrix[0]):
        coord_matrix, index_to_app = first_cycle(i, coord_matrix, index_to_check)
        for i in range(0, len(index_to_app)):
            index_to_check.append(index_to_app[i])
        index.append(index_to_app)
        
        return "four", index
    else:
        return "three", index
# END

def empty_matrix(rows, col):

    matrix = [float(0)] * col
    for i in range(col):
        matrix[i] = [float(0)] * rows
    return np.array(matrix)


def in_and_fin(matrix, n):
    matrix_new_in, matrix_new_fin = [[], [], []], [[], [], []]
    for i in range(0, 3):
        for k in range(0, n):
            matrix_new_in[i].append(matrix[i][0][k])
            matrix_new_fin[i].append(matrix[i][len(matrix[0])-1][k])
    return np.array(matrix_new_in), np.array(matrix_new_fin)


def split_array(matrix2, number_to_split):

    matrix2_new = [[], [], []]
    for i in range(0, len(matrix2)):
        matrix2_new[i] = np.split(matrix2[i], len(matrix2[i])/number_to_split)
    return np.array(matrix2_new)


def extract_bullet_target(matrix3, num_1, num_2):

    num_tot = num_1 + num_2
    matrix3 = split_array(matrix3, num_tot)
    matrix_1, matrix_2 = [[], [], []], [[], [], []]

    for l in range(0, len(matrix3[0])):
        n = 0
        for k in range(0, num_1):
            matrix_1[0].append(matrix3[0][l][n])
            matrix_1[1].append(matrix3[1][l][n])
            matrix_1[2].append(matrix3[2][l][n])
            n += 1
        for k in range(0, num_2):
            matrix_2[0].append(matrix3[0][l][n])
            matrix_2[1].append(matrix3[1][l][n])
            matrix_2[2].append(matrix3[2][l][n])
            n += 1
    return np.array(matrix_1), np.array(matrix_2)


def get_distance(x1, x2, y1, y2, z1, z2):

    return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)


def get_matrix_distance(matrix4):

    distance_matrix = []
    for i in range(0, len(matrix4[0])-1):
        for j in range(i+1, len(matrix4[0])):
            distance_matrix.append(get_distance(matrix4[0][i], matrix4[0][j], matrix4[1][i], matrix4[1][j], matrix4[2][i], matrix4[2][j]))
    return np.array(distance_matrix)


def reaction_type(coord_2_in, coord_2_fin, coord_1_in, coord_1_fin):

    flag_t, flag_p, flag_tp = 0, 0, 0
    dist_target_in = get_matrix_distance(coord_2_in)
    dist_target_fin = get_matrix_distance(coord_2_fin)
    dist_bullett_in = get_matrix_distance(coord_1_in)
    dist_bullett_fin = get_matrix_distance(coord_1_fin)
    dist_bullet_target = []
    for i in range(0, len(coord_2_fin[0])):
        for j in range(0, len(coord_1_fin[0])):
            dist_bullet_target.append(get_distance(coord_2_fin[0][i], coord_1_fin[0][j], coord_2_fin[1][i], coord_1_fin[1][j], coord_2_fin[2][i], coord_1_fin[2][j]))
    for i in range(0, len(dist_target_fin)):
        for j in range(len(dist_target_in)):
            if abs(dist_target_fin[i]/dist_target_in[i]) > cut_off:
                flag_t = 1
    for i in range(0, len(dist_bullett_fin)):
        for j in range(0, len(dist_bullett_in)):
            if abs(dist_bullett_fin[i]/dist_bullett_in[i]) > cut_off:
                flag_p = 1
    for i in range(0, len(dist_bullet_target)):
        if dist_bullet_target[i] <= cut_off:
            flag_tp = 1

    with open("out", "a") as f:
        if flag_p == 1 or flag_t == 1 or flag_tp == 1:
            f.write("Y\n")
            return "Y"
        elif flag_p == 0 and flag_t == 0 and flag_tp == 0:
            f.write("N\n")
            return "N"


def en_pot_dftbplus(n):

    lines = []
    lines2 = []
    with open("md_{}.out".format(n), "r") as f:
        for line in f:
            if "Kinetic" in line:
                lines.append(line)
            elif "Pot" in line:
                lines2.append(line)
    with open("temp_cin", "w") as f:
        for i in range(len(lines)):
            f.write("{}".format(lines[i]))
    with open("temp_pot", "w") as f:
        for i in range(len(lines2)):
            f.write("{}".format(lines2[i]))

    kin = np.loadtxt("temp_cin", dtype=np.float, usecols=[5])  # eV
    pot = np.loadtxt("temp_pot", dtype=np.float, usecols=[4])  # eV
    with open("en_pot{}".format(n), "w") as f:
        for i in range(len(pot)):
            f.write("{}\t{}\n".format(i, pot[i]))
    with open("kin_{}".format(n), "w") as f:
        for i in range(len(kin)):
            f.write("{}\t{}\n".format(i, kin[i]))

    return kin, pot


def en_pot_gamess(nf, n):

    new_line2 = []
    new_line = []
    with open(nf, "r") as f:
        for line in f:
            if "POT." in line:
                new_line.append(line)
            elif "KIN." in line:
                new_line2.append(line)
    with open("temp_enpot", "w") as f:
        for i in range(len(new_line)):
            f.write("{}".format(new_line[i]))
    with open("temp_cin", "w") as f:
        for i in range(len(new_line2)):
            f.write("{}".format(new_line2[i]))
    kin = np.loadtxt("temp_cin", usecols=[2])
    pot = np.loadtxt("temp_enpot", usecols=[2])
    with open("en_pot{}".format(n), "w") as f:
        for i in range(len(pot)):
            f.write("{}\t{}\n".format(i, pot[i]))
    with open("kin_{}".format(n), "w") as f:
        for i in range(len(kin)):
            f.write("{}\t{}\n".format(i, kin[i]))
    return kin, pot  # kcal/mol


def conversion_matrix(matrix_to_convert, conversion_factor):
    matrix_converted = [[], [], []]
    for l in range(0, 3):
        for i in range(0, len(matrix_to_convert[0])):
            matrix_converted[l].append(matrix_to_convert[l][i] * conversion_factor)
    return np.array(matrix_converted)


def analysis_energy(program, number_products, flag_reaction, b_impact, number_traj, label_b, label_t, label_1, label_2, vb_i, vt_i, vf_1, vf_2, kin_in, kin_fin):
    with open("file_analisi", "a") as f:
        f.write("Traj: {}\n".format(number_traj))
        f.write("Reaction, {}: ".format(flag_reaction))
        for i in range(0, len(label_b)):
            f.write("{}".format(label_b[i]))
        f.write(" + ")
        for i in range(0, len(label_t)):
            f.write("{}".format(label_t[i]))
        f.write(" -> ")

    vcmx_t, vcmy_t, vcmz_t = cm_finder(vt_i, label_t)
    vcmx_b, vcmy_b, vcmz_b = cm_finder(vb_i, label_b)
    
    mass_array_t, mass_array_b = [], []
    for i in range(0, len(label_b)):
        mass_array_b.append(mass[label_b[i]])
    for i in range(0, len(label_t)):
        mass_array_t.append(mass[label_t[i]])

    if number_products == "two":
        vcmx_1, vcmy_1, vcmz_1 = cm_finder(vf_1, label_1)
        vcmx_2, vcmy_2, vcmz_2 = cm_finder(vf_2, label_2)

        mass_array_1, mass_array_2 = [], []
        for i in range(0, len(label_1)):
            mass_array_1.append(mass[label_1[i]])
        for i in range(0, len(label_2)):
            mass_array_2.append(mass[label_2[i]])

        trasl_kin_in = 0.5 * sum(mass_array_b) * (vcmx_b**2 + vcmy_b**2 + vcmz_b**2) + 0.5 * sum(mass_array_t) * (vcmx_t**2 + vcmy_t**2 + vcmz_t**2)
        trasl_kin_fin = 0.5 * sum(mass_array_1) * (vcmx_1**2 + vcmy_1**2 + vcmz_1**2) + 0.5 * sum(mass_array_2) * (vcmx_2**2 + vcmy_2**2 + vcmz_2**2)
        with open("file_analisi", "a") as f:
            for i in range(0, len(label_1)):
                f.write("{}".format(label_1[i]))
            f.write(" + ")
            for i in range(0, len(label_2)):
                f.write("{}".format(label_2[i]))
            f.write("\n")

    elif number_products == 'one':
        vcmx1, vcmy1, vcmz1 = cm_finder(vf_1, label_1)
        mass_array_1 = []
        for i in range(0, len(label_1)):
            mass_array_1.append(mass[label_1[i]])
        
        trasl_kin_in = 0.5 * sum(mass_array_b) * (vcmx_b**2 + vcmy_b**2 + vcmz_b**2) + 0.5 * sum(mass_array_t) * (vcmx_t**2 + vcmy_t**2 + vcmz_t**2)
        trasl_kin_fin = 0.5 * sum(mass_array_1) * (vcmx1**2 + vcmy1**2 + vcmz1**2)
        with open("file_analisi", "a") as f:
            for i in range(0, len(label_1)):
                f.write("{}".format(label_1[i]))
            f.write("\n")

    if program == 'DFTB+':
        trasl_kin_in = trasl_kin_in * conversion_Jmol_eVmol/avogadro
        trasl_kin_fin = trasl_kin_fin * conversion_Jmol_eVmol/avogadro
        unit = 'eV'
    elif program == 'GAMESS':
        trasl_kin_in = trasl_kin_in * conversion_Jmol_kcalmol
        trasl_kin_fin = trasl_kin_fin * conversion_Jmol_kcalmol
        unit = 'kcal/mol'

    rotovib_in = kin_in - trasl_kin_in
    rotovib_fin = kin_fin - trasl_kin_fin
    with open("file_analisi", "a") as f:
        f.write("Numero traiettoria: {0}\n".format(number_traj))
        f.write("Fattore di impatto [Angstrom] = {}\n\nEnergia cinetica traslazionale + rotovibrazionale, {}\nIniziale: {:.3e}\nFinale: {:.3e}\n\n "
                "Energia cinetica traslazionale, {}:\nIniziale: {:.3e}\nFinale: {:.3e}\n\n "
                "Energia rotovibrazionale, {}:\niniziale: {:.3e}\n"
                "finale:{:.3e}\n\n"
                .format(b_impact, unit, kin_in, kin_fin,  unit, trasl_kin_in, trasl_kin_fin, unit, rotovib_in, rotovib_fin, ))
        f.write("\n")
        for i in range(100):
            f.write("_")
        f.write("\n\n\n")

    if flag_reaction == 'N':
        with open("bi_etrasf_tot_N", "a") as f:
            f.write("{}\t{}\n".format(b_impact, abs((kin_fin-kin_in)/kin_in)*100))
        with open("bi_etrasf_trasl_N", "a") as f:
            f.write("{}\t{}\n".format(b_impact, abs((trasl_kin_fin-trasl_kin_in)/trasl_kin_in)*100))
        with open("bi_etrasf_rotovib_N", "a") as f:
            f.write("{}\t{}\n".format(b_impact, abs((rotovib_fin-rotovib_in)/rotovib_in)*100))
    else:
        with open("bi_etrasf_tot_Y", "a") as f:
            f.write("{}\t{}\n".format(b_impact, abs((kin_fin-kin_in)/kin_in)*100))
        with open("bi_etrasf_trasl_Y", "a") as f:
            f.write("{}\t{}\n".format(b_impact, abs((trasl_kin_fin-trasl_kin_in)/trasl_kin_in)*100))
        with open("bi_etrasf_rotovib_Y", "a") as f:
            f.write("{}\t{}\n".format(b_impact, abs((rotovib_fin-rotovib_in)/rotovib_in)*100))


def conversion_array(array_to_convert, conversion_factor):
    array_converted = []
    for i in range(array_to_convert.size):
        array_converted.append(array_to_convert[i] * conversion_factor)
    return np.array(array_converted)


def angoli(flag, v, label, b_traj):
    vcmx, vcmy, vcmz = cm_finder(v, label)
    vcm_scal = np.sqrt(vcmx**2 + vcmy**2 + vcmz**2)
    theta_2 = math.acos(vcmx/vcm_scal)
    theta_2 = math.degrees(theta_2)
    theta_2 = 180-theta_2
    if flag == 'N':
        with open("angle_vs_b", "a") as f:
            f.write("{}\t{}\n".format(b_traj, theta_2))
    if flag == 'Y':
        with open("map_3d", "a") as f:
            f.write("{}\t{}\t{}\n".format(vcmx, vcmy, vcmz))


def cross_section(num_ang, b):

    y_or_n = np.loadtxt("out", dtype=np.str, usecols=[0])

    prob = to_prob(y_or_n)
    bi = list(dict.fromkeys(b))  # b singoli e non ripetuti

    prob = np.array(prob)
    prob = np.split(prob, num_ang)  # prob divise per angolo

    prob_4_bi = [0]*len(bi)
    for i in range(len(bi)):
        prob_4_bi[i] = [0]*num_ang
    prob_4_bi = np.array(prob_4_bi)
    j = 0
    while j < num_ang:
        for i in range(len(bi)):
            prob_4_bi[i][j] = prob[j][i]
        j += 1

    probability = []
    for i in range(len(bi)):
        count = 0
        for j in range(num_ang):
            if prob_4_bi[i][j] == 1:
                count += 1
        probability.append(float(count/num_ang))
    dx = bi[1] - bi[0]
    fun = []
    for i in range(len(probability)):
        fun.append(bi[i]*probability[i])

    area2 = np.trapz(fun, bi, dx)

    with open("file_analisi", "a") as f:
        for i in range(100):
            f.write("_")
        f.write("\nCROSS SECTION")
        f.write("\n{:.4e}[Angstrom^2]".format(area2*2*math.pi))
    with open("graph_p_b", "w") as f:
        for i in range(len(bi)):
            f.write("{} {}\n".format(bi[i], probability[i]))


def to_prob(v):

    prob = []
    for i in range(len(v)):
        if v[i] == "Y":
            prob.append(1)
        else:
            prob.append(0)
    return prob


def main(num_tot_traj, num_atoms, num_p, program, n_angoli):

    num_t = num_atoms - num_p
    b = np.loadtxt("summ_traj", dtype=float, usecols=[3])
    for number_ist in range(0, num_tot_traj):
        print("Analysis traj: {}".format(number_ist))
        if program == 'GAMESS':
            name_traj = "geometry_{}.trj".format(number_ist)
            with open(name_traj, "r") as f:
                linee = f.read().splitlines()

            with open("temporary.xyz", "w") as f:  # Temporary file with all coordinates in Angstrom
                for i in range(0, len(linee)):
                    if linee[i] == '----- QM PARTICLE COORDINATES FOR $DATA GROUP -----':
                        i += 1
                        for j in range(0, num_atoms):
                            f.write("{}\n".format(linee[i]))
                            i += 1
            atoms_label = np.loadtxt("temporary.xyz", dtype=np.str, usecols=[0])
            label = []
            for i in range(num_atoms):
                label.append(atoms_label[i])
            coord_total = [[], [], []]
            col = 0
            for i in range(2, 5):
                coord_total[col] = np.loadtxt("temporary.xyz", dtype=float, usecols=[i])
                col += 1
            with open("velocity_gamess", "w") as f:  # New temporary file with all the velocities, in Bohr/ps
                for i in range(0, len(linee)):
                    if linee[i] == "TVELQM(1)=     ! QM ATOM TRANS. VELOCITIES (BOHR/PS) !":
                        i += 1
                        for j in range(num_atoms):
                            f.write("{}\n".format(linee[i]))
                            i += 1
            velocity_total = [[], [], []]
            col = 0
            for i in range(0, 3):
                velocity_total[col] = np.loadtxt("velocity_gamess", dtype=float, usecols=[i])
                col += 1
            with open("traj_{}.xyz".format(number_ist), "w") as f:
                r = 0
                c = 0
                while r < len(coord_total[0]):
                    j = 0
                    f.write("{}\nMD iter: {}\n".format(num_atoms, c))
                    c += 1
                    while j < num_atoms:
                        f.write("{} {} {} {} {} {} {}\n".format(atoms_label[r], coord_total[0][r], coord_total[1][r], coord_total[2][r], velocity_total[0][r], velocity_total[1][r], velocity_total[2][r]))
                        j += 1
                        r += 1
                    r = r
            kinetic_energy_arr, pot_energy_array = en_pot_gamess(name_traj, number_ist)
            velocity_total = conversion_matrix(velocity_total, conversion_bohrps_ms)

        elif program == 'DFTB+':
            name_traj = 'NVE_{}.xyz'.format(number_ist)
            kinetic_energy_arr, pot_energy_array = en_pot_dftbplus(number_ist)
            with open(name_traj, "r") as f:
                linee = f.read().splitlines()
            coord_total = [[], [], []]
            velocity_total = [[], [], []]
            col1 = 0
            col2 = 5
            r = 0
            label = np.loadtxt(name_traj, dtype=np.str, skiprows=len(linee)-num_atoms, usecols=[0])
            with open("temporary", "w") as f:  # crea un file con label, x, y, z, electrons, vx, vy, vz
                while r < len(linee):
                    if num_atoms >= 10:
                        if linee[r] == "   {}".format(num_atoms):
                            r += 2
                            for i in range(0, num_atoms):
                                f.write("{}\n".format(linee[r]))
                                r += 1
                    else:
                        if linee[r] == "    {}".format(num_atoms):
                            r += 2
                            for i in range(0, num_atoms):
                                f.write("{}\n".format(linee[r]))
                                r += 1

            for i in range(1, 4):
                coord_total[col1] = np.loadtxt("temporary", dtype=float, usecols=[i])  # Angstrom
                velocity_total[col1] = np.loadtxt("temporary", dtype=float, usecols=[col2])  # Angstrom/ps
                col1 += 1
                col2 += 1
            velocity_total = conversion_matrix(velocity_total, conversion_angstromps_ms)

        coord_bullett, coord_target = extract_bullet_target(coord_total, num_p, num_t)  
        velocity_bullett, velocity_target = extract_bullet_target(velocity_total, num_p, num_t)

        coord_bullett = split_array(coord_bullett, num_p)
        velocity_bullett = split_array(velocity_bullett, num_p)
        coord_target = split_array(coord_target, num_t)
        velocity_target = split_array(velocity_target, num_t)

        velocity_target_in, velocity_target_fin = in_and_fin(velocity_target, num_t)
        coord_target_in, coord_target_fin = in_and_fin(coord_target, num_t)
        coord_bullett_in, coord_bullett_fin = in_and_fin(coord_bullett, num_p)
        velocity_bullett_in, velocity_bullett_fin = in_and_fin(velocity_bullett, num_p)

        coord_total = split_array(coord_total, num_atoms)
        coord_total_in, coord_total_fin = in_and_fin(coord_total, num_atoms)

        label_bull, label_targ = [], []
        for i in range(0, num_p):
            label_bull.append(label[i])
        for i in range(num_p, num_atoms):
            label_targ.append(label[i])
    
        flag = reaction_type(coord_target_in, coord_target_fin, coord_bullett_in, coord_bullett_fin)

        number_prod, index_prod = coordinate_check(coord_total_fin)

        velocity_total = split_array(velocity_total, num_atoms)
        velocity_total_in, velocity_total_fin = in_and_fin(velocity_total, num_atoms)

        xcm1, ycm1, zcm1 = cm_finder(coord_total_in, label)
        xcm2, ycm2, zcm2 = cm_finder(coord_total_fin, label)
        dist_cm =  np.sqrt((xcm2-xcm1)**2 + (ycm2-ycm1)**2 + (zcm2-zcm1)**2)
        with open("distance_cm_f_i", "a") as f:
            f.write("{}\n".format(dist_cm))

        
        if number_prod == 'two':
            velocity_product_one = [[], [], []]
            velocity_product_two = [[], [], []]
            coord_product_one = [[], [], []]
            coord_product_two = [[], [], []]
            label_one = []
            label_two = []

            for k in range(0, 3):
                for j in range(0, len(index_prod[0])):
                    velocity_product_one[k].append(velocity_total_fin[k][index_prod[0][j]])
                    coord_product_one[k].append(coord_total_fin[k][index_prod[0][j]])
                
                for l in range(0, len(index_prod[1])):
                    velocity_product_two[k].append(velocity_total_fin[k][index_prod[1][l]])
                    coord_product_two[k].append(coord_total_fin[k][index_prod[1][l]])
            for i in range(0, len(index_prod[0])):
                label_one.append(label[index_prod[0][i]])
            for i in range(0, len(index_prod[1])):
                label_two.append(label[index_prod[1][i]])
        
            if program == 'GAMESS':
                analysis_energy(program, number_prod, flag, b[number_ist], number_ist, label_bull, label_targ, label_one, label_two,velocity_bullett_in, velocity_target_in, velocity_product_one, velocity_product_two, kinetic_energy_arr[1], kinetic_energy_arr[len(kinetic_energy_arr)-1])

            else:
                analysis_energy(program, number_prod, flag, b[number_ist], number_ist, label_bull, label_targ, label_one, label_two,velocity_bullett_in, velocity_target_in, velocity_product_one, velocity_product_two, kinetic_energy_arr[0], kinetic_energy_arr[len(kinetic_energy_arr)-1])

            if len(label_one) > len(label_two):
                angoli(flag, velocity_product_one, label_one, b[number_ist])
            elif len(label_two) > len(label_one):
                angoli(flag, velocity_product_two, label_two, b[number_ist])
            else:
                angoli(flag, velocity_product_two, label_two, b[number_ist])  
        elif number_prod == 'one':
            if program == 'GAMESS':
                analysis_energy(program, number_prod, flag, b[number_ist], number_ist, label_bull, label_targ, label, "XX", velocity_bullett_in, velocity_target_in, velocity_total_fin, "XX", kinetic_energy_arr[1], kinetic_energy_arr[len(kinetic_energy_arr)-1])
            else:
                analysis_energy(program, number_prod, flag, b[number_ist], number_ist, label_bull, label_targ, label, "XX", velocity_bullett_in, velocity_target_in, velocity_total_fin, "XX", kinetic_energy_arr[0], kinetic_energy_arr[len(kinetic_energy_arr)-1])

                
    cross_section(n_angoli, b)
