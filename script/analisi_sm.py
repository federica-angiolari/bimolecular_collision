#    This script analyze fragmentation path 
#    Read text file in man dir for furhter information
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


cut_off = 2.00  # Angstrom
atoms_el = {"N": 5.00, "C": 4.00, "O": 6.00, "H": 1.00, "F": 7.00, "Br": 7.00}
mass = {"H": 1.007825032, "C": 12.000, "N": 14.003074004, "O": 15.994914619, "F": 18.998403163, "Br": 78.9183376}


def f7(seq):

    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def num_lines_file(name_file):

    lines_number = 0
    with open(name_file, "r") as f:
        for line in f:
            lines_number += 1
    return lines_number


def empty_matrix(rows, col):

    matrix = [float(0)] * col
    for i in range(col):
        matrix[i] = [float(0)] * rows
    return np.array(matrix)


def get_distance(x1, x2, y1, y2, z1, z2):

    return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)


def count_charge(atomi, elettroni):

    charge = []
    for i in range(len(atomi)):
        if atomi[i] != 1e20:
            charge.append(atoms_el[atomi[i]] - elettroni[i])
    return charge


def first_cycle(coord, electrons, label):

    electron_bond, index_bond, frammento = [], [], []
    index_bond.append(0)

    number_atoms = len(coord[0])
    for i in range(1, number_atoms):
        temporary_val = get_distance(coord[0][0], coord[0][i], coord[1][0], coord[1][i], coord[2][0], coord[2][i])
        if temporary_val <= 2.0: 
            index_bond.append(i)
    continue_index = second_cycle(coord, index_bond)
    for l in range(0, len(continue_index)):
        index_bond.append(continue_index[l])
    index_bond = f7(index_bond) 
    for i in range(0, len(index_bond)):
        frammento.append(label[index_bond[i]])
        electron_bond.append(electrons[index_bond[i]])
    for i in range(0, len(index_bond)):  
        coord[0][index_bond[i]] = 1e20
        coord[1][index_bond[i]] = 1e20
        coord[2][index_bond[i]] = 1e20
        label[index_bond[i]] = 1e20
        electrons[index_bond[i]] = 1e20
    final_charge = count_charge(frammento, electron_bond)
    with open("Frammenti", "a") as f:
        f.write("\n")
        for i in range(0, len(frammento)):
            f.write("{}".format(frammento[i]))
        f.write("\t{}".format(round(sum(final_charge), 3)))
    return coord, label, electrons


def second_cycle(coord, index):

    temporary_distance = []
    for j in range(len(index)):  
        for i in range(len(coord[0])):  
            if i not in index:
                temporary = get_distance(coord[0][index[j]], coord[0][i], coord[1][index[j]], coord[1][i], coord[2][index[j]], coord[2][i])
                if temporary != 0.00:  
                    temporary_distance.append(temporary)
                if temporary <= cut_off and temporary != 0.00:
                    index.append(i)  
    if all(i > cut_off for i in temporary_distance):  
        return index
    else:
        return second_cycle(coord, index)  


def check_coord_label_electrons(coord, label, electrons):

    coord1, label1, electrons1 = [[], [], []], [], []
    for i in range(len(coord[0])):
        if coord[0][i] != 1e20:
            coord1[0].append(coord[0][i])
            coord1[1].append(coord[1][i])
            coord1[2].append(coord[2][i])
            label1.append(label[i])
            electrons1.append(electrons[i])
    return coord1, label1, electrons1


def massa_function(stringa):

    m = []
    for i in range(len(stringa)):
        m.append(mass[stringa[i]])
    return sum(m)


def frammenti_to_list(framm):

    n = list(framm)
    return massa_function(n)


def spettro(mass_tot):

    with open("Frammenti", "r") as f:
        lines = f.read().splitlines()
    l_fin = []  
    for i in range(len(lines)):
        if lines[i] == '':
            i += 1
        elif 'Traj' in lines[i]:
            i += 1
        else:
            l_fin.append(lines[i])
    with open("temporary_file", "w") as f:
        for i in range(len(l_fin)):
            f.write("{}\n".format(l_fin[i]))
    frammenti = np.loadtxt("temporary_file", dtype=np.str, usecols=[0])  
    cariche_parziali = np.loadtxt("temporary_file", dtype=float, usecols=[1]) 
    frammenti_pos, frammenti_neg, frammenti_neutr = [], [], []
    for i in range(0, len(cariche_parziali)):
        if round(cariche_parziali[i],2) > 0.65:
            frammenti_pos.append(frammenti[i])
        elif round(cariche_parziali[i], 2) < 0.00:
            frammenti_neg.append(frammenti[i])
        else:
            frammenti_neutr.append(frammenti[i])
    n_positive = len(frammenti_pos)
    n_negative = len(frammenti_neg)
    n_neutral = len(frammenti_neutr)

    col1 = []
    temporary_var = 0.000
    while temporary_var <= mass_tot:
        col1.append(round(temporary_var, 3))
        temporary_var += 0.001
    massa_frammenti_tot = []
    for i in range(0, len(frammenti)):
        massa_frammenti_tot.append(round(frammenti_to_list(frammenti[i]),3))
    col2_tot = []
    for i in range(0, len(col1)):
        col2_tot.append(massa_frammenti_tot.count(col1[i]))

    max_col2_tot = max(col2_tot)
    for i in range(0, len(col2_tot)):
        col2_tot[i] = col2_tot[i]/max_col2_tot
    with open("spectrum", "w") as f:
        for i in range(0, len(col2_tot)):
            f.write("{} {}\n".format(col1[i], col2_tot[i]*100))

    massa_frammenti_pos, massa_frammenti_neg, massa_frammenti_neutr = [], [], []
    if n_positive != 0:
        for i in range(0, len(frammenti_pos)):
            massa_frammenti_pos.append(round(frammenti_to_list(frammenti_pos[i]),3))  
        col2_pos = []
        for i in range(0, len(col1)):
            col2_pos.append(massa_frammenti_pos.count(col1[i]))
        max_col2_pos = max(col2_pos)
        for i in range(0, len(col2_pos)):
            col2_pos[i] = col2_pos[i]/max_col2_pos
        with open("spectrum_cation", "w") as f:
            for i in range(0, len(col1)):
                f.write("{} {}\n".format(col1[i], col2_pos[i]*100))
    if n_negative != 0:
        for i in range(0, len(frammenti_neg)):
            massa_frammenti_neg.append(round(frammenti_to_list(frammenti_neg[i]),3))
        col2_neg = []
        for i in range(0, len(col1)):
            col2_neg.append(massa_frammenti_neg.count(col1[i]))
        max_col2_neg = max(col2_neg)
        for i in range(0, len(col2_neg)):
            col2_neg[i] = col2_neg[i]/max_col2_neg
        with open("spectrum_anion", "w") as f:
            for i in range(0, len(col1)):
                f.write("{} {}\n".format(col1[i], col2_pos[i]*100))
    if n_neutral != 0:
        for i in range(0, len(frammenti_neutr)):
            massa_frammenti_neutr.append(round(frammenti_to_list(frammenti_neutr[i]),3))
        col2_neutr = []
        for i in range(0, len(col1)):
            col2_neutr.append(massa_frammenti_neutr.count(col1[i]))
        max_col2_neutr = max(col2_neutr)
        for i in range(0, len(col2_neutr)):
            col2_neutr[i] = col2_neutr[i]/max_col2_neutr
        with open("spectrum_neutral", "w") as f:
            for i in range(0, len(col1)):
                f.write("{} {}\n".format(col1[i], col2_neutr[i]*100))

    frammenti_list, indici_to_del = [], []
    for i in range(len(frammenti)):
        frammenti_list.append(list(frammenti[i]))
    for i in range(len(frammenti_list)-1):
        for j in range(i, len(frammenti_list)):
            if i != j:
                if frammenti_list[i] == frammenti_list[j]:
                    indici_to_del.append(i)
    indici_to_del = f7(indici_to_del)
    occurrance = []
    frammenti_old = []
    for i in range(len(frammenti)):
        frammenti_old.append(frammenti[i])
    for i in range(len(indici_to_del)):
        frammenti[indici_to_del[i]] = 'xxx'
    cariche_new, frammenti_new = [], []
    for i in range(len(frammenti)):
        if frammenti[i] != 'xxx':
            cariche_new.append(cariche_parziali[i])
            frammenti_new.append(frammenti[i])
    for i in range(0, len(frammenti_new)):
        occurrance.append(frammenti_old.count(frammenti_new[i]))
    massa_new_frammenti = []
    for i in range(0, len(frammenti_new)):
        massa_new_frammenti.append(frammenti_to_list(frammenti_new[i]))


    massa_new_frammenti, frammenti_new, cariche_new, occurrance = bubblesort(massa_new_frammenti, frammenti_new, cariche_new, occurrance)
    with open("table", "w") as f:
        f.write("      Fragments  \t  Average charge\t        Mass      \t      Number       \n\n")
        for i in range(0, len(frammenti_new)):
            if len(list(frammenti_new[i])) > 5:
                f.write("\t{}\t\t{}\t\t\t{}\t\t{}\n".format(frammenti_new[i], cariche_new[i], massa_new_frammenti[i], occurrance[i]))
            else:
                f.write("\t{}\t\t\t{}\t\t\t{}\t\t{}\n".format(frammenti_new[i], cariche_new[i], massa_new_frammenti[i], occurrance[i]))


def bubblesort(v1, v2, v3, v4):
    for i in range(0, len(v1)-1):
        for j in range(i+1, len(v1)):
            if v1[i] <= v1[j]:
                v1[i], v1[j] = v1[j], v1[i]
                v2[i], v2[j] = v2[j], v2[i]
                v3[i], v3[j] = v3[j], v3[i]
                v4[i], v4[j] = v4[j], v4[i]
    return v1, v2, v3, v4

def main(num_target, num_traj):

    for k in range(0, num_traj):
        print("Analysis Traj: {}".format(k))
        name_file = str("NVE_{}.xyz".format(k))
        with open("Frammenti", "a") as f:
            f.write("\nTraj: {}\n".format(k))
        total_lines_traj = num_lines_file(name_file)
        coord_matrix = empty_matrix(num_target, 3) 
        index_col = 0
        for i in range(1, 4):
            coord_matrix[index_col] = np.loadtxt(name_file, dtype=float, skiprows=total_lines_traj-num_target, usecols=[i])
            index_col += 1
        label = np.loadtxt(name_file, dtype=np.str, skiprows=total_lines_traj - num_target, usecols=[0])
        mass_tot = 0
        for i in range(len(label)):
            mass_tot += mass[label[i]]
        electrons = np.loadtxt(name_file, dtype=float, skiprows=total_lines_traj - num_target, usecols=[4])
        while True:
            coord_matrix, label, electrons = first_cycle(coord_matrix, electrons, label)
            coord_matrix, label, electrons = check_coord_label_electrons(coord_matrix, label, electrons)
            if all(i == 1e20 for i in label):
                break
    spettro(mass_tot)
