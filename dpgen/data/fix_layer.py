import numpy as np

def layer(jdata, poscar_in, poscar_out):
    with open(poscar_in, 'r') as file_obj:
        lines = [line.strip() for line in file_obj.readlines()]

    total = sum(list(map(int, lines[6].split())))
    cord_str = lines[8: 8 + total]
    cord_cartes = [list(map(float, cord.split())) for cord in lines[8: 8 + total]]

    cord_sort = sorted(cord_cartes,key=lambda x:x[2])
    dis_small = cord_sort[0][2]
    dis_big = cord_sort[-1][2]
    max_dis = dis_big - dis_small
    lay_space = max_dis / jdata["total_layers"]

    if len(jdata["fix_layers"]) == 1:
        with open(poscar_out, 'w') as file_out:
            head = lines[:8]
            head.insert(7, 'Selective dynamics')
            for each in head:
                file_out.write(f"{each}\n")

            for idx, cord in enumerate(cord_cartes):
                threshold = dis_small + (jdata["fix_layers"][0]-1) * lay_space
                if cord[2] <= threshold:
                    cord_str[idx] = f"{cord_str[idx]}        F    F    F\n"
                else:
                    cord_str[idx] = f"{cord_str[idx]}        T    T    T\n"

            for each in cord_str:
                file_out.write(each)

    elif len(jdata["fix_layers"]) == 2:
        with open(poscar_out, 'w') as file_out:
            head = lines[:8]
            head.insert(7, 'Selective dynamics')
            for each in head:
                file_out.write(f"{each}\n")

            for idx, cord in enumerate(cord_cartes): 
                threshold_min = dis_small + (jdata["fix_layers"][0]-1) * lay_space
                threshold_max = dis_small + (jdata["fix_layers"][1]-1) * lay_space
                if threshold_min <= cord[2] <= threshold_max:
                    cord_str[idx] = f"{cord_str[idx]}        F    F    F\n"
                else:
                    cord_str[idx] = f"{cord_str[idx]}        T    T    T\n"

            for each in cord_str:
                file_out.write(each)