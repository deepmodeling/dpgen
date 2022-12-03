import dpdata
import numpy as np

default_config={
  "GLOBAL": {
    "PROJECT": "DPGEN"
  },
  "FORCE_EVAL": {
    "METHOD": "QS",
    "STRESS_TENSOR": "ANALYTICAL",
    "DFT": {
      "BASIS_SET_FILE_NAME": "./cp2k_basis_pp_file/BASIS_MOLOPT",
      "POTENTIAL_FILE_NAME": "./cp2k_basis_pp_file/GTH_POTENTIALS",
      "CHARGE": 0,
      "UKS": "F",
      "MULTIPLICITY": 1,
      "MGRID": {
        "CUTOFF": 400,
        "REL_CUTOFF": 50,
        "NGRIDS": 4
      },
      "QS": {
        "EPS_DEFAULT": "1.0E-12"
      },
      "SCF": {
        "SCF_GUESS": "ATOMIC",
        "EPS_SCF": "1.0E-6",
        "MAX_SCF": 50
      },
      "XC": {
        "XC_FUNCTIONAL": {
          "_": "PBE"
        }

      }
    },
    "SUBSYS": {
        "CELL":{
            "A": "10 .0 .0",
            "B": ".0 10 .0",
            "C": ".0 .0 10"
            },
      "COORD": {"@include": "coord.xyz"},
      "KIND": {
          "_": ["H","C","N"],
          "POTENTIAL": ["GTH-PBE-q1","GTH-PBE-q4", "GTH-PBE-q5"],
          "BASIS_SET": ["DZVP-MOLOPT-GTH","DZVP-MOLOPT-GTH","DZVP-MOLOPT-GTH"]
          }
      },
    "PRINT": {
      "FORCES": {
        "_": "ON"
    },
      "STRESS_TENSOR":{
          "_": "ON"
          }
  }
}
}

def update_dict(old_d, update_d):
    """
    a method to recursive update dict
    :old_d: old dictionary
    :update_d: some update value written in dictionary form
    """
    import collections.abc
    for k, v in update_d.items():
        if (k in old_d and isinstance(old_d[k], dict) and isinstance(update_d[k], collections.abc.Mapping)):
            update_dict(old_d[k], update_d[k])
        else:
            old_d[k] = update_d[k]

def iterdict(d, out_list, flag=None):
    """
    :doc: a recursive expansion of dictionary into cp2k input
    :k: current key
    :v: current value
    :d: current dictionary under expansion
    :flag: used to record dictionary state. if flag is None,
    it means we are in top level dict. flag is a string.
    """
    for k,v in d.items():
        k=str(k) # cast key into string
        #if value is dictionary
        if isinstance(v, dict):
            # flag == None, it is now in top level section of cp2k
            if flag==None :
                out_list.append("&"+k)
                out_list.append("&END "+k)
                iterdict(v, out_list, k)
            # flag is not None, now it has name of section
            else:
                index = out_list.index("&END " + flag)
                out_list.insert(index, "&"+k)
                out_list.insert(index+1,"&END "+k )
                iterdict(v, out_list, k)
        elif isinstance(v, list):
#            print("we have encountered the repeat section!")
            index = out_list.index("&"+flag)
            # delete the current constructed repeat section
            del out_list[index:index+2]
            # do a loop over key and corresponding list
            k_tmp_list = []
            v_list_tmp_list = []
            for k_tmp, v_tmp in d.items():
                k_tmp_list.append(str(k_tmp))
                v_list_tmp_list.append(v_tmp)
            for repeat_keyword in zip(*v_list_tmp_list):
                out_list.insert(index,"&" + flag)
                out_list.insert(index+1, "&END " + flag)
                for idx, k_tmp in enumerate(k_tmp_list):
                    if k_tmp == "_":
                        out_list[index] = "&" + flag + " " + repeat_keyword[idx]
                    else:
                        out_list.insert(index+1, k_tmp+" "+repeat_keyword[idx])

            break

        else:
            v=str(v)
            if flag==None :
                out_list.append(k+" "+v)
                print (k,":",v)
            else:
                if k == "_":
                    index = out_list.index("&" + flag)
                    out_list[index] = ("&" + flag + " " + v)

                else:
                    index = out_list.index("&END "+flag)
                    out_list.insert(index, k+" "+v)


def make_cp2k_input(sys_data, fp_params):
    #covert cell to cell string
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3,3])
    cell_a = np.array2string(cell[0,:])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1,:])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2,:])
    cell_c = cell_c[1:-1]

    #get update from user
    user_config=fp_params
    #get update from cell
    cell_config={"FORCE_EVAL":{
        "SUBSYS":{
            "CELL":{
                "A": cell_a,
                "B": cell_b,
                "C": cell_c
                }
            }
        }
            }
    update_dict(default_config, user_config)
    update_dict(default_config, cell_config)
    #output list
    input_str = []
    iterdict(default_config, input_str)
    string="\n".join(input_str)
    return string




def make_cp2k_xyz(sys_data):
    #get structral information
    atom_names = sys_data['atom_names']
    atom_types = sys_data['atom_types']

    #write coordinate to xyz file used by cp2k input
    coord_list = sys_data['coords'][0]
    u = np.array(atom_names)
    atom_list = u[atom_types]
    x = '\n'
    for kind, coord in zip(atom_list, coord_list) :
        x += str(kind) + ' ' + str(coord[:])[1:-1] + '\n'
    return x



def make_cp2k_input_from_external(sys_data, exinput_path):
    # read the input content as string
    with open(exinput_path, 'r') as f:
        exinput = f.readlines()

    # find the ABC cell string
    for line_idx, line in enumerate(exinput):
        if 'ABC' in line:
            delete_cell_idx = line_idx
            delete_cell_line = line

    # remove the useless CELL line
    exinput.remove(delete_cell_line)

    # insert the cell information
    # covert cell to cell string
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3,3])
    cell_a = np.array2string(cell[0,:])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1,:])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2,:])
    cell_c = cell_c[1:-1]

    exinput.insert(delete_cell_idx, 'A  ' + cell_a + '\n')
    exinput.insert(delete_cell_idx+1, 'B  ' + cell_b + '\n')
    exinput.insert(delete_cell_idx+2, 'C  ' + cell_c + '\n')

    return ''.join(exinput)


