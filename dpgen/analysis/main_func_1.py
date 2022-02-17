#!/usr/bin/env python3

from dpdata import LabeledSystem

from dpgen.analysis.analysis_utils import select_logs, process_bar

###########################################################
def data_collections():
    while True:
        print("\n")
        print("========================================================")
        print("                    Main Menu 1                        ")
        print("                                                       ")
        print(" ( 11)  Collect all data in current dir from OUTCAR file")
        print("\n Tips: Input q or -10 to exit program                ")
        print("========================================================")
        
        jsel = input("\n Please input the menu index: ")
        if jsel == str(-10):
            break
        elif jsel == str(11):
            root_dir = input("\n Please input the root dir name: ")
            str_in   = input("\n Collect (1) all frames (e.g., AIMD task) or (2) final frame (e.g., Optimization task)? ")
            out_dir  = input("\n Please input the out dir name: ")
            total_frames = -1
            if str_in == str(1):
                total_frames = 0
            collect_all_data_in_current_dir_from_all_OUTCAR_file(root_dir=root_dir, total_frames=total_frames, out_dir=out_dir)
###########################################################


###########################################################
def collect_all_data_in_current_dir_from_all_OUTCAR_file(root_dir: str, total_frames: int, out_dir: str) -> LabeledSystem:
    all_outcar = select_logs(root_dir, "OUTCAR")
    print(" Totally find %d OUTCAR in current dir."%len(all_outcar))        

    final_sys = None
    counter = 0
    for outcar in all_outcar:
        counter += 1
        labeled_system = LabeledSystem(outcar, fmt="vasp/outcar")
        if total_frames == -1:
            frames = [-1, ]
        else:
            frames = range(labeled_system.get_nframes())
        for iframe in frames:
            if final_sys is None:
                final_sys = labeled_system[iframe]
            else:
                final_sys.append(labeled_system[iframe])
        process_bar(counter, len(all_outcar), " Read %s"%outcar)
    
    final_sys.to_deepmd_raw(out_dir)
    final_sys.to_deepmd_npy(out_dir)
###########################################################
