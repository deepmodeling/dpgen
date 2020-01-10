import subprocess
import os
def stat_time(target_folder, 
            param_file = 'param.json', 
            verbose = True, 
            mute = False):
    train_dirs = subprocess.run([f"ls -d -1 {target_folder}/iter.??????/00.train/", ],
        shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n')
    for dir in train_dirs:
        abs_dir = os.path.abspath(dir)
        stage = os.path.basename(os.path.dirname(dir))
        train_time_logs = subprocess.run([f"grep -H --text 'wall time' {dir}/???/train.log", ],
            shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n')
        upload_task_dir_num = subprocess.run([f"ls -1 -d {dir}/??? |wc -l", ], 
            shell=True, stdout=subprocess.PIPE).stdout.decode().strip('\n')
        total_core_sec = float(0)
        
        # assume training on single GPU
        paral_cores = 1
        finished_task_file_num = len(train_time_logs)
        # gpu_type_set = set([])
        for log in train_time_logs:
            # log example : 
            #   .//iter.000000/00.train//003/train.log:# DEEPMD: wall time: 7960.265 s
            # print(log.split(':'))
            file_path, text1, text2, wall_time = log.split(':') # pylint: disable=unused-variable
            abs_file_path = os.path.abspath(file_path)
            # stage=='00.train'
            
            wall_time_sec = float(wall_time.strip('s').strip(' '))
            total_core_sec += wall_time_sec * paral_cores

            # r'd\nja\1lgje' leading 'r' means dont treat '\' as  Escape character
            # gpu_type = subprocess.run([fr"grep -e 'physical GPU' {abs_file_path} |sed -n -E -e 's|^.*name: (.*), pci.*|\1|p'", ],
            #     shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n').pop()
            # gpu_type_set.add(gpu_type)
        
        total_core_hour = total_core_sec * paral_cores / 3600 
        print(f"{stage}:{abs_dir}"
            f"paral_cores:{paral_cores}"
            f":upload_task_dir_num:{upload_task_dir_num}"
            f":finished_task_file_num:{finished_task_file_num}"
            f":total_core_hour:{total_core_hour:.3f}")
    
    model_devi_dirs = subprocess.run([f"ls -d -1 {target_folder}/iter.??????/01.model_devi/", ],
        shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n')
    # print(model_devi_dirs)
    for dir in model_devi_dirs:
        abs_dir = os.path.abspath(dir)
        stage = os.path.basename(os.path.dirname(dir))
        # print(dir)
        model_devi_time_logs = subprocess.run([f"grep -H --text 'wall time' {dir}/task.*/log.lammps", ],
            shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n')
        upload_task_dir_num = subprocess.run([f"ls -1 -d {dir}/task.* |wc -l", ], 
            shell=True, stdout=subprocess.PIPE).stdout.decode().strip('\n')
        total_core_sec = float(0)
        finished_task_file_num = len(model_devi_time_logs)
        # assume model_devi lammps job running on GPUs , set paral_cores==1
        paral_cores = 1 
        for log in model_devi_time_logs:
            # log example:
            #   .//iter.000002/01.model_devi//task.018.000075/log.lammps:Total wall time: 0:00:39
            # print(log)
            file_path, text1, hour, min, sec = log.split(':') # pylint: disable=unused-variable
            abs_file_path =  os.path.abspath(file_path)
            wall_time_sec = 3600*int(hour) + 60*int(min) + 1*int(sec)
            total_core_sec += wall_time_sec * paral_cores
        total_core_hour = total_core_sec / 3600

        print(f"{stage}:{abs_dir}"
            f":paral_cores:{paral_cores}"
            f":upload_task_dir_num:{upload_task_dir_num}"
            f":finished_task_file_num:{finished_task_file_num}"
            f":total_core_hour:{total_core_hour:.3f}")

    fp_dirs = subprocess.run([f"ls -d -1 {target_folder}/iter.??????/02.fp/", ],
        shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n')
    for dir in fp_dirs:
        abs_dir = os.path.abspath(dir)
        stage = os.path.basename(os.path.dirname(dir))
        fp_time_logs = subprocess.run([f"grep -H --text 'CPU time' {dir}/task.*/OUTCAR", ],
            shell=True,stdout=subprocess.PIPE).stdout.decode().strip().split('\n')
        upload_task_dir_num = subprocess.run([f"ls -1 -d {dir}/task.* |wc -l", ], 
            shell=True, stdout=subprocess.PIPE).stdout.decode().strip('\n')
        total_core_sec = float(0)
        finished_task_file_num = len(fp_time_logs)
        for log in fp_time_logs:
            # log example:
            #   .//iter.000002/02.fp//task.018.000048/OUTCAR:                  Total CPU time used (sec):      288.395
            file_path, text1, sec = log.split(':')
            abs_file_path = os.path.abspath(file_path)
            wall_time_sec = float(sec)
            paral_cores = subprocess.run([fr"head -n 1000 {abs_file_path} | grep 'running on' | sed -n -E -e 's|running on\s+([0-9]+)+\s.*|\1|p' ", ],
                shell=True,stdout=subprocess.PIPE).stdout.decode().strip()
            total_core_sec += wall_time_sec * int(paral_cores)
        total_core_hour = total_core_sec /3600

        print(f"{stage}:{abs_dir}"
            f":paral_cores:{paral_cores}"
            f":upload_task_dir_num:{upload_task_dir_num}"
            f":finished_task_file_num:{finished_task_file_num}"
            f":total_core_hour:{total_core_hour:.3f}")

if __name__=='__main__':
    stat_time(target_folder="./")