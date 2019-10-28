#!/bin/bash

model_devi_paral_cores=1

if [[ -a time.log ]]
then
	rm time.log
fi
for train_dir in `ls -d -1 iter.??????/00.train/`;do
sec=0
tothour=0
upload_task_dir_num=0
recycle_task_file_num=0
# echo $train_dir
upload_task_dir_num=$(ls -1 -d $train_dir/??? |wc -l)
if [[ -a train_time.log ]]
then
	rm train_time.log
fi
grep -H --text 'wall time'  $train_dir/???/train.log > train_time.log
recycle_task_file_num=$(wc -l < train_time.log)
 while read line; do
mysec=$(echo "$line" |cut -d: -f4 |sed 's/s\| //g')
sec=$(echo "$mysec + $sec" | bc)
   done < train_time.log
# echo $hour:$min:$sec
tothour=$(echo "scale=3; $sec/3600"|bc)
echo "00.train:$(realpath $train_dir):paral_cores:GPUV100:upload_task_dir_num:$upload_task_dir_num:recycle_task_file_num:$recycle_task_file_num:total core hour:$tothour" | tee -a time.log
done

for model_devi_dir in `ls -d -1 iter.??????/01.model_devi/`;do
sec=0
min=0
hour=0
tothour=0
upload_task_dir_num=0
recycle_task_file_num=0
# echo $model_devi_dir
upload_task_dir_num=$(ls -1 -d $model_devi_dir/task.* |wc -l)
if [[ -a model_devi_time.log ]]
then
	rm model_devi_time.log
fi
grep -H --text 'wall'  $model_devi_dir/task.*/log.lammps > model_devi_time.log
recycle_task_file_num=$(wc -l < model_devi_time.log)
 while read line; do
mysec=$(echo "$line" |cut -d: -f5)
sec=$(echo "$mysec + $sec" | bc)
mymin=$(echo "$line" |cut -d: -f4)
min=$(echo "$mymin + $min" | bc)
myhour=$(echo "$line" |cut -d: -f3)
hour=$(echo "$myhour + $hour" | bc)
   done < model_devi_time.log
# echo $hour:$min:$sec
tothour=$(echo "scale=3; ($hour*3600+$min*60+$sec)*$model_devi_paral_cores/3600"|bc)
echo "01.model_devi:$(realpath $model_devi_dir):paral_cores:$model_devi_paral_cores:upload_task_dir_num:$upload_task_dir_num:recycle_task_file_num:$recycle_task_file_num:total core hour:$tothour" | tee -a time.log
done

for fp_dir in `ls -d -1 iter.??????/02.fp/`;do
core_sec=0
tothour=0
upload_task_dir_num=0
recycle_task_file_num=0
# echo $fp_dir
upload_task_dir_num=$(ls -1 -d $fp_dir/task.* |wc -l)
if [[ -a fp_time.log ]]
then
	rm fp_time.log
fi
grep -H --text 'CPU time' $fp_dir/task.*/OUTCAR > fp_time.log
recycle_task_file_num=$(wc -l < fp_time.log)
 while read line;do
mysec=$(echo "$line" |cut -d: -f3 |sed 's| ||g')
file_name=$(echo "$line" | cut -d: -f1)
fp_paral_cores=$(grep 'total cores' $file_name |grep -o '[0-9]*')
core_sec=$(echo "$mysec * $fp_paral_cores + $core_sec" | bc)
   done < fp_time.log
tothour=$(echo "scale=3; $core_sec/3600"|bc)
echo "02.fp:$(realpath $fp_dir):paral_cores:$fp_paral_cores:upload_task_dir_num:$upload_task_dir_num:recycle_task_file_num:$recycle_task_file_num:total core hour:$tothour" | tee -a time.log
done
wc -l iter.??????/02.fp/*out> candi_fail_accu.log
