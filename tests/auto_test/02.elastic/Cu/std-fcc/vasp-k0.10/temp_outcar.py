import numpy as np
import os,re,glob

def select_outcar(keyword):
    task=glob.glob(keyword+"-*")
    for ii in task:
        #os.chdir(ii)
        ret ="new OUTCAR\n"
        outcar=os.path.join(ii,'OUTCAR')
        with open(outcar, 'r') as fp:
            lines = fp.read().split('\n')
            for ii in lines:
                if ('Elapsed time (sec):' in ii) or ('free  energy   TOTEN' in ii) or ('ions per type' in ii) or ('volume of cell' in ii) or ('in kB' in ii):

                    ret=ret+ii+"\n"
        with open(outcar,'w') as new_fp:
            new_fp.write(str(ret))
    print(task)


select_outcar("dfm")


