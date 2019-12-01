import numpy as np
import os,re,glob

def select_outcar(keyword):
    task=glob.glob(keyword+"-*")
    for ii in task:
        #os.chdir(ii)
        ret ="new OUTCAR\n"
        new_outcar=os.path.join(ii,'new_OUTCAR')
        outcar=os.path.join(ii,'OUTCAR')
        count=0
        with open(outcar, 'r') as fp:
            lines = fp.read().split('\n')
            for ii in lines:
                if ('Elapsed' in ii) or ('free  energy   TOTEN' in ii) or ('ions per type' in ii) or ('volume of cell' in ii):
                    ret=ret+ii+"\n"
                if ('direct lattice vectors' in ii):
                    count = 42
                if count>0:
                    ret = ret+ii+"\n"
                    count =count-1
        with open(outcar,'w') as new_fp:
            new_fp.write(str(ret))
    print(task)


select_outcar("struct")


