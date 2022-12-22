import numpy as np
import os

def _parse_calypso_input(var,input_path):

        if os.path.basename(input_path) != 'input.dat':
            input_path = os.path.join(input_path,'input.dat')
        if not os.path.exists(input_path):
            raise FileNotFoundError(input_path)

        f = open(input_path,'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if var in line:
                variable = line.split('=')[1].strip()
                return variable

def _parse_calypso_dis_mtx(numberofspecies,input_path):
        try:
            f = open(input_path,'r')
        except:
            f = open(os.path.join(input_path,'input.dat'),'r')
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            if '@DistanceOfIon' in line: 
                dis = []
                for i in range(int(numberofspecies)):
                    line = f.readline()
                    dis.append(line.split())
                f.close()
                break
        dis = np.array(dis)
        dis = dis.reshape((1,int(numberofspecies)**2))
        return dis[0][np.argmin(dis)]
