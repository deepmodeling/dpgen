import numpy as np
import requests

def voigt_to_stress(inpt) :
    ret = np.zeros((3,3))
    ret[0][0] = inpt[0]
    ret[1][1] = inpt[1]
    ret[2][2] = inpt[2]
    ret[0][1] = inpt[3]
    ret[0][2] = inpt[4]
    ret[1][2] = inpt[5]
    ret[2][0] = ret[0][2]
    ret[1][0] = ret[0][1]
    ret[2][1] = ret[1][2]
    return ret

def insert_data(task,task_type,username,file_name):
    assert task in ['eos','elastic','surf']
    assert task_type in ['vasp','deepmd']
    #check the corresponding of expr_type and data_type
    url='http://115.27.161.2:5000/insert_test_data?username=%s&expr_type=%s&data_type=%s' % (username,task_type,task)
    res = requests.post(url, data=open(file_name).read())
    print('Successful upload!')
