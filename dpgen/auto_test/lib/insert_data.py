#coding=utf-8
# ls 0{1,2,5}*/*/*/vasp-k0.10/result|python insert_data.py
import sys
import requests
for line in sys.stdin:
  expr_type, element, structure, data_type, result = line.strip().split('/')
  expr_type = expr_type.split('.')[-1]
  data_type = data_type.split('-')[0].strip()
  res = requests.post('http://115.27.161.2:5000/insert_test_data?username=chenweijie&expr_type=%s&data_type=%s' % (data_type, expr_type), data=open(line.strip()).read())
