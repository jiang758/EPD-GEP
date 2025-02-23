# 开发人员：九重！
# 开发时间：2022-03-3019:25
# 开发名称:生成数据.py
# 开发应用:PyCharm
import math
import random

file=open("test.txt",'a')


num=10
v_num=3
file.write(str(num)+"  "+str(v_num)+"\n")
for i in range(num):
    x1=round(random.uniform(0,10),5)
    x2=round(random.uniform(0,10),5)
    x3=round(random.uniform(0,10),5)
    while not (-x1+x2*x3>=0 and x2>0):
        x1=round(random.uniform(0,10),5)
        x2=round(random.uniform(0,10),5)
        x3=round(random.uniform(0,10),5)
    function=x1+math.log(x2) + math.log(-x1+x2*x3)
    y=round(function,6)
    file.write(str(x1)+"  "+str(x2)+"  "+str(x3)+"  "+str(y)+"\n")