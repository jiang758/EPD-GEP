
import pandas as pd
def read_xls(path):

    df = pd.read_excel(path,sheet_name="Sheet1",usecols=[2,3,4],nrows=75)
    Inputs=df.values
    print(Inputs)
    print(Inputs.tolist())
    Inputs=Inputs.tolist()
    df = pd.read_excel(path,sheet_name="Sheet1",usecols=[14],nrows=75)
    Outputs=df.values
    # print(type(Outputs.tolist()))
    # print(Outputs)
    # print(list(Outputs.tolist()))
    Outputs=Outputs.tolist()
    num_data=75 # 数据个数
    num_variable=3  # 自变量个数
    return num_data, num_variable, (Inputs,Outputs)

if __name__ == '__main__':
    read_xls("./data.xls")
    # print(read_xls("./data.xls"))
