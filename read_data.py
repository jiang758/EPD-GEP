def read_data(path=".\\F1_0_training_data.txt"):
    with open(path) as f:
        num_data, num_variable = f.readline().split()
        num_data, num_variable = int(num_data), int(num_variable)

        Inputs=[]
        Outputs=[]
        for line in f:
            Input=[]
            Output=[]
            for i in range(len(line.split())):
                if i<num_variable:
                    Input.append(float(line.split()[i]))
                else:
                    Output.append(float(line.split()[i]))
            Inputs.append(Input)
            Outputs.append(Output)

        return num_data, num_variable, (Inputs,Outputs),


if __name__ == '__main__':
    print(read_data())
