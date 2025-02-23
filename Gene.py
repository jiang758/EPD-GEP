
import random
import math
import Define

class Gene():

    def __init__(self, genome, homeotic_flag,DC_flag):
        '''
        genome:         基因的 Genome
        homeotic_flag:  Flag值 True:为homeotic基因 False: 为传统基因
        head            基因头部
        tail            基因尾部
        '''
        self.genome = genome
        self.homeotic_flag = homeotic_flag
        self.DC_flag = DC_flag
        self.head = []
        self.tail = []
        self.DC = []
        self.RandomSet=[]
        self.nfunctions=[]

    def selectFunction(self):
        total_arity = sum([arity for _, arity in self.nfunctions.items()])
        probabilities = {func: arity / total_arity for func,arity in self.nfunctions.items()}

        selected_func = random.choices(list(self.nfunctions.keys()), weights=list(probabilities.values()))[0]
        return selected_func

    def initRand(self, headLength, numGenes ,numRandom=10,RandomRangeStart=-1,RandomRangeEnd=1,RandomPRECISION=2):

        max_arity = self.genome.functions["max_arity"]
        functions = list(self.genome.functions.keys())
        functions.pop(functions.index("max_arity")) # functions中剔除max_arity
        arithmetic_set = {}
        for key in self.genome.arithmetic_set:
            if key != "max_arity":
                arithmetic_set[key] = self.genome.arithmetic_set[key][1]
        self.nfunctions = arithmetic_set

        # 产生随机数集合
        if self.DC_flag:
            # 产生DC域
            DC_set = list(range(numRandom))
            self.DC=[random.choice(DC_set) for i in range(headLength * (max_arity - 1) + 1)]
            # 产生DC的随机数集合
            self.RandomSet=[round(random.uniform(RandomRangeStart,RandomRangeEnd),RandomPRECISION) for i in range (numRandom)]

        # 产生传统基因和连接基因
        if self.homeotic_flag:
            terminals = list(range(numGenes))       #生成homeotic的terminals、
            if numGenes > list(self.genome.genicTerminals)[-1]:
                self.genome.genicTerminals = range(numGenes)
        else:
            terminals = self.genome.terminals
        #self.head = [random.choice(functions + terminals) for i in range(headLength)]
        #half_head_length = headLength // 10000
        #self.head = [self.selectFunction() if i < half_head_length else random.choice(functions + terminals) for i in
                    # range(headLength)]
        #self.tail = [random.choice(terminals) for _ in range(headLength * (max_arity - 1) + 1)]


        self.head = [self.selectFunction() for i in range(headLength)]
        self.tail = [random.choice(terminals) for i in range(headLength * (max_arity - 1) + 1)]
        return self
        #基因首个为运算符
        # self.head[0] = random.choice(functions)
        #return self

    def initRand2(self, headLength, numGenes, numRandom=10, RandomRangeStart=-1, RandomRangeEnd=1, RandomPRECISION=2):

        max_arity = self.genome.functions["max_arity"]
        functions = list(self.genome.functions.keys())
        functions.pop(functions.index("max_arity"))  # functions中剔除max_arity

        # 产生随机数集合
        if self.DC_flag:
            # 产生DC域
            DC_set = list(range(numRandom))
            self.DC = [random.choice(DC_set) for i in range(headLength * (max_arity - 1) + 1)]
            # 产生DC的随机数集合
            self.RandomSet = [round(random.uniform(RandomRangeStart, RandomRangeEnd), RandomPRECISION) for i in
                              range(numRandom)]

        # 产生传统基因和连接基因
        if self.homeotic_flag:
            terminals = list(range(numGenes))  # 生成homeotic的terminals、
            if numGenes > list(self.genome.genicTerminals)[-1]:
                self.genome.genicTerminals = range(numGenes)
        else:
            terminals = self.genome.terminals
        self.head = [random.choice(functions + terminals) for i in range(headLength)]
        self.tail = [random.choice(terminals) for i in range(headLength * (max_arity - 1) + 1)]
        # 基因首个为运算符
        # self.head[0] = random.choice(functions)
        return self

    def eval(self, inputs):
        '''
        对精英进行赋值运算评估
        :param inputs:      参数输入值
        :return:
        '''
        elementLayers = self.orderStack()
        evalStack = self.tranformRPN(elementLayers)
        return self.evalStack(evalStack, inputs)


    #例子【'+':2,'-':2,'*':2,'a':0,'a':0,'a':0,'a':0】
    #首先，取出一个符号'+',即arity初始为：arity=1
    #然后，保存该符号的arity数值。由'+':2得arity=2，同时将取出得'+'在elementLayers保存
    #循环，此时arity=2，取出2个符号即'-','*'。
    #根据'-','*'得到max_arity数值，max_arity=4
    #以此类推直到max_arity=0。
    #这其实就模拟树，一层一层得取出。
    #【【'+':2】，【'-':2,'*':2】,【'a':0,'a':0,'a':0,'a':0】】
    def orderStack(self):
        elementLayers = []
        gene = self.head + self.tail
        arity = 1
        while arity > 0:
            elementList = []
            for i in range(arity):
                elementList.append(gene.pop(0))
            arity = 0
            for i_element in elementList:
                arity += self.genome.symbols()[i_element]
            elementLayers.append(elementList)
        # print(elementLayers)
        # 树结构
        return elementLayers

    #用递归实现逆波兰式的装转变
    def tranformRPN(self, elementLayers, layerIndex = 0):
        '''
        :param elementLayers:    orderStack输出的层次结构
        :param layerIndex:       几层用于递归调用
        :return:
        '''
        stack = []
        for i in range(self.genome.symbols()[elementLayers[layerIndex][0]]):
            stack =  stack + self.tranformRPN(elementLayers, layerIndex + 1)
        stack.append(elementLayers[layerIndex].pop(0))
        return stack

    #inputs=['a':1,'b':3]已赋值，取值运算
    def evalStack(self, stack, inputs):
        '''
        :param stack:       为逆波兰式
        :param inputs:      参数输入值
        :return:returnStack     保存最后的运算结果
        '''
        DC=self.DC
        RandomSet=self.RandomSet
        countDC=0
        returnStack = []
        for symbol in stack:
            if symbol=='?':
                returnStack.append(RandomSet[DC[countDC]])
                countDC+=1
                pass
            elif symbol in inputs.keys():
                returnStack.append(inputs[symbol])
            else:
                isreturn=self.genome.functions[symbol][0](returnStack)
                if isreturn=='error':
                    return math.inf
        return returnStack[0]

    def simplifyGene(self,inputs):
        elementLayers = self.orderStack()
        evalStack = self.tranformRPN(elementLayers)
        return self.simplifyEvalStack(evalStack, inputs)

    def simplifyEvalStack(self,stack,inputs):
        DC=self.DC
        RandomSet=self.RandomSet
        countDC=0
        returnStack = []
        for symbol in stack:
            if self.homeotic_flag:
                if symbol in self.genome.genicTerminals:
                    returnStack.append(inputs[symbol])
                else:
                    self.genome.simplify_all_function_set[symbol](returnStack)
            else:
                if symbol=='?':
                    returnStack.append(RandomSet[DC[countDC]])
                    countDC+=1
                    pass
                elif symbol in self.genome.terminals:
                    returnStack.append(symbol)
                else:
                    self.genome.simplify_all_function_set[symbol](returnStack)
        return returnStack[0]

    def replicate(self):
        '''
        创建新的基因，并对基因内的元素进行复制
        :return: newGene    返回新的基因
        '''
        newGene = Gene(self.genome, self.homeotic_flag,self.DC_flag)
        newGene.head = self.head[:]
        newGene.tail = self.tail[:]
        newGene.DC =self.DC[:]
        newGene.RandomSet=self.RandomSet[:]

        return newGene

    def printGene(self):
        if self.DC_flag:
            print(''.join(self.head+self.tail),"",self.DC,"", self.RandomSet,end="")
        else:
            if self.homeotic_flag:
                print(self.head+self.tail,end=" ")
            else:
                print(''.join(self.head+self.tail),end="" )
