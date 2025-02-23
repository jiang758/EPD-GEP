
import random
import numpy as np
from sympy import *
from Gene import Gene

class Chromosome():
    def __init__(self,homeotic_flag,DC_flag,genome,link_genome):
        '''
        genes           基因：用于保存传统基因
        homeotics       连接基因：用于保存连接基因
        '''
        self.genes = []
        self.homeotics = []

        self.homeotic_flag=homeotic_flag
        self.DC_flag=DC_flag

        self.genome=genome
        self.link_genome=link_genome
        self.fitness = 0  # 添加fitness属性，并初始化为0

    #染色体初始化（基因、连接函数）
    def initRand(self, numGenes, numHomeotics, headLength, homeoticHeadLength
                 ,numRandom, RandomRangeStart, RandomRangeEnd,RandomPRECISION):
        '''
        初始化染色体
        :param numGenes:        基因个数
        :param numHomeotics:    连接基因个数
        :param headLength:      基因头部长度
        :param homeoticHeadLength:      连接基因头部长度
        :param genome:                  基因 Genome
        :param link_genome:             连接基因 Genome
        :return self:
        genes           基因：用于保存传统基因
        homeotics       连接基因：用于保存连接基因
        '''
        #self.genes = [Gene(genome=self.genome,homeotic_flag=False,DC_flag=self.DC_flag)
         #                 .initRand(headLength=headLength,numGenes=numGenes,numRandom=numRandom,
          #                          RandomRangeStart=RandomRangeStart,RandomRangeEnd=RandomRangeEnd,
           #                         RandomPRECISION=RandomPRECISION) for i in range(numGenes)]
        self.genes = [Gene(genome=self.genome, homeotic_flag=False, DC_flag=self.DC_flag)
                          .initRand(headLength=headLength, numGenes=numGenes, numRandom=numRandom,
                                    RandomRangeStart=RandomRangeStart, RandomRangeEnd=RandomRangeEnd,
                                    RandomPRECISION=RandomPRECISION) if i < numGenes *1.0
                      else Gene(genome=self.genome, homeotic_flag=False, DC_flag=self.DC_flag)
            .initRand2(headLength=headLength, numGenes=numGenes, numRandom=numRandom,
                       RandomRangeStart=RandomRangeStart, RandomRangeEnd=RandomRangeEnd,
                       RandomPRECISION=RandomPRECISION) for i in range(numGenes)]

        if  self.homeotic_flag:
            self.homeotics = [Gene(genome=self.link_genome, homeotic_flag=True,DC_flag=False)
                                  .initRand(headLength=homeoticHeadLength,numGenes=numGenes,numRandom=numRandom,
                                            RandomRangeStart=RandomRangeStart,RandomRangeEnd=RandomRangeEnd,
                                            RandomPRECISION=RandomPRECISION) for i in range(numHomeotics)]
        return self

        # 染色体评估

    def eval(self, inputs):
        '''
        :param inputs:      样本输入值
        :return:            returnValues
        '''
        returnValues = []
        #基因赋值运算、结果保存在returnValues

        for gene in self.genes:
            #将终结值赋予对应的数值,如 a:1 ,b:2
            evalInputs = {terminal: value for (terminal, value) in list(zip(gene.genome.terminals, inputs))}
            value=gene.eval(evalInputs)
            returnValues.append(gene.eval(evalInputs))

        #连接函数运算、最终结果保存在returnValues
        homeoticReturnValues = []
        if len(self.homeotics) > 0:
            outputs = {index: output for (index, output) in list(zip(range(0, len(self.genes)), returnValues))}
            for homeotic in self.homeotics:
                homeoticReturnValues.append(homeotic.eval(outputs))
            returnValues = homeoticReturnValues
        else:
            sum=0
            for value in returnValues:
                sum+=value
            returnValues=[sum]
        return returnValues

    def replicate(self):
        '''
        基因复制
        newGenes        保存新的基因
        newHomeotics    保存新的连接基因
        newChromosome   创建新的染色体，并赋值
        :return: newChromosome 返回新的染色体
        '''
        newGenes = []
        newHomeotics = []
        # 基因复制
        for gene in self.genes:
            newGenes.append(gene.replicate())
        # 连接基因复制
        for homeotic in self.homeotics:
            newHomeotics.append(homeotic.replicate())
        # 创建新的染色体
        newChromosome = Chromosome(self.homeotic_flag,self.DC_flag,self.genome,self.link_genome)
        newChromosome.genes = newGenes
        newChromosome.homeotics = newHomeotics
        return newChromosome

    def mutation(self, rate, homeoticRate,DcRate,randomSetRate,RandomRangeStart,RandomRangeEnd,RandomPRECISION):
        '''

        :param rate:    基因突变率
        :param homeoticRate: 连接基因突变率
        :param DcRate:  Dc域突变率
        :param randomSetRate:   随机数突变率
        :param RandomRangeStart: 随机数开始
        :param RandomRangeEnd: 随机数结束
        :param RandomPRECISION: 随机数精度
        :return:
        '''


        for gene in self.genes:
            head = gene.head
            tail = gene.tail
            functions = list(gene.genome.functions.keys())
            functions.pop(functions.index("max_arity"))
            # 基因头部突变
            for i in range(len(head)):
                if random.random() < rate:
                    head[i] = random.choice(functions + gene.genome.terminals)

            # 基因尾部突变
            for i in range(len(tail)):
                if random.random() < rate:
                    tail[i] = random.choice(gene.genome.terminals)
            # Dc域：随机常数突变
            if self.DC_flag:
                DcSpecific=gene.DC
                for i in range(len(DcSpecific)):
                    if random.random() < DcRate:
                        DC_set = list(range(len(gene.RandomSet)))
                        DcSpecific[i]=random.choice(DC_set)
                randomSet = gene.RandomSet
                for i in range(len(randomSet)):
                    if random.random() < randomSetRate:
                        randomSet[i]=round(random.uniform(RandomRangeStart,RandomRangeEnd),RandomPRECISION)

        #连接基因的突变为 rate = homeoticRate*rate
        if self.homeotic_flag:
            rate = homeoticRate
            for homeotic in self.homeotics:
                head = homeotic.head
                tail = homeotic.tail
                functions = list(homeotic.genome.functions.keys())
                functions.pop(functions.index("max_arity"))
                # 连接基因头部突变，突变的选择有些不同：functions + list(range(len(self.genes)))
                for i in range(len(head)):
                    if random.random() < rate:
                        head[i] = random.choice(functions + list(range(len(self.genes))))
                # 连接基因尾部部突变
                for i in range(len(tail)):
                    if random.random() < rate:
                        tail[i] = random.randint(0, len(self.genes) - 1)


    def ISTransposition(self, rate,elementsLength,homeoticRate,homeoticelementsLength,DcRate,DcElementsLength):
        '''
                IS转座是从种群中随机选择一个染色体，然后随机从IS转座长度中选择IS长度，
        然后在染色体中选择IS长度的基因片段，并随机选择基因，插入到除基因首元素之外的头部部分中。
        :param rate: 基因IS概率
        :param elementsLength: 基因IS长度
        :param homeoticRate: 连接基因IS概
        :param homeoticelementsLength: 连接基因IS长度
        :param DcRate: Dc域IS概率
        :param DcElementsLength: Dc域IS长度
        :return:
        '''

        # 出现IS转座
        if random.random() < rate:
            seq = []
            for gene in self.genes:
                seq = seq + gene.head + gene.tail

            #   随机选择一个基因
            gene = random.choice(self.genes)
            headLength=len(gene.head)
            length = random.choice(elementsLength)

            #   选择截取一小断保存在 ISSeq 中
            start = random.randint(0, len(seq) - 1)
            stop = start + length
            ISSeq = seq[start:stop]

            #   选择一个插入点，进行插入，并且去掉超过长度的部分
            insertPoint = random.randint(1, headLength - 1)
            for i in range(len(ISSeq)):
                gene.head.insert(insertPoint + i, ISSeq[i])
            gene.head = gene.head[:headLength]

        #连接基因的IS转座为 rate = homeoticRate*rate
        if self.homeotic_flag:
            rate = homeoticRate
            if random.random() < rate:
                seq = []
                for homeotic in self.homeotics:
                    seq = seq + homeotic.head + homeotic.tail
                #   随机选择一个基因
                homeotic = random.choice(self.homeotics)
                headLength = len(homeotic.head)
                length=random.choice(homeoticelementsLength)

                #   选择截取一小断保存在 ISSeq 中
                start = random.randint(0, len(seq) - 1)
                stop = start + length
                ISSeq = seq[start:stop ]
                insertionPoint = random.randint(1, headLength - 1)

                #   选择一个插入点，进行插入，并且去掉超过长度的部分
                for i in range(len(ISSeq)):
                    homeotic.head.insert(insertionPoint + i, ISSeq[i])
                homeotic.head = homeotic.head[:headLength]

        if self.DC_flag:
            seq = []
            rate = DcRate
            if random.random() < rate:
                for gene in self.genes:
                    seq = seq+gene.DC
                #   随机选择一个基因
                gene = random.choice(self.genes)
                TailLength=len(gene.DC)
                length = random.choice(DcElementsLength)

                #   选择截取一小断保存在 ISSeq 中
                start = random.randint(0, len(seq) - 1)
                stop = start + length
                ISSeq = seq[start:stop]

                #   选择一个插入点，进行插入，并且去掉超过长度的部分
                insertPoint = random.randint(0, TailLength - 1)
                for i in range(len(ISSeq)):
                    gene.DC.insert(insertPoint + i, ISSeq[i])
                gene.DC = gene.DC[:TailLength]


    def RISTransposition(self, rate,elementsLength, homeoticRate,homeoticelementsLength):
        '''
        :param rate: 基因RIS转座概率
        :param elementsLength:  基因RIS转座长度
        :param homeoticRate:    连接基因RIS转座概率
        :param homeoticelementsLength: 连接基因RIS转座长度
        :return:
        '''
        # 进行RIS转座

        if random.random() < rate:
            seq = []
            for gene in self.genes:
                seq = seq + gene.head + gene.tail

            #   随机选择一个基因
            gene = random.choice(self.genes)
            headLength = len(gene.head)
            length=random.choice(elementsLength)

            functions = list(gene.genome.functions.keys())
            functions.pop(functions.index("max_arity"))

            #   选择截取一小断保存在 RISSeq 中
            start = random.randint(0, len(seq) - 1)

            #直到发现一个函数为止
            while seq[start] not in functions:
                if start == len(seq) - 1:
                    break
                start += 1
            stop = start + length
            RISSeq = seq[start:stop + 1]

            #进行插入，并且去掉超过长度的部分
            for i in range(len(RISSeq)):
                gene.head.insert(i, RISSeq[i])
            gene.head = gene.head[:headLength]

        #连接基因的RIS转座为 rate = homeoticRate*rate
        if self.homeotic_flag:
            rate = homeoticRate
            if random.random() < rate:
                seq = []
                for homeotic in self.homeotics:
                    seq = seq + homeotic.head + homeotic.tail

                #   随机选择一个基因
                homeotic = random.choice(self.homeotics)
                headLength = len(homeotic.head)
                length=random.choice(homeoticelementsLength)

                functions = list(homeotic.genome.functions.keys())
                functions.pop(functions.index("max_arity"))

                #   选择截取一小断保存在 RISSeq 中
                start = random.randint(0, len(seq) - 1)
                while seq[start] not in functions:
                    if start == len(seq) - 1:
                        break
                    start += 1
                stop = start + length

                #进行插入，并且去掉超过长度的部分
                RISSeq = seq[start:stop ]
                for i in range(len(RISSeq)):
                    homeotic.head.insert(i, RISSeq[i])
                homeotic.head = homeotic.head[:headLength]


    def geneTransposition(self, rate):
        '''
        基因转座仅仅改变基因在同一染色体中的位置
        :param rate:            基因转座概率
        :param homeoticRate:    连接基因转座概率
        :return:
        '''

        if random.random() < rate:
            self.genes[random.randint(0,len(self.genes)-1)] = random.choice(self.genes).replicate()



    def onePointRecombination(self, rate, otherChromosome):
        '''
        进行单点重组的时候，父代染色体相互配对并在相同的位置切断，两个染色体相互交换重组点之后的部分
        :param rate:                单点重组概率
        :param homeoticRate:        连接基因单点重组概率
        :param otherChromosome:     染色体
        :return:
        seq             当前染色体将所有基因进行连接保存
        otherSeq        不同于当前染色体将所有基因进行连接保存
        recombinationPoint       重组结点
        '''
        # 进行单点重组

        if random.random() < rate:
            seq = []
            otherSeq = []
            for gene in self.genes:
                seq = seq + gene.head + gene.tail
            for otherGene in otherChromosome.genes:
                otherSeq = otherSeq + otherGene.head + otherGene.tail
            if len(seq) != len(otherSeq):
                return
            # 进行单点重组（交换）
            recombinationPoint = random.randint(0, len(seq) - 1)
            seq[recombinationPoint:], otherSeq[recombinationPoint:] = otherSeq[recombinationPoint:], seq[recombinationPoint:]

            # 单点重组后,当前染色体进行保存
            for gene in self.genes:
                gene.head = seq[:len(gene.head)]
                del seq[0:len(gene.head)]
                gene.tail = seq[:len(gene.tail)]
                del seq[0:len(gene.tail)]

            # 单点重组后,不同于当前染色体进行保存
            for otherGene in otherChromosome.genes:
                otherGene.head = otherSeq[:len(otherGene.head)]
                del otherSeq[0:len(otherGene.head)]
                otherGene.tail = otherSeq[:len(otherGene.tail)]
                del otherSeq[0:len(otherGene.tail)]

        #连接基因的单点重组为 rate = homeoticRate*rate

    def twoPointRecombination(self, rate, otherChromosome):
        '''
        进行两点重组的时候，父代染色体相互配对，在染色体中随机选择两个点，将染色体切断。
        两个染色体相互交换重组点之间的部分，形成两个新的子代染色体
        :param rate:                两点重组概率
        :param homeoticRate:        连接基因两点重组概率
        :param otherChromosome:     不同于当前的染色体
        :return:
        seq             当前染色体将所有基因进行连接保存
        otherSeq        不同于当前染色体将所有基因进行连接保存
        recombinationPoint       重组结点
        otherPoint               另一个重组结点
        '''
        # 进行两点重组

        if random.random() < rate:
            seq = []
            otherSeq = []
            for gene in self.genes:
                seq = seq + gene.head + gene.tail
            for otherGene in otherChromosome.genes:
                otherSeq = otherSeq + otherGene.head + otherGene.tail
            if len(seq) != len(otherSeq):
                return
            # 进行两点重组（交换）
            recombinationPoint = random.randint(0, len(seq) - 1)
            otherPoint = random.randint(recombinationPoint, len(seq) - 1)
            seq[recombinationPoint:], otherSeq[recombinationPoint:] = otherSeq[recombinationPoint:], seq[recombinationPoint:]
            seq[:otherPoint], otherSeq[:otherPoint] = otherSeq[:otherPoint], seq[:otherPoint]

            # 两点重组后,当前染色体进行保存
            for gene in self.genes:
                gene.head = seq[:len(gene.head)]
                del seq[0:len(gene.head)]
                gene.tail = seq[:len(gene.tail)]
                del seq[0:len(gene.tail)]

            # 两点重组后,不同于当前染色体进行保存
            for otherGene in otherChromosome.genes:
                otherGene.head = otherSeq[:len(otherGene.head)]
                del otherSeq[0:len(otherGene.head)]
                otherGene.tail = otherSeq[:len(otherGene.tail)]
                del otherSeq[0:len(otherGene.tail)]

        #连接基因的两点重组为 rate = homeoticRate*rate

    def geneRecombination(self, rate, otherChromosome):
        '''
        在 GEP 的第三种重组中，两个染色体中的整个基因相互交换，
        形成的两个子代染色体含有来自两个父体的基因。
        :param rate:            基因重组概率
        :param homeoticRate:    连接基因重组概率
        :param otherChromosome: 不同于当前的染色体
        :return:
        '''

        # 进行基因重组
        if random.random() < rate:
            if len(self.genes) != len(otherChromosome.genes):
                return
            recombinationPoint = random.randint(0, len(self.genes) - 1)
            self.genes[recombinationPoint:recombinationPoint+1], otherChromosome.genes[recombinationPoint:recombinationPoint+1] = otherChromosome.genes[recombinationPoint:recombinationPoint+1], self.genes[recombinationPoint:recombinationPoint+1]


    def createNeutralGene(self,inputsOutputs,neutralValue,numRandom,RandomRangeStart,RandomRangeEnd,RandomPRECISION):
        returnMean=-1
        numGenes=len(self.genes)
        headLength=len(self.genes[0].head)
        while(returnMean!=neutralValue):
            gene =Gene(self.genome,homeotic_flag=False,DC_flag=self.DC_flag).initRand(headLength=headLength,numGenes=numGenes,
                                                                                        numRandom=numRandom,RandomRangeStart=RandomRangeStart,                                                                    RandomRangeEnd=RandomRangeEnd,RandomPRECISION=RandomPRECISION)
            returnValues = []
            for inputs in inputsOutputs[0]:
                evalInputs = {terminal: value for (terminal, value) in list(zip(gene.genome.terminals, inputs))}
                returnValues.append(gene.eval(evalInputs))
            returnMean=np.mean(returnValues)
        return gene


    def updateNeutralGene(self,neutralGene):
        numGenes=len(self.genes)
        self.genes.append(neutralGene)
        if self.homeotic_flag:
            numHomeotics=len(self.homeotics)
            homeoticHeadLength=len(self.homeotics[0].head)
            self.homeotics = [Gene(genome=self.link_genome, homeotic_flag=True,DC_flag=False)
                                  .initRand(headLength=homeoticHeadLength,numGenes=numGenes,) for i in range(numHomeotics)]

        # newGenes = []
        # newHomeotics = []
        # # 基因复制
        # for gene in self.genes:
        #     newGenes.append(gene.replicate())
        # # 连接基因复制
        # for homeotic in self.homeotics:
        #     newHomeotics.append(homeotic.replicate())
        # # 创建新的染色体
        # newChromosome = Chromosome(self.homeotic_flag,self.DC_flag,self.genome,self.link_genome)
        # newChromosome.genes = newGenes
        # newChromosome.homeotics = newHomeotics
        # return newChromosome


  #def printChromosome(self):
     #   for gene in self.genes:
      #     gene.printGene()
       # for homeotic in self.homeotics:
        #    homeotic.printGene()
           #print()

    def printChromosome(self):
        for gene in self.genes:
            print(gene.head + gene.tail, end=':\n')
        for homeotic in self.homeotics:
            print(homeotic.head + homeotic.tail, end=';\n')
            print()


    def simplifyChromosome(self):

        returnSimplify = []
        for gene in self.genes:
            returnSimplify.append(gene.simplifyGene(inputs=None))
        #连接函数运算、最终结果保存在returnValues
        homeoticReturnValues = []
        if len(self.homeotics) > 0:
            outputs = {index: output for (index, output) in list(zip(range(0, len(self.genes)), returnSimplify))}
            for homeotic in self.homeotics:
                homeoticReturnValues.append(homeotic.simplifyGene(inputs=outputs))
            returnSimplify = homeoticReturnValues
        else:
            sum_str=''
            for simplify_str in  returnSimplify:
                sum_str="(" + sum_str + "+" + str(simplify_str) + ")"
            returnSimplify=[sum_str]
        print (simplify(returnSimplify[0]))
        # return simplify(returnSimplify[0])



