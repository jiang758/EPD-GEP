import random
import numpy as np

import matplotlib.pyplot as plt

from Chromosome import Chromosome
import Define


class Environment:

    def __init__(self):

        '''
        population          种群：用于保存染色体
        bestChromosome      最佳染色体：保存当代最佳染色体
        bestfitness         最佳染色体的适应度：保证当代最佳染色体的适应度
        homeoticRate        连接基因突变率
        mutationRate        突变
        ISTranspositionRate         IS
        RISTranspositionRate        RIS
        geneTranspositionRate       基因转座
        onePointRecombinationRate   单点重组
        twoPointRecombinationRate   两点重组
        geneRecombination           基因重组
        '''

        self.population = []
        self.bestChromosome = []
        self.bestfitness = []
        self.mutationRate = 0
        self.ISTranspositionRate = 0
        self.ISElementLength = 0
        self.RISTranspositionRate = 0
        self.RISElementLength = 0
        self.geneTranspositionRate = 0
        self.onePointRecombinationRate = 0
        self.twoPointRecombinationRate = 0
        self.geneRecombinationRate = 0

        self.homeoticMutationRate = 0
        self.homeoticISTranspositionRate = 0
        self.homeoticISElementLength = 0

        self.randomSetRate = 0
        self.DcMutationRate = 0
        self.DcISTranspositionRate = 0
        self.DcISElementLength = 0

        self.numRandom = 0
        self.RandomRangeStart = 0
        self.RandomRangeEnd = 0
        self.RandomPRECISION = 0

        self.neutral_flag = False
        self.TriggerFitnessCount = 20
        self.neutralValue = 0

    def setRates(self, mutationRate=0, ISTranspositionRate=0, ISElementLength=[], RISTranspositionRate=0,
                 RISElementLength=[], geneTranspositionRate=0,
                 onePointRecombinationRate=0, twoPointRecombinationRate=0, geneRecombinationRate=0,
                 homeoticMutationRate=0, homeoticISTranspositionRate=0,
                 homeoticISElementLength=[], homeoticRISTranspositionRate=0, homeoticRISElementLength=[],
                 randomSetRate=0, DcMutationRate=0, DcISTranspositionRate=0, DcISElementLength=[], numRandom=0,
                 RandomRangeStart=-1, RandomRangeEnd=1, RandomPRECISION=2, TriggerFitnessCount=20, neutralValue=0):
        '''
        homeoticRate        连接基因突变率
        mutationRate        突变
        ISTranspositionRate         IS
        RISTranspositionRate        RIS
        geneTranspositionRate       基因转座
        onePointRecombinationRate   单点重组
        twoPointRecombinationRate   两点重组
        geneRecombination           基因重组
        '''
        self.mutationRate = mutationRate
        self.ISTranspositionRate = ISTranspositionRate
        self.ISElementLength = ISElementLength
        self.RISTranspositionRate = RISTranspositionRate
        self.RISElementLength = RISElementLength
        self.geneTranspositionRate = geneTranspositionRate
        self.onePointRecombinationRate = onePointRecombinationRate
        self.twoPointRecombinationRate = twoPointRecombinationRate
        self.geneRecombinationRate = geneRecombinationRate

        self.homeoticMutationRate = homeoticMutationRate
        self.homeoticISTranspositionRate = homeoticISTranspositionRate
        self.homeoticISElementLength = homeoticISElementLength
        self.homeoticRISTranspositionRate = homeoticRISTranspositionRate
        self.homeoticRISElementLength = homeoticRISElementLength

        self.randomSetRate = randomSetRate
        self.DcMutationRate = DcMutationRate
        self.DcISTranspositionRate = DcISTranspositionRate
        self.DcISElementLength = DcISElementLength

        self.numRandom = numRandom
        self.RandomRangeStart = RandomRangeStart
        self.RandomRangeEnd = RandomRangeEnd
        self.RandomPRECISION = RandomPRECISION

        self.TriggerFitnessCount = TriggerFitnessCount
        self.neutralValue = neutralValue

    def init(self, populationSize, numGenes, numHomeotics, headLength, homeoticHeadLength, genome, link_genome,
             homeotic_flag=False, DC_flag=False, neutral_flag=False):
        '''
        初始化种群
        populationSize      种群大小
        numGenes            染色体基因个数
        numHomeotics        染色体连接基因个数
        headLength          基因头部长度
        homeoticHeadLength  连接基因头部长度
        genome              基因 Genome类
        link_genome         连接 link_genome 类
        '''
        self.neutral_flag = neutral_flag

        self.homeotic_flag = homeotic_flag
        self.DC_flag = DC_flag
        self.genome = genome
        self.link_genome = link_genome
        self.numHomeotics = numHomeotics
        self.headLength = headLength
        self.homeoticHeadLength = homeoticHeadLength


        self.population = [
            Chromosome(homeotic_flag, DC_flag, genome, link_genome).initRand(numGenes, numHomeotics, headLength,
                                                                             homeoticHeadLength
                                                                             , self.numRandom, self.RandomRangeStart,
                                                                             self.RandomRangeEnd, self.RandomPRECISION)
            for i in range(populationSize)]

    def generate_new_chromosome(self):

        newChromosome = Chromosome(self.homeotic_flag, self.DC_flag, self.genome, self.link_genome)
        newChromosome.initRand(6,self.numHomeotics, self.headLength, self.homeoticHeadLength, self.numRandom, self.RandomRangeStart,
                                                                             self.RandomRangeEnd, self.RandomPRECISION)
        return newChromosome



    def run(self, inputsOutputs, evalFunction, generationCount, M, DataCount,absoluteError_flagz, absoluteError_flag, roulette_flag):
        '''
                :param inputsOutputs: 保存输入与输出集合
                :param evalFunction: 适应度函数
                :param generationCount: 执行代数
                :param M: 选择范围
                :param DataCount: 数据个数
                :param absoluteError_flag:  绝对误差or相对误差
                :param roulette_flag:       赌盘or锦标赛
                :return:
                '''

        ResultBestFitness = float("-inf")
        ResultBestChromosome = []

        generation = 0
        AccFitnessCount = 0
        PreBestFitness = 0

        DX = 0
        zz = 0
        SSE = 0
        SST = 0
        X = []
        Y = []
        max_unchanged_count =0.2* generationCount  # 最大不更新次数为总迭代次数的5%
        while True:
            sum_fitness = 0
            NowBestFitness = float("-inf")
            NowBestChromosome = []
            generation += 1

            for i in range(len(self.population)):
                C_valueList = []
                i_chromosome = self.population[i]
                # 计算每个染色体对应的 C_value
                # print(inputsOutputs[0])
                for inputs in inputsOutputs[0]:
                    C_vlaue = i_chromosome.eval(inputs)
                    C_valueList.append(C_vlaue)
                # print(C_valueList)
                # evalFunction 进行适应度评估
                fitness = evalFunction(C_valueList, inputsOutputs[1], maxFitness_flag=False, M=M, DataCount=DataCount,
                                       absoluteError_flagz=absoluteError_flagz, absoluteError_flag=absoluteError_flag)
                # 每个染色体对应适应度
                # print("[", i, "] ", "   Fitness=", fitness)
                sum_fitness = sum_fitness + fitness
                # 保存当代种群最佳适应度的值,与对应的染色体
                if fitness > NowBestFitness:
                    NowBestFitness = fitness
                    NowBestChromosome = i_chromosome
                    # 当代最优所在种群位置
                    # Nowi=i
                    # print(Nowi)
                    # 保存当代种群最佳适应度对应建模值
                    NowBestC_Valuelist = C_valueList
                    # print(NowBestC_Valuelist)
                # 完美解停止执行
                if abs(fitness - evalFunction([], [], maxFitness_flag=True, M=M, DataCount=DataCount)) < 1e-5:
                    print("*" * 46, "  完美解  ", "*" * 46)
                    print("第", generation, "代", "BestFitness=", fitness)
                    i_chromosome.printChromosome()
                    i_chromosome.simplifyChromosome()
                    return i_chromosome
                i_chromosome.fitness = fitness  # 将适应度值赋给fitness属性
                self.population[i] = i_chromosome

            # 打印当代最佳染色体
            print("第", generation, "代，BestFitness=", NowBestFitness)
            #NowBestChromosome.printChromosome()
            #NowBestChromosome.simplifyChromosome()
            print("")

            x = generation
            y = NowBestFitness
            # X,Y轴坐标
            X.append(x)
            Y.append(y)

            if ResultBestFitness < NowBestFitness:
                ResultBestFitness = NowBestFitness
                ggeneration = generation
                ResultBestChromosome = NowBestChromosome
                ResultNowBestC_Valuelist = NowBestC_Valuelist
            if generation < ggeneration:
                ggeneration = generation

            # 繁衍代数大于generationCount退出, 并输出最佳解
            if generation >= generationCount:
                print(generationCount, "代后最佳解  ", " ", "最佳解在第", ggeneration, "代")
                print("result_fitness=  ", ResultBestFitness)
                ResultBestChromosome.printChromosome()
                ResultBestChromosome.simplifyChromosome()
                # for i in range(DataCount):
                #     zz+= np.array(inputsOutputs[1][i])
                # zz = zz / DataCount
                for i in range(DataCount):
                    dx = (np.array(inputsOutputs[1][i]) - np.array(ResultNowBestC_Valuelist[i])) * (
                                np.array(inputsOutputs[1][i]) - np.array(ResultNowBestC_Valuelist[i]))
                    DX = dx + DX
                    print("真实值", inputsOutputs[1][i], " ", "建模值：", ResultNowBestC_Valuelist[i])
                #     SSE+=(np.array(inputsOutputs[1][i]) - np.array(ResultNowBestC_Valuelist[i]))**2
                #     SST+=(np.array(inputsOutputs[1][i]) - zz)**2
                # R=1-SSE/SST
                # print(SSE,SST,R)
                DX = DX / DataCount
                SX = pow(DX, 0.5)
                print("方差：", DX, "\n标准差:", SX)
                # 散点图
                plt.scatter(X, Y, s=0.35)
                plt.xlabel('generationCount')
                plt.ylabel('NowBestFitness')
                plt.title('STATE')
                # plt.xticks(np.arange(0, 1000, step=50))
                plt.show()
                return ResultBestChromosome

            # 赌盘选择，复制，修饰
            selected = self.select(sum_fitness, NowBestChromosome, roulette_flag)
            self.replicate(selected)

            # 产生中性基因
            if self.neutral_flag:
                if (NowBestFitness > PreBestFitness * 0.90 and NowBestFitness < PreBestFitness * 1.10):
                    AccFitnessCount += 1
                    print("ACC:", AccFitnessCount, "  LenGenes:", len(self.population[0].genes))
                    if AccFitnessCount >= self.TriggerFitnessCount:
                        self.addNeutralGene(inputsOutputs, self.neutralValue)
                        AccFitnessCount = 0
                else:
                    AccFitnessCount = 0
            PreBestFitness = NowBestFitness
            self.modify()
            self.population[0] = NowBestChromosome

            # 根据最好适应度值是否更新进行染色体替换
            if generation - ggeneration >= max_unchanged_count:
                # 对适应度值进行排序
                self.population = sorted(self.population, key=lambda x: x.fitness, reverse=True)
                replace_count = int(0.2 * len(self.population))  # 需要替换的染色体数量
                #if NowBestFitness <= 0.7 * M * DataCount:
                    # 替换适应度值排序的最后replace_count个染色体
                for i in range(replace_count):
                    self.population[-(i + 1)] = self.generate_new_chromosome()  # 生成新的染色体
                #else:
                    # 替换适应度值小于0.7*70%*M*数据个数的染色体
                    #no77md = True
                   # for i in range(len(self.population)):
                  #      if self.population[i].fitness <= 0.7 * 0.7 * M * DataCount:
                    #        no77md = False
                    #        self.population[i] = self.generate_new_chromosome()
                            # 生成新的染色体
                   # if no77md:
                        # 替换适应度值排序的最后replace_count个染色体
                    #    for i in range(replace_count):
                     #       self.population[-(i + 1)] = self.generate_new_chromosome()  # 生成新的染色体


    def select(self, sum_fitness, NowBestChromosome, roulette_flag=True):

        '''
        :param sum_fitness:         总的适应度
        :param NowBestChromosome:   用于保存当代最佳染色体
        :return: selected           用于保存赌盘选择出的染色体
        :self.population[i] = (i_chromosome, fitness) 其中保存了染色体的适应度
        '''


        if roulette_flag:
            accumulator = 0  # 设置累加器
            roulette = []  # 设置轮盘
            selected = []

            # 制作赌盘，如果percentage  区间为0的
            for i in range(len(self.population)):
                percentage = self.population[i][1] / sum_fitness * 10
                if percentage != 0:
                    roulette.append(percentage + accumulator)
                    accumulator = accumulator + percentage
                else:
                    roulette.append(0)
            # 一共选择 len(self.population)-1 条染色体体，最后一条保存最佳染色体。
            for i in range(len(self.population) - 1):
                # 产生随机数
                ranNumber = random.uniform(0, 10)
                # 进行选择染色体
                for i_pos in range(len(self.population)):
                    if ranNumber <= roulette[i_pos]:
                        selected.append(self.population[i_pos])
                        break
            # 最后一条保存最佳染色体
            selected.append(NowBestChromosome)
            return selected
        else:
            selected = []
            for i in range(len(self.population) - 1):
                s1 = random.randint(0, len(self.population) - 1)
                s2 = random.randint(0, len(self.population) - 1)
                while (s1 == s2):
                    s2 = random.randint(0, len(self.population) - 1)
                if (self.population[s1].fitness < self.population[s2].fitness):
                    selected.append(self.population[s2])
                else:
                    selected.append(self.population[s1])
            selected.append(NowBestChromosome)
            return selected

    def replicate(self, chromosomes):
        '''
        进行染色体复制，保证更新在 self.population
        :param chromosomes:       赌盘后选择种群染色体
        :return:
        '''
        self.population = []
        for chromosome in chromosomes:
            self.population.append(chromosome.replicate())

    def modify(self):

        '''
        染色体修饰
        mutation:           突变
        ISTransposition     IS转座
        RISTransposition    RIS转座
        geneTransposition   基因转座
        onePointRecombination   单点重组
        twoPointRecombination   两点重组
        geneRecombination       基因重组
        :return:
        '''

        for i in range(len(self.population)):

            chromosome = self.population[i]

            chromosome.mutation(self.mutationRate, self.homeoticMutationRate, self.DcMutationRate, self.randomSetRate,
                                self.RandomRangeStart,
                                self.RandomRangeEnd, self.RandomPRECISION)

            chromosome.ISTransposition(self.ISTranspositionRate, self.ISElementLength, self.homeoticISTranspositionRate,
                                       self.homeoticISElementLength,
                                       self.DcISTranspositionRate, self.DcISElementLength)

            chromosome.RISTransposition(self.RISTranspositionRate, self.RISElementLength,
                                        self.homeoticRISTranspositionRate, self.homeoticRISElementLength)

            chromosome.geneTransposition(self.geneTranspositionRate)

            # 选择于i不同的染色体的，otherIndex
            otherIndex = i
            while otherIndex == i:
                otherIndex = random.randint(0, len(self.population) - 1)
            chromosome.onePointRecombination(self.onePointRecombinationRate, self.population[otherIndex])

            otherIndex = i
            while otherIndex == i:
                otherIndex = random.randint(0, len(self.population) - 1)
            chromosome.twoPointRecombination(self.twoPointRecombinationRate, self.population[otherIndex])

            otherIndex = i
            while otherIndex == i:
                otherIndex = random.randint(0, len(self.population) - 1)
            chromosome.geneRecombination(self.geneRecombinationRate, self.population[otherIndex])

    def addNeutralGene(self, inputsOutputs, neutralValue):
        chromosome = self.population[0]
        neutralGene = chromosome.createNeutralGene(inputsOutputs, neutralValue, numRandom=self.numRandom,
                                                   RandomRangeStart=self.RandomRangeStart,
                                                   RandomRangeEnd=self.RandomRangeEnd,
                                                   RandomPRECISION=self.RandomPRECISION)
        for chromosome in self.population:
            chromosome.updateNeutralGene(neutralGene=neutralGene)

    def printChromosomes(self, generation):
        '''
        打印染色体
        :param generation:   迭代的代数
        :return: None
        '''
        print("generation: ", generation)
        for Chromosome in self.population:
            Chromosome.printChromosome()

    @staticmethod
    def roulette_wheel_selection(weights):
        total = sum(weights)
        cum_weights = [sum(weights[:i + 1]) / total for i in range(len(weights))]
        rand_val = random.random()
        for i, weight in enumerate(cum_weights):
            if rand_val <= weight:
                return i

    def init_population(self, populationSize, numGenes, headLength, genome):
        function_weights = [func.operand_count for func in genome.functions]  # 假设这是函数的权重列表

        self.population = []
        for _ in range(populationSize):
            chromosome = Chromosome()  # 假设的染色体初始化
            # 对染色体头部进行轮盘赌选择
            for _ in range(headLength):
                selected_func_index = Environment.roulette_wheel_selection(function_weights)
                chromosome.add_gene_head(selected_func_index)  # 假设的方法来添加基因头部
            # 对染色体的其他部分进行简单初始化
            chromosome.initialize_rest()  # 假设的方法来初始化染色体的其余部分
            self.population.append(chromosome)
