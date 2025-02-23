import random
import time
import numpy as np
import math

import Define
from Genome import Genome
from Environment import Environment
from read_data import read_data

# Define.SEED=5
Define.SEED = 8


def evalFunction(C_value, T_value, maxFitness_flag=False, M=100, DataCount=10, absoluteError_flagz=True,absoluteError_flag=True):
    '''
    评估函数,返回fitness值
    :param C_value:                 表示染色体,样本输出值
    :param T_value:                 表示染色体,样本真实值
    :param maxFitness_flag:         True 返回为最大适应度值
    :param M:                       选择范围 M 常数
    :param DataCount:
    :param absoluteError_flagz:
    :param absoluteError_flag:
    :return:
    '''
    '''
    :param C_value:         表示染色体,样本输出值
    :param T_value:         表示染色体,样本真实值
    :param maxFitness_flag: 用于是否为最大适应度值
    :param M:               M 常数
    :param DataCount:       数据的个数
    :return:
    '''
    if maxFitness_flag:
        return M * DataCount
    sum_fitness = 0
    sumresult = 0
    for i in range(len(C_value)):
        if not math.isfinite(C_value[i][0]):
            sum_fitness += 0
            pass
        if absoluteError_flagz:
            if absoluteError_flag:
                result = M - abs(C_value[i][0] - T_value[i][0])
                if result < 0: result = 0
                sum_fitness = result + sum_fitness
            else:
                result = (M - abs((C_value[i][0] - T_value[i][0]) / (T_value[i][0]+0.1)) * 100)
                if result < 0: result = 0
                sum_fitness = result + sum_fitness
        else:
            if absoluteError_flag:
                result=(C_value[i][0] - T_value[i][0])*(C_value[i][0] - T_value[i][0])
                if result < 0: result = 0
                sumresult=result+sumresult
                sumresult = pow(sumresult/DataCount,0.5)
                sum_fitness=1000*1/(1+sumresult)
            else:
                result=abs(C_value[i][0] - T_value[i][0])
                sumresult = result + sumresult
                sumresult=sumresult/DataCount
                sum_fitness = 1000 * 1 / (1 + sumresult)
    fitness = sum_fitness
    return fitness


if __name__ == '__main__':
    startTime = time.perf_counter()
    # ————————————————参数设置 Begin————————————————#
    generationCount = 4000  #迭代次数                                                                                                                                                                                                      000  # 迭代的代数3000
    populationSize = 20   # 一个群染色体的个数1600
    # DataCount = 10              #样本数据的个数（自动生成）
    headLength =6    # 基因的头部长度10
    numGenes =3  # 一个染色体中基因的个数5
    M = 100  # 选择访问M常数
    absoluteError_flagz = True #True不带根，False带根号（FT均方根差，FF绝对均差）
    absoluteError_flag = False  # True，True 是绝对误差，True,False是相对误差

    # ----连接基因 start----#
    homeotic_flag = False  # 是否有连接基因,False:用加法进行连接,True
    numHomeotics = 1  # 一个染色体中连接基因的个数
    homeoticHeadLength = 5  # 连接基因的头部长度
    homeoticMutationRate = 0.044  # 连接基因的突变率
    homeoticISTranspositionRate = 0.1  # 连接基因IS转座率
    homeoticISElementLength = [1, 2, 3]  # 连接基因IS元素长度
    homeoticRISTranspositionRate = 0.1  # 连接基因RIS转座率
    homeoticRISElementLength = [1, 2, 3]  # 连接基因RIS转座率

    # ----DC域参数 start----#
    DC_flag = False  # 是否有DC域
    numRandom = 10  # 随机数的个数
    RandomRangeStart = -1  # 随机数区间起点
    RandomRangeEnd = 1  # 随机数区间终点
    RandomPRECISION = 3  # 随机数精确度(保留几位小数)
    randomSetRate = 0.044  # 随机数突变率
    DcMutationRate = 0.044  # DC域的突变率
    DcISTranspositionRate = 0.1  # DC域特殊的IS转座率
    DcISElementLength = [1, 2, 3]  # DC域特殊的IS元素长度

    # ----中性基因start----#
    neutral_flag = False
    TriggerFitnessCount = 1000
    neutralValue = 0

    # ----变异概率参数 start----#
    mutationRate = 0.2  # 突变率0.044
    ISTranspositionRate = 0.2  # IS转座率0.1
    ISElementLength = [1, 2, 3]  # IS元素长度
    RISTranspositionRate =0.2  # RIS转座率0.1
    RISElementLength = [1, 2, 3]  # RIS元素长度
    geneTranspositionRate = 0.1  # 基因转座率
    onePointRecombinationRate = 0.3  # 单点重组率0.3
    twoPointRecombinationRate = 0.3  # 两点重组率0.3
    geneRecombinationRate = 0.2  # 基因重组率0.1

    # ----函数、终结符定义 start----#
    roulette_flag = False  # True 赌盘，False 锦标赛
    ###############增加读文件#################
    DataCount, num_terminators, inputsOutputs = read_data(".\\1.txt")  # 根据txt导入DataCount数据个数，参数个数
    print(DataCount, num_terminators)
    all_terminators = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'n', 'm', 'o', 'p', 'r']
    terminals = all_terminators[0:num_terminators]

    # 有DC域加入 “？” 终结符
    if DC_flag:
        terminals.append('?')
    genome = Genome()
    genome.functions.update(Genome.arithmetic_set)  # 基因的函数集合
    genome.terminals = terminals  # 基因的终结符
    link_genome = Genome()
    link_genome.functions.update(Genome.linker_set)  # 连接基因的函数集合
    # ————————————————参数设置 End————————————————#

    # ————————————————种群初始化 Begin————————————————#
    environment = Environment()
    environment.setRates(mutationRate=mutationRate, ISTranspositionRate=ISTranspositionRate,
                         ISElementLength=ISElementLength, RISTranspositionRate=RISTranspositionRate,
                         RISElementLength=RISElementLength,
                         geneTranspositionRate=geneTranspositionRate,
                         onePointRecombinationRate=onePointRecombinationRate,
                         twoPointRecombinationRate=twoPointRecombinationRate,
                         geneRecombinationRate=geneRecombinationRate,
                         homeoticMutationRate=homeoticMutationRate,
                         homeoticISTranspositionRate=homeoticISTranspositionRate,
                         homeoticISElementLength=homeoticISElementLength,
                         homeoticRISTranspositionRate=homeoticRISTranspositionRate,
                         homeoticRISElementLength=homeoticRISElementLength,
                         randomSetRate=randomSetRate, DcMutationRate=DcMutationRate,
                         DcISTranspositionRate=DcISTranspositionRate,
                         DcISElementLength=DcISElementLength, numRandom=numRandom, RandomRangeStart=RandomRangeStart,
                         RandomRangeEnd=RandomRangeEnd,
                         RandomPRECISION=RandomPRECISION, TriggerFitnessCount=TriggerFitnessCount,
                         neutralValue=neutralValue)

    environment.init(populationSize=populationSize, numGenes=numGenes, numHomeotics=numHomeotics,
                     headLength=headLength, homeoticHeadLength=homeoticHeadLength, genome=genome,
                     link_genome=link_genome, homeotic_flag=homeotic_flag, DC_flag=DC_flag, neutral_flag=neutral_flag
                     )

    result = environment.run(inputsOutputs, evalFunction, generationCount=generationCount, M=M, DataCount=DataCount,
                             absoluteError_flagz=absoluteError_flagz,absoluteError_flag=absoluteError_flag, roulette_flag=roulette_flag)
    # ————————————————Run end————————————————#
    endTime = time.perf_counter()
    print("Time:", endTime - startTime)
