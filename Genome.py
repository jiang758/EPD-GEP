
import random
import math

class Genome():
    '''对函数符运算进行定义、终结符定义'''
    random.seed(10)
    def __init__(self):
        '''定义函数、终结符
        functions   保存函数 如：*，/，+，—
        terminals   保存终结符 如：a，b
        genicTerminals 保存连接基因的终结符 如：0，1，2
        '''
        self.functions = {"max_arity": 0}
        self.terminals = []
        self.genicTerminals = range(1)
        self.PRECISION=0

    #定义 add ,sub ,mul ,div 函数的运算
    #进行保护，当除数为小于1e-6时，将除数设置为1

    def add(inputs):
        result=inputs.pop(-2)+ inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'

    def sub(inputs):
        result=inputs.pop(-2)- inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'

    def mul(inputs):
        result=inputs.pop(-2)* inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'
    def div(inputs):
        try:
            result=inputs.pop(-2)/ inputs.pop(-1)
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'
    def pow(inputs):
        try:
            result= inputs.pop(-1)**2
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'
    def reverse(inputs):
        try:
            result= 1/inputs.pop(-1)
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'


    def ln(inputs):
        try:
            result = math.log(inputs.pop(-1))
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'



    # 三角函数
    def sin(inputs):
        result=math.sin(inputs.pop(-1))
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'
    def cos(inputs):
        result=math.cos(inputs.pop(-1))
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'

    def tan(inputs):
        try:
            result=math.tan(inputs.pop(-1))
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'
    def exp(inputs):
        try:
            result=math.exp(inputs.pop(-1))
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'

    #开根号、
    def sqrt(inputs):
        try:
            result = math.sqrt(inputs.pop(-1))
            if math.isfinite(result):
                inputs.append(result)
            else:
                return 'error'
        except:
            return 'error'

    # #常数
    # def t(inputs):
    #     try:
    #         result = random.random()
    #         if math.isfinite(result):
    #             inputs.append(result)
    #         else:
    #             return 'error'
    #     except:
    #         return 'error'
    # #常数
    # def s(inputs):
    #     try:
    #         result = 2.5
    #         if math.isfinite(result):
    #             inputs.append(result)
    #         else:
    #             return 'error'
    #     except:
    #         return 'error'
    #布尔运算
    def NOT(inputs):
        result= not inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'
    def AND(inputs):
        result= inputs.pop(-2) and inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'
    def OR(inputs):
        result= inputs.pop(-2) or inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'
    def XOR(inputs):
        result= inputs.pop(-2) is not inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'
    def XNOR(inputs):
        result= inputs.pop(-2) is inputs.pop(-1)
        if math.isfinite(result):
            inputs.append(result)
        else:
            return 'error'

    #运算符特征【'+':2,'-':2,'a':0,'b':0】
    def symbols(self):
        symbols = self.functions.copy()
        symbols.pop("max_arity", None)
        for symbol in symbols.keys():
            symbols[symbol] = symbols[symbol][1]
        symbols.update({terminal: 0 for terminal in self.terminals + list(self.genicTerminals)})
        return symbols

    '''传统基因函数集合、连接基因函数集合：
    对于：{"+": (add, 2),"max_arity": 2}
    +是运算符，add保存是对应的运算，2是需要参数个数
    max_arity 是保存最大需要的参数个数，2是最大需要参数个数'''

    # arithmetic_set = {"+": (add, 2), "-": (sub, 2), "*": (mul, 2), "/": (div, 2),"ln":(ln,1),
    #                   "Q": (sqrt,1),"exp": (exp,1),"sin": (sin,1),"cos": (cos,1),"tan": (tan,1),
    #     "max_arity": 2}
    arithmetic_set = {"+": (add, 2), "-": (sub, 2), "*": (mul, 2), "/": (div, 2),"ln":(ln,1),
                      "Q": (sqrt,1),"exp": (exp,1),"P":(pow,1),"R":(reverse,1),"cos":(cos,1),"sin":(sin,1),"tan":(tan,1),
                      "max_arity": 2
                      # ,"t":(t,0),"s":(s,0)
                      }
    linker_set = {"+": (add, 2), "-": (sub, 2), "*": (mul, 2), "/": (div, 2),"ln":(ln,1),
                  "Q": (sqrt,1),"exp": (exp,1),"P":(pow,1),"R":(reverse,1),"cos":(cos,1),"sin":(sin,1),"tan":(tan,1),
                  "max_arity": 2
                  # ,"t":(t,0),"s":(s,0)
                  }
    # linker_set = {"+": (add, 2),"max_arity": 2}
    #定义 add ,sub ,mul ,div 化简函数的运算
    def stirng_add(inputs):
        inputs.append("(" + str(inputs.pop(-2)) + "+" + str(inputs.pop(-1)) + ")")
    def string_sub(inputs):
        inputs.append("(" + str(inputs.pop(-2)) + "-" + str(inputs.pop(-1)) + ")")
    def string_mul(inputs):
        inputs.append("(" + str(inputs.pop(-2)) + "*" + str(inputs.pop(-1)) + ")")
    def string_div(inputs):
        inputs.append("(" + str(inputs.pop(-2)) + "/" + str(inputs.pop(-1)) + ")")

    def string_ln(inputs):
        inputs.append("(" + "log"+"("+str(inputs.pop(-1))+")"+")")

    def string_sin(inputs):
        inputs.append("(sin("+str(inputs.pop(-1))+"))")

    def string_cos(inputs):
        inputs.append("(cos("+ str(inputs.pop(-1))+"))")

    def string_tan(inputs):
        inputs.append("(tan("+ str(inputs.pop(-1))+"))")

    def string_exp(inputs):
        inputs.append("(exp("+ str(inputs.pop(-1))+"))")

    def string_sqrt(inputs):
        inputs.append("(" + str(inputs.pop(-1)) + "**(1/2))")

    def string_pow(inputs):
        inputs.append("(" + str(inputs.pop(-1)) + "**2)")

    def string_reverse(inputs):
        inputs.append("(" + "1" + "/" + str(inputs.pop(-1)) + ")")

    # def string_t(inputs):
    #     inputs.append(random.random())
    #
    # def string_s(inputs):
    #     inputs.append(2.5)

    simplify_all_function_set = {"+": stirng_add, "-": string_sub,
                                 "*": string_mul, "/": string_div,
                                 "cos":string_cos,"sin":string_sin,
                                 "ln":string_ln,"Q":string_sqrt,
                                 "tan":string_tan,"exp":string_exp,
                                 "P":string_pow,"R":string_reverse
                                 # "t":string_t,"s":string_s
                                 }
