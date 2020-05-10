# ##################################################################################################
#  Disclaimer                                                                                      #
#  This file is a python3 translation of AutoDockTools (v.1.5.7)                                   #
#  Modifications made by Valdes-Tresanco MS (https://github.com/Valdes-Tresanco-MS)                #
#  Tested by Valdes-Tresanco-MS and Valdes-Tresanco ME                                             #
#  There is no guarantee that it works like the original distribution,                             #
#  but feel free to tell us if you get any difference to correct the code.                         #
#                                                                                                  #
#  Please use this cite the original reference.                                                    #
#  If you think my work helps you, just keep this note intact on your program.                     #
#                                                                                                  #
#  Modification date: 10/5/20 18:04                                                                #
#                                                                                                  #
# ##################################################################################################

import numpy


class ConstrainedParameterSet:

    def __init__(self):
        # conDict records constraints between parameters
        # key is parm1Name,  name of parameter1 ,
        # value is list of triples: (parm2Name, func, args)
        # changes in parm1 cause parm2 to be updated by func
        # eg: to enforce
        # self.center = 2*self.offset
        # so that changes in self.offset force changes in self.center
        # self.conDict['offset'] = [('center', numpy.multiply, ('offset, 2.0'))]
        # to enforce the reciprocal constraint
        # so that changes in self.center force changes in self.offset
        # self.conDict['center'] = [('offset', numpy.multiply, ('center,0.5'))]
        # suitable functions are numpy.divide and numpy.multiply
        # DO SOMETHING SO THIS DOESN'T GET into an endless loop
        #
        self.conDict = {}
        # rangeDict provides methods of checking validity
        # of a given value for key
        # possible types methods include:
        # list of discrete values, tuple defining valid range, a type
        # or a function which returns 1 if value is valid or 0 if not
        self.rangeDict = {}
        # values can be integers, floats, strings....???
        self.typesList = [type(1), type(1.0), type('a')]

    def tie(self, parm1Name, parm2Name, func, args):
        # eg:
        # changes in self.center force  changes in self.offset
        # self.tie('center','offset',numpy.multiply,'center,2.0')
        if parm1Name not in self.conDict:
            self.conDict[parm1Name] = []
        self.conDict[parm1Name].append((parm2Name, func, args))
        # is this at all clear?
        # cD = self.conDict
        # cD[parm1Name] = cD.get(parm1Name, []).append((parm2Name, func, args))

    def updateConstraints(self, parmName):
        # called when parm changes to update other linked parms
        conList = self.conDict.get(parmName, None)
        if not conList:
            # do nothing + return
            print('no constraints on ', parmName)
            return
        # conList has tuples (func, args)
        # eg: sample value in conList
        # (parm2Name, numpy.multiply, (parmName, 0.5))
        for parm2Name, func, args in conList:
            # FIX THIS:
            # to update self.parm2Name:
            #   need to get (self.center, 0.5) from args='center, 0.5'
            setattr(self, parm2Name, func(*eval('self.' + args)))

    def untie(self, parm1Name, parm2Name, func, args):
        # eg:
        # g.untie('center','offset',numpy.multiply,'center,2.0')
        conList = self.conDict.get(parm1Name, None)
        if not conList:
            print('no constraints on ', parm1Name)
            return "ERROR"
        if (parm2Name, func, args) not in conList:
            print('(%s,%s,%s) not in %s constraints' % (parm2Name, func, args, parm1Name))
            return "ERROR"
        self.conDict[parm1Name].remove((parm2Name, func, args))

    def setRange(self, parm, range):
        # range can be list, interval tuple, type or validation func
        # FIX THIS: do some range validation
        self.rangeDict[parm] = range

    def validateValue(self, parm, value):
        rangeD = self.rangeDict
        if parm not in rangeD:
            # nothing specified for this parm
            return value
        range = rangeD[parm]
        if type(range) == list:
            if value in range:
                return value
            else:
                return "ERROR: value not in range list"
        elif type(range) == tuple:
            if value >= range[0] and value <= range[1]:
                return value
            else:
                return "ERROR: value not in range interval"
        elif range in self.typesList:
            if type(value) == range:
                return value
            else:
                return "ERROR: value not specified type"
        else:
            # only thing left is validation function
            ok = list(range(*value))
            if ok:
                return value
            else:
                return "ERROR: value failed validation func"

    def update(self, parm, value):
        check = self.validateValue(parm, value)
        if check != value:
            print('failed validation:\n', check)
            return "ERROR"
        self.updateConstraints(parm)

    def fix(self, parmName):
        # ????????????????????????????
        # this method makes self.parmName constant
        # by removing any constraints which force it to change
        for k, v in list(self.conDict.items()):
            for triple in v:
                if triple[0] == parmName:
                    self.conDict[k].remove(triple)


class GeneralRegularGridDescriptor(ConstrainedParameterSet):
    keywords = ['center',
                'offset',
                'length',
                'nbGridPoints',
                'gridSpacing'
                ]

    def __init__(self, **kw):

        ConstrainedParameterSet.__init__(self)

        for k in self.keywords:
            setattr(self, k, kw.get(k, numpy.array((0., 0., 0.))))

        consDict = kw.get('consDict', None)
        if consDict:
            for k, v in list(consDict.items()):
                self.tie(k, v[0], v[1], v[2])

        rangeDict = kw.get('rangeDict', None)
        if rangeDict:
            for k, v in list(rangeDict.items()):
                self.setRange(k, v)

        fixed = kw.get('fixed', None)
        if fixed:
            # fixed should be a list of parameters to fix
            # same problem of type(parm)...a string???
            for p in fixed:
                self.fix(p)
