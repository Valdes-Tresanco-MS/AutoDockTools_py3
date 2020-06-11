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
#  Modification date: 10/05/20, 6:51 p. m.                                                         #
#                                                                                                  #
# ##################################################################################################

#############################################################################
# 
# Histogram class: 
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr> 
# last revision: 1999-7-6
#
# HistogramRI class written by Ruth Huey
# HistogramRI adds a reverse Index to Histogram class 
#
#############################################################################


import numpy


class Histogram:
    """Histogram in one variable

    Constructor: Histogram(|data|, |bins|, |range|=None)

    Arguments:

    |data| -- a sequence of data points

    |bins| -- the number of bins into which the data is to be sorted

    |range| -- a tuple of two values, specifying the lower and
               the upper end of the interval spanned by the bins.
               Any data point outside this interval will be ignored.
               If no range is given, the smallest and largest
               data values are used to define the interval.

    The number of points in a bin can be obtained by indexing the
    histogram with the bin number. Application of len() yields the
    number of bins. A histogram thus behaves like a sequence of
    numbers.
    """

    def __init__(self, data, nbins, range=None):
        if range is None:
            self.min = numpy.minimum.reduce(data)
            self.max = numpy.maximum.reduce(data)
        else:
            self.min, self.max = range
        self.min = self.min + 0.
        self.max = self.max + 0.
        self.bin_width = (self.max - self.min) / nbins
        if self.bin_width == 0:
            print('range is 0 so set bin_width to 1.')
            self.bin_width = 1.
        self.array = numpy.zeros((nbins, 2), numpy.float)
        self.array[:, 0] = self.min + self.bin_width * (numpy.arange(nbins) + 0.5)
        self.addData(data)

    def __len__(self):
        return self.array.shape[0]

    def __getitem__(self, index):
        return self.array[index]

    def __getslice__(self, first, last):
        return self.array[first:last]

    def addData(self, data):
        """Add the values in |data| (a sequence of numbers) to the
        originally supplied data. Note that this does not affect the
        default range of the histogram, which is fixed when the
        histogram is created.
        """
        n = (len(data) + 999) / 1000
        for i in range(n):
            self._addData(data[1000 * i:1000 * (i + 1)])

    def _addData(self, data):
        data = numpy.array(data, numpy.float)
        data = numpy.repeat(data, numpy.logical_and(numpy.less_equal(data, self.max),
                                                    numpy.greater_equal(data, self.min)))
        data = numpy.floor((data - self.min) / self.bin_width).astype(numpy.int)
        self.rIdata = data
        nbins = self.array.shape[0]
        histo = numpy.add.reduce(numpy.equal(numpy.arange(nbins)[:, None], data), -1)
        histo[-1] = histo[-1] + numpy.add.reduce(numpy.equal(nbins, data))
        self.array[:, 1] = self.array[:, 1] + histo

    def normalize(self, norm=1.):
        "Scales all counts by the same factor such that their sum is |norm|."
        self.array[:, 1] = norm * self.array[:, 1] / numpy.add.reduce(self.array[:, 1])


class HistogramRI(Histogram):
    """ This class is based on Histogram class developed by K.Hinsen. It adds
    a method, 'createReverseIndex', which builds a list of data points in each
    bin"""

    def __init__(self, data, nbins, range=None):
        Histogram.__init__(self, data, nbins, range)
        self.data = data

    def createReverseIndex(self):
        # first build the data in the correct form:
        nbins = self.array.shape[0]
        data = numpy.array(self.data, numpy.float)
        data = numpy.repeat(data, numpy.logical_and(numpy.less_equal(data, self.max),
                                                    numpy.greater_equal(data, self.min)))
        data = numpy.floor((data - self.min) / self.bin_width).astype(numpy.int)
        self.cRI_data = data
        # check for pt outside of nbins
        # sometimes max point falls outside of right-most bin
        # first stuff works with python2.0 second with both
        try:
            self.cRI_data = numpy.putmask(data, numpy.greater_equal(data, nbins), nbins - 1)
        except:
            for i in range(data.shape[0]):
                if data[i] >= nbins: data[i] = nbins - 1
        self.reverseIndex = []
        for i in range(nbins):
            newentry = numpy.nonzero(numpy.equal(data, i))[0].tolist()
            # print i,'th entry =', newentry
            self.reverseIndex.append(newentry)


if __name__ == '__main__':
    data = numpy.arange(500)
    # data = numpy.arange(5000)
    # h = Histogram(data, 10)
    hRI = HistogramRI(data, 10)
