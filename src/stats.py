#!/usr/bin/python 
# Begin stats.py
"""Module for statistics.

Written by Jesse Bloom."""
#
import re, math, os, string, random 
try:
    import numarray.fft
except ImportError:
    pass
#----------------------------------------------------------------------------
class StatsError(Exception):
    """Statistics error."""
#--------------------------------------------------------------------------
def DistWidth(xlist, frac, accuracy = 1.0e-6):
    """Calculates the width of a distribution.

    Call is: 'w = DistWidth(xlist, frac, [accuracy = 1.0e-6]'
    'xlist' is a list of numeric data.  Any values of 'None' are discarded.
    'frac' is the fraction of the distribution that we want to lie within
        certain bounds.  The returned value 'w' is the width we must go
        on either side of the mean of the distribution in order to have
        'frac' of the data points lie between mean - w  and mean + w.
    'accuracy' is the resolution within which 'w' is determined.  Default
        is 1.0e-6."""
    if not 0.0 <= frac <= 1.0:
        raise StatsError("Invalid 'frac' of %r." % frac)
    xdata = []
    for x in xlist:
        if x != None:
            xdata.append(x)
    xdata.sort()
    step = math.sqrt(Variance(xdata)) # initial step size
    n = len(xdata)
    mean = Mean(xdata)
    nwithin = frac * n # number of points to make 'frac'
    # make an initial guess assuming a symmetric distribution
    w = xdata[int(n / 2)] - xdata[int(n * (1 - frac) / 2)]
    #-----------------------------------------------------
    def NumWithin(min, max, xdata):
        """Returns number of entries in 'xdata' between 'min' and 'max'.
        
        'xdata' is sorted."""
        for i in range(len(xdata)):
            if xdata[i] >= min:
                break
        ilow = i
        for i in range(ilow, len(xdata)):
            if xdata[i] >= max:
                break
        return i - ilow
    #-----------------------------------------------------
    if NumWithin(mean - w, mean + w, xdata) > nwithin:
        toomany = True
    else:
        toomany = False
    while True:
        if toomany:
            w -= step 
            if NumWithin(mean - w, mean + w, xdata) <= nwithin:
                if step < accuracy:
                    break
                else:
                    toomany = not toomany
                    step /= 2.0
        else:
            w += step 
            if NumWithin(mean - w, mean + w, xdata) >= nwithin:
                if step < accuracy:
                    break
                else:
                    toomany = not toomany
                    step /= 2.0
    return w
#--------------------------------------------------------------------------
def List_of_Tuples_to_Tuple_of_Lists(xlist):
    """Converts a list of tuples to a tuple of lists.

    Takes as input a list 'xlist'.  Each element of 'list' must be a tuple
        of the same length.  Returns a tuple of lists, with the same number
        of lists as elements in the tuples in 'list', and with the first list
        containing all of the first elements in the tuples in 'list', et cetera."""
    if not isinstance(xlist, list):
        raise StatsError("Error, the argument is not a list.")
    n = None
    for tup in xlist:
        if not isinstance(tup, tuple):
            raise StatsError("All list elements are not tuples.")
        if n == None:
            n = len(tup)
        elif len(tup) != n:
            raise StatsError("All tuples do not have the same length.")
    list_of_lists = []
    for i in range(n):
        ilist = []
        for tup in xlist:
            ilist.append(tup[i])
        assert len(ilist) == len(xlist) 
        list_of_lists.append(ilist)
    assert len(list_of_lists) == n
    return tuple(list_of_lists)
#--------------------------------------------------------------------------
def Write_Columns(filename, *data_lists):
    """Writes data lists as columns to a file.

    'filename' is the name of the file to which the data should be written.
    An arbitrary number of lists can be specified.  Writes each list as
        a column, in the order they are given, to 'filename'.
    A list can be specified alone as a list, or in a 2-tuple
        as '(title, list)' where title is a title printed in a title line 
        (beginning with a '#') for that list.  If one list has a title,
        all must have them.
    Finally, any lists that are passed as the 2-tuple
        '("LIST_OF_LISTS", l_of_l)' interprets l_of_l to be lists
        of lists or 2-tuples, each of which is handled as if it were
        a separate argument."""
    if len(data_lists) == 0:
        raise StatsError("No data lists were specified.")
    i = 0
    data_lists = list(data_lists)
    while i < len(data_lists):
        if isinstance(data_lists[i], tuple) and data_lists[i][0] == 'LIST_OF_LISTS':
            l_of_l = data_lists[i][1]
            data_lists.pop(i)
            data_lists += l_of_l
        else:
            i += 1
    n = None
    if isinstance(data_lists[0], list):
        have_lists = True
    elif isinstance(data_lists[0], tuple):
        have_lists = False
    else:
        raise StatsError("First 'data_lists' value is neither list nor tuple.")
    for xentry in data_lists:
        if (not have_lists) and (isinstance(xentry, tuple) and len(xentry) == 2):
            if not isinstance(xentry[0], str):
                raise StatsError("First entry of tuple is not a string.")
            xlist = xentry[1]
        elif have_lists:
            xlist = xentry
        else:
            raise StatsError("Not a valid data lists entry.")
        if not isinstance(xlist, list):
            raise StatsError("Specified list is not a list.")
        if n == None:
            n = len(xlist)
        elif n < len(xlist):
            n = len(xlist)
    file = open(filename, 'w')
    if not have_lists:
        file.write('#')
        for itup in range(len(data_lists) - 1):
            file.write(data_lists[itup][0] + '\t')
        file.write(data_lists[len(data_lists) - 1][0] + '\n')
    for i in range(n):
        for itup in range(len(data_lists) - 1):
            if have_lists:
                ilist = data_lists[itup]
            else:
                ilist = data_lists[itup][1]
            try:
                file.write(str(ilist[i]) + '\t')
            except IndexError:
                file.write('NA\t')
        if have_lists:
            ilist = data_lists[len(data_lists) - 1]
        else:
            ilist = data_lists[len(data_lists) - 1][1]
        try:
            file.write(str(ilist[i]) + '\n')
        except IndexError:
            file.write('NA\n')
    file.close()
#---------------------------------------------------------------------------
def Make_Bins(data, binsize = None, numbins = None, numbins_10_90 = None, min_max = None, gnuplot_file = None):
    """Makes histogram bins for a data set.

    On call, 'data' is a list of numbers specifying the data.
    Any values of 'None' in the data list are discarded.
    The data is binned either so all bins are of a specified size, 
        or that there are a specified number of bins.  If all bins are
        all the same size, 'binsize' is this size. If there are
        to be a specified number of bins, 'numbins' is this value.  
        Only one of 'numbins' or 'binsize' should be specified.
    For both the 'binsize' and 'numbins' options, the user can specify 
        the minimum of the smallest bins and the maximum of 
        the largest bin in the 'min_max' (example, 'min_max = (0, 4)'
        means the smallest bin has minimum 0, and the largest bin
        bin has maximum 4.  If 'binsize' is chosen, the maximum
        value is adjusted upward to make it consistent with the bin sizes.
    'numbins_10_90' is like 'numbins', but the bin size is chosen so that
        there are 'numbins_10_90' bins between the 10th and 90th 
        percentiles of the data.
    Upper bin boundaries are less than, lower boundaries are >=.
    The function returns a 4-tuple: (bincounts, bincenter, binmin, binmax)
        'bincounts' is a list with 'bincounts[i]' being the number of
        of elements in bin i (0 <= i < numbins). 'bincenter[i]' is
        the center of bin 'i', 'binmin[i]' is the lower limit
        of bin i, and 'binmax[i]' is the upper limit of bin 'i'.
    If there is a problem, the function returns 'None'.
    If 'gnuplot_file' is set to a value other than 'None', 
        it should be a string representing a filename.  This
        This file is then created, with the data in gnuplot format."""
    assert isinstance(data, list)
    xdata = [] # copy of data we use here
    for x in data:
        if isinstance(x, (int, float)):
            xdata.append(x)
        elif x == None:
            pass
        else:
            raise StatsError("Invalid data entry of %r." % x)
    xdata.sort()  # sort the list from smallest to biggest
    small = xdata[0] # smallest value
    big = xdata[len(xdata) - 1] # biggest value
    if min_max != None: # set the minimum and maximum
        assert isinstance(min_max, tuple)
        assert len(min_max) == 2
        assert min_max[0] <= min_max[1]
        if small < min_max[0]:
            raise StatsError("Smallest data value of %r is smaller than the specified minimum of %r." % (small, min_max[0]))
        else:
            small = min_max[0]
        if big > min_max[1]:
            raise StatsError("Largest data value of %r is larger than the specified minimum of %r." % (large, min_max[1]))
        else:
            big = min_max[1]
    big += (big - small) / 1000000.0 # add a small number to big to account for rounding problems
    assert small <= big
    if binsize: # bins of specified size
        if numbins != None:
            raise StatsError("Both binsize and numbins are set.")
        assert isinstance(binsize, (int, float))
        assert binsize > 0
    elif numbins:  # specified number of bins, compute the binsize
        assert isinstance(numbins, int)
        assert numbins > 0
        binsize = (big - small) / float(numbins)
    elif numbins_10_90:
        x10 = xdata[int(len(xdata) / 10.0)]
        x90 = xdata[int(9.0 * len(xdata) / 10.0)]
        binsize = (x90 - x10) / float(numbins_10_90)
    else:
        raise StatsError("No bin size or number specification is set.")
    # Do the binning
    bincounts = []
    bincenter = []
    binmin = []
    binmax = []
    ibin = index = 0
    bincounts.append(0)
    binmin.append(small)
    binmax.append(binmin[ibin] + binsize)
    bincenter.append((binmin[ibin] + binmax[ibin]) / 2.0)
    xlength = len(xdata)
    while True:  # keep adding bins until the break statement
        while index < xlength and xdata[index] < binmax[ibin]: # add data counts
            bincounts[ibin] += 1
            index += 1
        if index == xlength: # we have binned all of the data
            break
        else: # next bin
            ibin += 1
            bincounts.append(0)
            binmin.append(binmin[ibin - 1] + binsize)
            binmax.append(binmin[ibin] + binsize)
            bincenter.append((binmin[ibin] + binmax[ibin]) / 2.0)
    if numbins and len(bincounts) != numbins:
        raise StatsError("Wanted %r bins but got %r." % (numbins, len(bincounts)))
    assert len(bincounts) == len(binmax) == len(binmin) == len(bincenter)
    assert binmax[len(binmax) - 1] > xdata[len(xdata) - 1]
    # write data to a gnuplot file
    if gnuplot_file:
        assert isinstance(gnuplot_file, str)
        file = open(gnuplot_file, 'w')
        file.write('#bin_min\tbin_max\tbin_center\tcounts\n')
        for i in range(len(bincounts)):
            bmin = binmin[i]
            bmax = binmax[i]
            cent = str((bmin + bmax) / 2.0)
            bmin = str(bmin)
            bmax = str(bmax)
            c = str(bincounts[i])
            file.write(bmin + '\t' + bmax + '\t' + cent + '\t' + c + '\n')
        file.close()
    return (bincounts, bincenter, binmin, binmax)  # return the final tuple
#------------------------------------------------------------------------
def StatsSummary(numlist):
    """Returns summary one-variable statistics for a list of numbers.
    
    Call is: '(mean, median, sd, n) = StatsSummary(numlist)'
    Any entries of 'None' or '-' are removed.
    If there are one or zero data points, just returns the number of
        data points rather than the 4-tuple.
    Returns a 4-tuple giving the mean, median, standard deviation, and
        number of data points."""
    n = len(numlist) - numlist.count(None) - numlist.count('-')
    if n <= 1:
        return n
    return (Mean(numlist), Median(numlist), StandardDeviation(numlist), n)
#-------------------------------------------------------------------------
def Median(numlist):
    """Returns the median of a list of numbers.
    
    If any entries of the list are 'None' or '-', they are removed first."""
    assert isinstance(numlist, list)
    xlist = []  # make a copy of the list
    for x in numlist:
        if isinstance(x, (int, float)):
            xlist.append(x)
        elif x in [None, '-']:
            pass
        else:
            raise StatsError("Invalid value of %r in list." % x)
    if len(xlist) == 0:
        raise StatsError("Empty list.")
    xlist.sort()
    n = len(xlist)
    if n % 2 == 0: # even length list, average two middle entries
        med = xlist[n / 2] + xlist[n / 2 - 1]
        med = float(med) / 2.0
        return med
    else: # odd length list, get middle entry
        return xlist[n / 2]
#--------------------------------------------------------------------------
def Mean(numlist):
    """Returns the mean of a list of numbers.
    
    If any entries of the list are 'None' or '-', they are removed first."""
    mean = 0.0
    n = 0
    for x in numlist:
        if x in [None, '-']:
            continue
        if not isinstance(x, (int, float)):
            raise StatsError("Invalid entry of %r." % x)
        mean += x
        n += 1
    assert n == len(numlist) - numlist.count(None) - numlist.count('-')
    if n <= 0:
        raise StatsError("Empty list.")
    return mean / float(n)
#-------------------------------------------------------------------------
def Variance(numlist):
    """Returns the variance of a list of numbers.

    If any entries of the list are 'None' or '-', they are removed first."""
    sum = sum2 = 0.0
    n = 0
    for x in numlist:
        if x in [None, '-']:
            continue
        assert isinstance(x, int) or isinstance(x, float)
        sum += x
        sum2 += x * x
        n += 1
    try:
        mean2 = sum2 / float(n)
        mean = sum / float(n)
    except ZeroDivisionError:
        raise StatsError( "Empty list.")
    var = mean2 - mean * mean
    assert var >= -0.00000001
    return var 
#------------------------------------------------------------------------
def StandardDeviation(numlist):
    """Returns the standard deviation of a list of numbers.

    If any entries of the list are 'None' or '-', they are removed first."""
    return math.sqrt(Variance(numlist))
#--------------------------------------------------------------------------
def Read_Numeric_Columns(filename, columnlist):
    """Returns numbers in columns of the indicated file.

    'filename' is a file containing columns of numeric entries, with 
        the columns separated by whitespace.
    'columnlist' is a list of the columns for which we want to 
        get entries.  The columns range from 1 to n (NOT ZERO TO n).
    For all rows of the file that contain numeric entries in all of the 
        columns indicated in 'columnlist', creates lists of the entries in 
        the columns and returns them as a list: [columnX, columnY, ...] 
        where the columns in the list correspond to those given in 
        'columnlist'. Note that rows without numeric entries in all
        of the indicated columns are disregarded.
    Returns an empty list if nothing can be read."""
    columns = []
    if not os.path.isfile(filename):
        raise StatsError("Cannot find file %r." % filename)
    if not columnlist:
        raise StatsError("The column list is empty.")
    # Make sure all elements of 'columnlist' are valid integers
    for i in range(len(columnlist)):
        if (not isinstance(columnlist[i], int)) or (columnlist[i] < 1):
            raise StatsError("Entry %r is not a valid column index." % columnlist[i])
    # Get column entries
    file = open(filename, 'r')
    firstitems = True 
    digitmatch = re.compile(r"^\-{0,1}[\d\.]+$")
    while True: # read lines from columns
        line = file.readline()
        if not line:
            break
        if line[0] == '#':
            continue
        items = string.split(line)
        row = []
        for column in columnlist:
            if column <= len(items):
                if digitmatch.search(items[column - 1]):
                    row.append(float(items[column - 1]))
        if len(row) == len(columnlist): # if we find an item for every column
            if firstitems:
                for column in columnlist:
                    columns.append([])
                firstitems = False 
            for column in range(len(row)):
                columns[column].append(row[column])
    file.close()
    return columns
#--------------------------------------------------------------------------
def Kendalls_Tau(xlist, ylist):
    """Calculates Kendall's tau non-parametric correlation.

    The input data is given in the two lists 'xdata' and 'ydata' which should be 
        of the same length.  If entry i of either list is 'None', this entry is
        disregarded in both lists.
    Returns Kendall's partial tau, the one-tailed P-value, and the number of
        data points as a tuple: (tau, P, N).
    Includes a correction for ties.
    Based on Gibbons, JD, "Nonparametric measures of association", 
        Sage University Papers, pg 15 (1983)."""
    if len(xlist) != len(ylist):
        raise StatsError("Data lists are of different length.")
    xdata = []
    ydata = []
    for i in range(len(xlist)):
        if xlist[i] != None and ylist[i] != None:
            xdata.append(xlist[i])
            ydata.append(ylist[i])
    assert len(xdata) == len(ydata)
    assert len(xdata) <= len(xlist) - xlist.count(None)
    assert len(ydata) <= len(ylist) - ylist.count(None)
    assert len(ydata) >= len(ylist) - xlist.count(None) - ylist.count(None)
    if len(xdata) == 0:
        raise StatsError("No valid data entries.")
    n = len(xdata)
    # compute the number of concordant and discordant pairs
    conc = disc = 0.0 # concordant and discordant pairs
    for i in range(n): # loop over all pairs
        xi = xdata[i]
        yi = ydata[i]
        for j in range(i + 1, n):
            xd = xi - xdata[j]
            yd = yi - ydata[j]
            prod = xd * yd
            if prod == 0.0: # this is a tie
                continue
            elif prod > 0.0:
                conc += 1
            else:
                disc += 1
    # compute the tie correction: sum(t * t - t)
    xcopy = []
    ycopy = []
    for i in range(n):
        xcopy.append(xdata[i])
        ycopy.append(ydata[i])
    xties = yties = 0.0
    while xcopy: 
        xi = xcopy[0]
        t = xcopy.count(xi)
        xties = xties + t * t - t
        while xcopy.count(xi) > 0:
            xcopy.remove(xi)
    while ycopy: 
        yi = ycopy[0]
        t = ycopy.count(yi)
        yties = yties + t * t - t
        while ycopy.count(yi) > 0:
            ycopy.remove(yi)
    # Compute tau
    n = float(n)
    denom = math.sqrt((n * n - n - xties) * (n * n - n - yties))
    try:
        tau = 2.0 * (conc - disc) / denom
    except ZeroDivisionError:
        raise StatsError("Too few entries: %r." % n)
    # Compute P-value
    z = 3.0 * tau * math.sqrt(n * (n - 1.0)) / math.sqrt(2.0 * (2.0 * n + 5.0))
    prob = Complementary_Error_Function(z / math.sqrt(2.0)) / 2.0
    return (tau, prob, int(n))
#--------------------------------------------------------------------------
def Kendalls_Partial_Tau(xdata, ydata, zdata):
    """Computes Kendall's partial tau of two variables controlling for a third.

    The correlation is between 'xdata' and 'ydata' controlling for 'zdata'.
        The data is given in lists that must be of the same length.
    Returns partial tau as a scalar number.
    Based on Gibbons JD, "Nonparametric measures of associations", 
        Sage University Papers, pg 49 (1983)."""
    if not len(xdata) == len(ydata) == len(zdata):
        raise StatsError("Data sets are of different lengths.")
    txy = Kendalls_Tau(xdata, ydata)[0]
    tyz = Kendalls_Tau(ydata, zdata)[0]
    txz = Kendalls_Tau(xdata, zdata)[0]
    partial_tau = (txy - txz * tyz) / math.sqrt((1 - txz * txz) * (1 - tyz * tyz))
    return partial_tau
#--------------------------------------------------------------------------
def LinearRegression(xdata, ydata):
    """Fits a line to a data set.

    Call is: '(a, b) = LinearRegression(xdata, ydata)'
    'xdata' and 'ydata' are lists of numbers, and should be of the
        same length.  If entry 'i' is 'None' for either list,
        this entry is disregarded for both lists.
    The function fits a line to the data, and returns the slope of
        that line as 'a' and the y-intercept as 'b'.  So the data
        is modeled by: y = a + b * x.
    Note that this function does no checking to see how good this fit
        is -- you should probably also check to make sure that your data
        is actually linear!."""
    if len(xdata) != len(ydata):
        raise StatsError( "Data sets are of different length.")
    xlist = []
    ylist = []
    for i in range(len(xdata)):
        if xdata[i] != None and ydata[i] != None:
            xlist.append(xdata[i])
            ylist.append(ydata[i])
    n = len(xlist)
    if n < 2:
        raise StatsError("We have less than two data points, so we can't fit a line.")
    sx = sy = sxx = sxy = 0.0
    for i in range(len(xdata)):
        x = xdata[i]
        y = ydata[i]
        s = n
        sx += x
        sy += y
        sxx += x * x
        sxy += x * y
    delta = s * sxx - sx * sx
    a = (sxx * sy - sx * sxy) / delta
    b = (s * sxy - sx * sy) / delta
    return (a, b)
#----------------------------------------------------------------------------
def PearsonCorrelation(xdata, ydata):
    """Computes the Pearson linear correlation between two data sets.
   
    Call is '(r, p, n) = PearsonCorrelation(xdata, ydata)'
    The input data is given in the two lists 'xdata' and 'ydata' which should
        be of the same length.  If entry i of either list is 'None', this
        entry is disregarded in both lists.
    Returns Pearson's correlation coefficient, the two-tailed P-value, 
        and the number of data points as a tuple '(r, p, n)'."""
    if len(xdata) != len(ydata):
        raise StatsError("Data sets are of different lengths.")
    xlist = []
    ylist = []
    xylist = []
    for i in range(len(xdata)):
        if xdata[i] != None and ydata[i] != None:
            xlist.append(xdata[i])
            ylist.append(ydata[i])
            xylist.append(xdata[i] * ydata[i])
    n = len(xlist)
    if n <= 1:
        raise StatsError("One or less data points: %r." % n)
    xmean = Mean(xlist)
    ymean = Mean(ylist)
    xymean = Mean(xylist)
    xsd = StandardDeviation(xlist)
    ysd = StandardDeviation(ylist)
    try:
        r = float(xymean - xmean * ymean) / (xsd * ysd)
    except ZeroDivisionError:
        raise StatsError("Standard deviation is zero.")
    if not (-1.0000000001 <= r <= 1.000000001):
        raise StatsError("Invalid correlation coefficient of %r." % r)
    z = math.fabs(r) * math.sqrt(n) / math.sqrt(2.0)
    p = Complementary_Error_Function(z)
    if not (0.0 <= p <= 1.0):
        raise StatsError("Invalid P-value of %r." % r)
    return (r, p, n)
#---------------------------------------------------------------------
def Complementary_Error_Function(z):
    """Calculates the error function of z.

    The complementary error function of z is defined as:
        erfc(z) = 2 / sqrt(pi) * integral(e^(t^2) dt) where the integral 
        is from z to infinity.
    Can be used to calculate cumulative normal probabilities: given a 
        distribution with mean m and standard deviation s,
        the probability of observing x > m  when x > 0 is: 
        P = 0.5 * erfc((x - m) / (s * sqrt(2)))
    Calculated according to Chebyshev fitting given by Numerical Recipes 
        in C, page 220-221."""
    x = math.fabs(z)
    t = 1.0 / (1.0 + 0.5 * x)
    ans = t * math.exp(-x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    if z > 0.0:
        return ans
    else:
        return 2.0 - ans
#-------------------------------------------------------------------
def Poisson(n, k):
    """Returns the Poisson probability of observing a number.

    'Poisson(n, k)' takes as input an integer n >= 0 and a 
        real number k >= 0.0.  Returns p, the probability of
        getting n counts when the average outcome is k, according
        to the Possion distribution.  Returns 'None' if there is
        an error."""
    p = math.exp(-k) * math.pow(k, n) / float(Factorial(n))
    assert 0.0 <= p <= 1.0, "Error, value of p is invalid probability: " + str(p)
    return p
#-------------------------------------------------------------------
def Factorial(n):
    """Returns the factorial of an integer."""
    x = 1
    for i in range(1, n + 1):
        x *= i
    return x
#-----------------------------------------------------------------
def NChooseK(n, k):
    """Computes nCrK.

    Call is: 'x = NChooseK(n, k)'
    'n' and 'k' are non-negative integers with 'n' >= 'k'.
    'x' is returned as nCrk, defined as: n! / ((n - k)! * k!)."""
    x = 1
    for i in range(max(k + 1, n - k + 1), n + 1):
        x *= i
    x /= Factorial(min(k, n - k))
    return x
#-------------------------------------------------------------------
def RandomGaussian(sigma):
    """Returns a number distributed according to a Gaussian distribution.

    Call is: 'x = RandomGaussian(sigma)'
    'x' is a number drawn from a Gaussian distribution with standard deviation
        'sigma'."""
    # generate two random numbers on unit circle
    r2 = 0.0
    while (r2 >= 1.0 or r2 == 0.0):
        x1 = 2.0 * random.random() - 1.0
        x2 = 2.0 * random.random() - 1.0
        r2 = x1**2 + x2**2
    fac = math.sqrt(-2.0 * math.log(r2) / r2)
    x1 *= fac * sigma
    return x1
#-------------------------------------------------------------------
def CumulativeDistribution(xlist, gnuplot_file = None):
    """Converts a list of numbers into a cumulative probability distribution.

    Call is: '(x, p) = CumulativeDistribution(xlist, [gnuplot_file = None])'
    This method computes the cumulative probability distribution.
    'xlist' is a list of numbers drawn from some distribution.
    The returned variable is the 2-tuple '(x, p)'.  'x' and 'p' are both
        lists.  'p[i]' is the probability that a number from the distribution
        has a value less then 'x[i]'.
    'gnuplot_file' specifies the name of the file we write this data to for
        plotting.  If it is 'None', the data is not written to a file;
        if it is a string, the data is written to a file with that name."""
    total = float(len(xlist))
    x = list(xlist)
    x.sort()
    p = [i / total for i in range(len(x))]
    if gnuplot_file != None:
        Write_Columns(gnuplot_file, x, p)
    return (x, p)
#-----------------------------------------------------------------
def CumDistFromProbDist(x, p, gnuplot_file = None):
    """Converts a probability distribution into a cumulative distribution.

    Call is: '(cx, cp) = CumDistFromProbDist(x, p, [gnuplot_file = None])'
    'x' and 'p' specify a probability distribution, with 'p[x[i]]' being
        the probability that a function has a value of 'x[i]', and
        with 'x[i]' < 'x[j]' if j > i.
    The returned variable is the 2-tuple '(cx, cp)'.  'cx' and 'cp' are both
        lists.  'cp[i]' is the probability that a number from the distribution
        has a value less then 'cx[i]'.
    'gnuplot_file' specifies the name of the file we write this data to for
        plotting.  If it is 'None', the data is not written to a file;
        if it is a string, the data is written to a file with that name."""
    cx = list(x)
    cp = []
    total = 0.0
    for pi in p:
        cp.append(total)
        total += pi
    if gnuplot_file != None:
        Write_Columns(gnuplot_file, cx, cp)
    return (cx, cp)
#-------------------------------------------------------------------
def RebinProbDist(x, p, old_dx, new_dx, gnuplot_file = None):
    """Rebins a probability distribution with larger bin sizes.

    Call is: '(xn, pn) = RebinProbDist(x, p, old_dx, new_dx, 
        [gnuplot_file = None])'
    'x' and 'p' specify a probability distribution, with 'p[x[i]]' being
        the probability that a function has a value of 'x[i]', and
        with 'x[i]' < 'x[j]' if j > i.  The 'x' values are the bin centers.
    'old_dx' are the widths of the bins in the probability distribution
        specified by 'x' and 'p'.  The bin sizes must be uniform.
    'new_dx' is the bin size for the new re-binned distribution.  It
        must be larger than 'old_dx'.
    The returned variable is the 2-tuple '(xn, pn)'.  'xn' and 'pn' 
        specify a new probability distribution with bin size of 
        'new_dx'.
    'gnuplot_file' specifies the name of the file we write this data to for
        plotting.  If it is 'None', the data is not written to a file;
        if it is a string, the data is written to a file with that name."""
    if old_dx > new_dx:
        raise StatsError("'old_dx' = %f is greater than 'new_dx' = %f" % (old_dx, new_dx))
    if len(p) != len(x):
        raise StatsError("'p' and 'x' differ in length.")
    xn = []
    pn = []
    maxx = x[0] - old_dx / 2.0 + new_dx
    psum = 0.0
    for i in range(len(x)):
        if x[i] > maxx: # new bin
            xn.append(maxx - new_dx / 2.0)
            maxx += new_dx
            pn.append(psum)
            psum = 0.0
        psum += p[i]
    xn.append(maxx - new_dx / 2.0)
    pn.append(psum)
    if gnuplot_file != None:
        Write_Columns(gnuplot_file, xn, pn)
    return (xn, pn)
#------------------------------------------------------------------
def ProbLessThan(x, x1, p1):
    """Returns the probability a variable is less than something.

    Call is: 'p = ProbLessThan(x, x1, p1)'
    'x1' and 'p1' specify a cumulative probability distribution such that
        'p1[i]' is the probability that a variable is <= 'x1[i]', and
        'x1[i] < x1[i + 1]'.
    'x' is a variable.
    'p' is the probability that a variable is <= 'x'.  Linearly
        inerpolates between specified points in the cumulative distribution."""
    for i in range(len(x1)):
        if x < x1[i]:
            break
    else:
        return 1.0
    if i == 0:
        return 0.0
    slope = float(p1[i] - p1[i - 1]) / (x1[i] - x1[i - 1])
    intercept = p1[i] - slope * x1[i]
    return x * slope + intercept 
#-------------------------------------------------------------------
def EstimateNSum(xlist, n, nestimate = None):
    """Estimates the sum of n independent random variables.

    Call is: 'sumlist = EstimateNSum(xlist, n, [nestimate = None])'
    'xlist' is a list of random variables drawn from some distribution.
    We want to compute the distribution of the sum of 'n' such independent
        random variables.
    'sumlist' is returned as a list of 'nestimate' variables each of which
        is the sum of 'n' randomly selected variables from 'xlist'.
        If 'nestimate' is 'None', it is set to 'len(xlist)'.
    This is not done in any clever statistical way -- just by random
        selection of variables."""
    if nestimate == None:
        nestimate = len(xlist)
    sumlist = []
    for i in range(nestimate):
        sum = 0.0
        for j in range(n):
            sum += random.choice(xlist)
        sumlist.append(sum)
    return sumlist
#------------------------------------------------------------------
def ProbDistribution(xlist, numbins_10_90 = 10, gnuplot_file = None):
    """Converts a list of numbers into a probability distribution.

    Call is: '(x, p, dx) = ProbDistribution(xlist, [numbins_10_90 = 20, gnuplot_file = None])'
    This method estimates the probability distribution from which some
        numbers are drawn.
    'xlist' is a list of numbers drawn from some distribution.
    The returned variable is the 3-tuple '(x, p, dx)'.  'x' and 'p'
        are both lists, while 'dx' is a number.  'p[i]' is the 
        probability that a variable in the distribution is
        >= 'x[i] - dx / 2', and <= 'x[i] + dx / 2'.
    'numbins_10_90' is used to specify the bin sizes, 'dx'.  The bin size
        is chosen such that there are 'numbins_10_90' bins between
        the 10th percentile and 90th percentile items in 'xlist'.
        If you want to specify the bin width rather than the number of bins,
        set 'numbins_10_90' to the 2-tuple '("BINSIZE", dx)' where 'dx' is
        a number representing the desired binsize.
    'gnuplot_file' specifies the name of the file we write this data to for
        plotting.  If it is 'None', the data is not written to a file;
        if it is a string, the data is written to a file with that name."""
    if isinstance(numbins_10_90, tuple) and numbins_10_90[0] == 'BINSIZE':
        (bincounts, bincenter, binmin, binmax) = Make_Bins(xlist, binsize = numbins_10_90[1])
    else:
        (bincounts, bincenter, binmin, binmax) = Make_Bins(xlist, numbins_10_90 = numbins_10_90)
    x = []
    p = []
    total = float(len(xlist))
    dx = binmax[0] - binmin[0]
    for i in range(len(bincounts)):
        p.append(bincounts[i] / total)
        x.append(bincenter[i])
    if gnuplot_file != None:
        Write_Columns(gnuplot_file, x, p)
    return (x, p, dx)
#----------------------------------------------------------------
def SumOfVariables(x, p, n, gnuplot_file = None):
    """Computes the sum of n independent variables from the same distribution.

    Call is: '(xsum, psum) = SumOfVariables(x, p, n)'
    'x' and 'p' specify a probability distribution, with 'p[x[i]]' being
        the probability that a function has a value of 'x[i]', and
        with 'x[i] - x[i + 1] = dx' where 'dx' is a constant > 0.
    'n' is the number of variables that are added.
    The returned 2-tuple, '(xsum, psum)', gives the probability distribution
        for a variable that is the sum of 'n' variables drawn from the 
        distribution specified by 'x' and 'p'.  'psum[xsum[i]]' is the
        probability that the sum has the value 'xsum[i]'.
    'gnuplot_file' specifies the name of the file we write this data to for
        plotting.  If it is 'None', the data is not written to a file;
        if it is a string, the data is written to a file with that name."""
    if not (isinstance(n, int) and n > 0):
        raise StatsError("Invalid 'n' of %r." % n)
    if not (len(x) == len(p) > 1):
        raise StatsError("'x' and 'p' differ in length or are too short.")
    dx = x[1] - x[0]
    minx = n * x[0]
    maxx = n * x[-1]
    px = [0.0] * n * len(x)
    for i in range(len(x)):
        px[i] = p[i]
    fx = numarray.fft.fft(px)
    fxn = numarray.array(fx)
    for i in range(1, n):
        fxn = fxn * fx
    pxn = numarray.fft.inverse_fft(fxn)
    psum = list(pxn.real)
    xsum = [minx + dx * i for i in range(n * len(x))]
    if gnuplot_file != None:
        Write_Columns(gnuplot_file, xsum, psum)
    return (xsum, psum)
#---------------------------------------------------------------
def KolmogorovSmirnovDistribution(data, cumfunc):
    """Is a set of data consistent with a distribution function?

    Call is: '(d, p) = KolmogorovSmirnov(data, cumfunc)'
    Based on the Numerical Recipes algorithm "ksone".
    'data' is a list of numbers drawn from the sample distribution.
    'cumfunc' is the cumulative distribution function of the
        of the distribution we think our sample is drawn from.
        'cumfunc(x)' is the probability that a point drawn from
        the sample distribution is less than or equal to 'x'.
        So 'cumfunc' should have a minimum value of 0 and a maximum
        value of one.
    'd' is returned as the maximum distance between the distributions.
    'p' is returned as the probability that we would observe a value
        of 'd' this large if the sample was drawn from the distribution."""
    datalist = list(data) # make a copy of the data
    datalist.sort() # sort from smallest to largest
    if not (0.0 <= cumfunc(datalist[0]) <= cumfunc(datalist[-1]) <= 1.0):
        raise StatsError("Invalid cumulative distribution function.")
    n = float(len(datalist)) # number of data points
    fo = d = 0.0
    for i in range(int(n)):
        fn = (i + 1) / n
        xi = datalist[i]
        ff = cumfunc(xi)
        dt = max(math.fabs(fo - ff), math.fabs(fn - ff))
        d = max(d, dt)
        fo = fn
    n = math.sqrt(n)
    p = termbf = 0.0
    fac = 2.0
    eps1 = 0.001
    eps2 = 1.0e-8
    alam = (n + 0.12 + 0.11 / n) * d
    a2 = -2.0 * alam * alam
    for j in range(1, 101):
        term = fac * math.exp(a2 * j * j)
        p += term
        if math.fabs(term) <= eps1 * termbf or math.fabs(term) <= eps2 * p:
            return (d, p)
        fac = - fac
    return (d, p)
#-------------------------------------------------------------------
def KolmogorovSmirnov(x1, x2):
    """Are two distributions different?

    Call is: '(d, p) = KolmogorovSmirnov(x1, x2)'
    Applies the Kolmogorov Smirnov test to determine if 'x1' and 'x2' are
        drawn from different distributions.
    'x1' and 'x2' are lists of numbers drawn from distributions 1 and 2.
    The returned value is the 2-tuple '(d, p)'.  'd' is the maximum distance
        between the two cumulative probability distributions.  'p' is 
        the probability that we would observe this value of 'd' if the
        numbers were drawn from the same underlying distribution.
    Based on Numerical Recipes, 14.3."""
    y1 = list(x1)
    y2 = list(x2)
    y1.sort()
    y2.sort()
    n1 = float(len(y1))
    n2 = float(len(y2))
    d = fn1 = fn2 = 0.0
    i1 = i2 = 1 
    # find 'd'
    while i1 <= n1 and i2 <= n2:
        dy = y1[i1 - 1] - y2[i2 - 1]
        if dy <= 0: 
            fn1 = i1 / n1
            i1 += 1
        if dy >= 0:
            fn2 = i2 / n2
            i2 += 1
        dt = math.fabs(fn2 - fn1)
        if dt > d:
            d = dt
    # compute effective 'n'
    n = math.sqrt(n1 * n2 / (n1 + n2))
    # compute 'p'
    eps1 = 0.001
    eps2 = 1.0e-8
    fac = 2.0
    p = termbf = 0.0
    lam = (n + 0.12 + 0.11 / n) * d
    a2 = -2.0 * lam * lam
    for j in range(1, 101):
        term = fac * math.exp(a2 * j * j)
        p += term
        if math.fabs(term) <= eps1 * termbf or math.fabs(term) <= eps2 * p:
            return (d, p)
        fac = - fac
        termbf = math.fabs(term)
    return (d, 1.0)
#------------------------------------------------------------------
def SumOfGaussiansProb(x, n, mu, sigma):
    """Gets the probability a number is the sum of Gaussian variables.

    Call is: 'p = SumOfGaussiansProb(x, n, mu, sigma)'
    This method computes the probability that the sum of 'n' random
        variables drawn from a Gaussian distribution with mean 'mu'
        and standard deviation 'sigma' has the value of 'x'.
    The returned probability 'p' is the probability that
        the sum of 'n' random variables from this Gaussian distribution
        would be some number 'y' such that:
        y >= |x - n * mu|."""
    if not (isinstance(n, int) and n > 0):
        raise StatsError("Invalid 'n' of %r." % n)
    d = math.fabs(x - n * mu) # distance from expected mean
    z = d / math.sqrt(2.0 * n * sigma * sigma)
    p = Complementary_Error_Function(z)
    return p
#--------------------------------------------------------------------
def StandardError(xlist):
    """Computes the standard error on the estimate of the mean.

    Call is: 'se = StandardError(xlist)'
    'xlist' is a list of numbers.
    'se' is returned as the standard error on the estimate of the mean."""
    m = Mean(xlist)
    se = 0.0
    for x in xlist:
        se += (x - m)**2
    se = se / float(len(xlist))
    se = math.sqrt(se)
    se /= math.sqrt(float(len(xlist)))
    return se
#-------------------------------------------------------------------
def QuotientError(x, x_error, y, y_error):
    """Propagates the error of a quotient.

    Call is: '(r, r_error) = QuotientError(x, x_error, y, y_error)'
    'r' is computed as the ratio 'x / y'.
    'r_error' is the error in 'r' given that the error in 'x' and 'y'
        are 'x_error' and 'y_error' respectively."""
    try:
        r = x / y
    except ZeroDivisionError:
        raise StatsError("'y' is zero.")
    r_error = math.sqrt(x_error**2 + r**2 * y_error**2) / float(y)
    return (r, r_error)
#----------------------------------------------------------------
def MinimizeLeastSquares(func, params, xvalues, yvalues, yerrors):
    """Finds the parameter value minimizing the least squares error.

    Call is: 'p = MinimizeLeastSquares(func, params, xvalues, yvalues,
        yerrors)'
    'func' is a function that takes as input two numbers and returns
        another single number: 'y = func(x, param)' where 'x' is the
        independent variable, 'param' is a function parameter, and 'y'
        is returned dependent variable.
    'params' is a list of possible parameter values.
    'xvalues', 'yvalues', and 'yerrors' are all lists of the same
        length.  'yvalues[i]' is the dependent variable corresponding
        to the independent variable 'xvalues[i]', and 'yerrors[i]' is
        the error on 'yvalues[i]'.
    'p' is returned as the single parameter value in 'params' that gives
        the smallest error when using 'func' to map 'xvalues' to 'yvalues'
        with the errors 'yerrors'.
    If a value of 'yerrors' is zero and if the function does not exactly
        hit this 'yvalues' point, an error is raised."""
    # error check on input variables
    if not (len(xvalues) == len(yvalues) == len(yerrors)):
        raise StatsError("Lists are of different lengths.")
    errors = []
    for param in params:
        e = 0.0
        for i in range(len(xvalues)):
            x = xvalues[i]
            y = yvalues[i]
            var = yerrors[i]**2
            fx = func(x, param)
            try:
                e += (fx - y)**2 / var
            except ZeroDivisionError:
                if fx == y:
                    continue
                else:
                    raise StatsError("Mismatch for zero error value.")
        errors.append(e)
    p = params[0]
    e = errors[0]
    for i in range(len(errors)):
        if errors[i] < e:
            e = errors[i]
            p = params[i]
    return p
#-------------------------------------------------------------------
def MakeGnuplotGraph(plotfile, outfile = None, extension = ".eps", fileformat = "post eps color", datafile = None, title = None, label_dict = {}, horizontal_lines = [], arrows = []):
    """Makes a gnuplot file.

    'plotfile' is the name of an existing gnuplot file that is used as
        the template for plotting.  For example, if the plotting
        specifications are in "plot.gnu", this should be the value
        of 'plotfile'.
    'outfile' is the name of the output file.  For example, "myplot.eps"
        if we are creating an EPS file with the name "myplot.eps"
        If 'outfile' is 'None', then it is just 'plotfile' with
        the extension changed to 'extension'.  If it is not 'None',
        we also create a new plotting file (".gnu") with the same
        base as 'outfile'.
    'extension' specifies the default extension of 'outfile'
    "fileformat" is the format of the output file. as a valid format
        string that can be passed to gnuplot as "set term".
    'datafile' specifies that we change the name of the data file
        in 'plotfile' to the value given by the string 'datafile'.
        If there are multiple data files, only replaces the name
        of the first one.
    'title' specifies that make the title of the plot into the
        value given by the 'title' string.
    'label_dict' is a dictionary.  If 'label_dict[tag] = label',
        then any labels in 'plotfile' with the value 'tag' are
        changed to have the value 'label'.
        If a tag indicated in 'label_dict' is not found, nothing happens.
    'horizontal_lines' is a list indicating that horizontal lines
        are drawn at all of the y-values in the list.
    'arrows' is a list indicating that we draw arrows.  Each entry in
        the list specifies an arrow.  The entries are 6-tuples where
        the first two items are the X and Y coordinates of the arrow's
        origin, the next two items are the X and Y coordinates of 
        the arrow's end, and the next item is the line type, and the
        next item is the line width."""
    if not os.path.isfile(plotfile):
        raise StatsError("Cannot find 'plotfile' %r." % plotfile)
    if outfile == None:
        (base, ext) = os.path.splitext(plotfile)
        outfile = "%s%s" % (base, extension)
    else:
        (base, ext) = os.path.splitext(outfile)
        file = open(plotfile)
        text = file.read()
        file.close()
        plotfile = "%s.gnu" % base
        if title: # replace title
            titlematch = re.compile('\nset title "([^"]*)" ')
            m = titlematch.search(text)
            oldtitle = m.group(1)
            text = text.replace('\nset title "%s"' % oldtitle, '\nset title "%s"' % title)
        if datafile: # replace data file
            datafilematch = re.compile('\nplot "([^"]*)" ')
            m = datafilematch.search(text)
            olddatafile = m.group(1)
            i = text.index('\nplot "%s" ' % olddatafile)
            text = "%s%s" % (text[ : i], text[i : ].replace(olddatafile, datafile))
        if horizontal_lines: # add horizontal lines
            l = [", f(x) = %f, f(x)" % y for y in horizontal_lines]
            l.append("\n#    EOF")
            s = ''.join(l)
            text = text.replace("\n#    EOF", s)
        if arrows: # add arrows
            l = ["set arrow from %f, %f to %f, %f lt %d lw %f\n" % tup for tup in arrows]
            s = ''.join(["unset arrow\n"] + l)
            text = text.replace("unset arrow\n", s)
        for (tag, label) in label_dict.iteritems():
            text = text.replace(tag, label)
        file = open(plotfile, 'w')
        file.write(text)
        file.close()
    (fi, foe) = os.popen4("gnuplot", 'w')
    fi.write('set term %s; set output "%s"; load "%s"\n' % (fileformat, outfile, plotfile))
    fi.close()
    oe = foe.read()
    foe.close()
    if oe:
        raise StatsError("Gnuplot error of %s." % oe)
#--------------------------------------------------------------------------
# End stats.py
