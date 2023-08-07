import numpy as np

######################
#   robust fitting   #
######################

av = 
mx = 

'''
  Sets up weighting groups, consisting of sets of views divided into rings if possible
  maxRings = maximum number of rings allowed
  minRes = minimum number of residuals in the group
  minTiltView = view number at minimum tilt
  ierr = return error value, 0 for success or 1 for failure

  For ordinary weighting per point, the arrays mean the following:
  ivStartWgtGroup has the starting view index of the weight group
  these view indexes increment through all the groups
  ipStartWgtView has the starting projection point index for each view subset in each
  group.  These indexes also increment through all the groups
  indProjWgtList has the true 2D point index for each projection point index
'''
def setupWeightGroups(maxRings, minRes, minTiltView):
    maxViewsForRings = [1, 10, 8, 6, 5, 5, 4, 4, 3, 3]
    irealRingList = np.zeros(mx.maxReal)
    distReal = np.zeros(mx.maxReal)

    nrealForViews = av.nrealPt
    if av.patchTrackModel:
        nrealForViews = av.numFullTracksUsed
    
    # Get distances from center and get sorted indexes to them
    # irealRingList are indexes from 0 not 1
    numRealUsed = 0
    for i in range(av.nrealPt):
        if (not(av.leavingOut and av.realLeftOut[i]) and 
            not(av.realInTestSet and av.realInTestSet[i])):
            distReal[numRealUsed] = (pow(av.xyz[i * 3], 2) + 
                                     pow(av.xyz[i * 3 + 1], 2)) ** 0.5
            irealRingList[numRealUsed] = numRealUsed
            numRealUsed += 1
    rsSortIndexedFloats(distReal, irealRingList)

    # Loop from largest number of rings down, first evaluate plausibility if
    # all points are present.  Here use the number of full tracks for patch tracking
    # NUM_RING_LOOP:
    for nring in range(maxRings, 0, -1):
        ierr = 1
        nrealPerRing = nrealForViews / nring
        if nrealPerRing == 0:
            continue
        neededViews = max(1, minRes / nrealPerRing)
        if neededViews > maxViewsForRings[nring - 1] and nring > 1:
            continue
        #
        # Try to set up groups of views of increasing sizes until one works
        # NUM_VIEW_LOOP:
        exitNumView = False
        for numViews in range(neededViews, av.nview + 1):
            if exitNumView:
                break
            numViewGroups = av.nview / numViews
            numExtra = av.nview % numViews
            iexStart = max(1, minTiltView - numExtra / 2)
            ngrpBeforeEx = (iexStart - 1) / numViews
            indProj = 1
            ivbase = 0
            indGroup = 1
            indView = 1
            #
            # Loop on the groups of views
            # GROUP_LOOP:
            exitGroup = False
            for igroup in range(1, numViewGroups + 1):
                if exitGroup:
                    break
                irbase = 0
                ninGroup = numViews
                if igroup > ngrpBeforeEx and igroup <= ngrpBeforeEx + numExtra:
                    ninGroup += 1
                #
                # loop on the rings in views; here base the number on the total number of real
                # points not full tracks, because each will be considered for whether they have
                # points in the view
                # RING_LOOP:
                for iring in range(1, nring + 1):
                    if exitGroup:
                        break
                    ninRing = numRealUsed / nring
                    if iring > nring - (numRealUsed % nring):
                        ninRing += 1
                    #
                    # This is one weight group, set the starting view index of it
                    av.ivStartWgtGroup[indGroup - 1] = indView
                    indGroup += 1
                    #
                    # loop on the views in the group; for each one, set the starting index
                    # in the projection list
                    for iv in range(ivbase + 1, ivbase + ninGroup + 1):
                        av.ipStartWgtView[indView - 1] = indProj
                        indView += 1
                        #
                        # Loop on the real points in the ring, and for each one on the given view,
                        # add its projection point index to the list
                        for ind in range(irbase + 1, irbase + ninRing + 1):
                            ireal = irealRingList[ind - 1] + 1
                            for iproj in range(av.irealStr[ireal - 1], av.irealStr[ireal]):
                                if (av.isecView[iproj - 1] == iv and 
                                    not(av.leavingOut and av.projLeftOut[iproj - 1])):
                                    av.indProjWgtList[indProj - 1] = iproj
                                    indProj += 1
                    #
                    # After each ring, increase the ring base index
                    irbase += ninRing
                    #
                    # But if there are too few in this group, make view groups bigger if
                    # possible for this number of rings; otherwise go on to try fewer rings;
                    # but if there is only one ring and it is up to all views, push on
                    numPts = indProj - av.ipStartWgtView[av.ivStartWgtGroup[indGroup - 2] - 1]
                    if numPts < minRes:
                        ierr = 1
                        if numViews >= maxViewsForRings[nring - 1] and nring > 1:
                            exitGroup = True
                            exitNumView = True
                            break
                        if nring > 1 or numViews < av.nview:
                            exitGroup = True
                            break
                # End of RING_LOOP
                if exitGroup:
                    break
                #
                # After each view group, increase the view base number
                ivbase += ninGroup
                ierr = 0
            # End of GROUP_LOOP
            if exitNumView:
                break
            #
            # If we got here with err 0, this setup fits constraints, finalize index lists
            if ierr == 0:
                av.ipStartWgtView[indView - 1] = indProj
                av.ivStartWgtGroup[indGroup - 1] = indView
                av.numWgtGroups = nring * numViewGroups
                if not(av.leavingOut):
                    print("\nStarting robust fitting with", av.numWgtGroups, 
                          "weight groups:", numViewGroups, 
                          "view groups in", nring, "rings\n")
                ierr = 0
                return ierr
        # End of NUM_VIEW_LOOP
    # End of NUM_RING_LOOP
    ierr = 1
    return ierr

'''
  Computes the weights given the current set of residuals and the weighting groups.
  distRes and work are temp arrays that need to be at least as big as the biggest view
  group, and iwork needs to be as big as number of points on view.
'''
def computeWeights(indAllReal, distRes, work, iwork):
    wholeTrack = av.patchTrackModel and av.obustByTrack
    minNonZero = 4
    tooManyZeroDelta = .02
    #
    # Loop on the groups of views
    for igroup in range(1, av.numWgtGroups + 1):
        ninGroup = 0
        numViews = av.ivStartWgtGroup[igroup] - av.ivStartWgtGroup[igroup - 1]
        #
        # loop on the views in the group and get the residuals, divide by the smoothed
        # median residual
        for indv in range(av.ivStartWgtGroup[igroup - 1], av.ivStartWgtGroup[igroup]):
            ninView = av.ipStartWgtView[indv] - av.ipStartWgtView[indv - 1]
            for i in range(1, ninView + 1):
                ind = av.indProjWgtList[i + av.ipStartWgtView[indv - 1] - 1 - 1]
                if wholeTrack:
                    #
                    # recompute mean residual for track
                    av.trackResid[ind - 1] = 0.
                    for j in range(av.irealStr[ind - 1], av.irealStr[ind]):
                        if not(av.leavingOut and av.projLeftOut[j - 1]):
                            av.trackResid[ind - 1] += (pow(av.xresid[j - 1], 2) + 
                                                       pow(av.yresid[j - 1], 2)) ** 0.5
                    av.trackResid[ind - 1] /= (av.irealStr[ind] - av.irealStr[ind - 1])
                    distRes[ninGroup + i - 1] = av.trackResid[ind - 1] / \
                        av.viewMedianRes[av.itrackGroup[indAllReal[ind - 1] - 1] - 1]
                else:
                    distRes[ninGroup + i - 1] = (pow(av.xresid[ind - 1], 2) + 
                                                 pow(av.yresid[ind - 1], 2)) ** 0.5 / \
                        av.viewMedianRes[av.isecView[ind - 1] - 1]
            ninGroup += ninView
        #
        # Get overall median and MADN and compute weights
        rmedian = rsFastMedian(distRes, ninGroup, work)
        rMADN = rsFastMADN(distRes, ninGroup, rmedian, work)
        ninGroup = 0
        for indv in range(av.ivStartWgtGroup[igroup - 1], av.ivStartWgtGroup[igroup]):
            ninView = av.ipStartWgtView[indv] - av.ipStartWgtView[indv - 1]
            numLow = weightsForView(rmedian, wholeTrack, ninView, indv, distRes, 
                                    ninGroup, rMADN, tooManyZeroDelta)
            
            maxSmall = ninView * av.smallWgtMaxFrac
            if numLow > maxSmall:
                #
                # If there are too many small weights, then find the first one that must be
                # given a weight above threshold and adjust the median to accomplish that
                for i in range(1, ninView + 1):
                    iwork[i - 1] = i - 1
                rsSortIndexedFloats(distRes[ninGroup:], iwork)
                
                dev = (1. - (1.1 * av.smallWgtThreshold) ** 0.5) ** 0.5
                adjMedian = distRes[ninGroup + iwork[ninView -  maxSmall - 1]] - \
                    dev * av.kfacRobust * rMADN
                numLow = weightsForView(adjMedian, wholeTrack, ninView, indv, distRes, 
                                        ninGroup, rMADN, tooManyZeroDelta)
                if numLow > maxSmall:
                    print("WARNING: Median adjustment for small weights bad\n")
            
            ninGroup += ninView
    
    if wholeTrack:
        return
    #
    # Now look at each real point and make sure it has enough non-zero weights, or
    # specifically points above the delta value; and if not, add delta to all weights
    for ind in range(1, av.nrealPt + 1):
        if ((av.leavingOut and av.ealLeftOut[ind - 1]) or 
            (av.realInTestSet and av.realInTestSet[ind - 1])):
            continue
        numAbove = 0
        for i in range(av.irealStr[ind - 1], av.irealStr[ind]):
            if (not(av.leavingOut and av.projLeftOut[i - 1]) and 
                av.weight[i - 1] > tooManyZeroDelta):
                numAbove += 1
        if numAbove < minNonZero:
            for i in range(av.irealStr[ind - 1], av.irealStr[ind]):
                if not(av.leavingOut and av.projLeftOut[i - 1]):
                    av.weight[i - 1] = min(1., av.weight[i - 1] + tooManyZeroDelta)

'''
  Compute the weights for one view in a weight group, or for tracks in a track group.
'''
def weightsForView(viewMedian, wholeTrack, ninView, indv, distRes, 
                   ninGroup, rMADN, tooManyZeroDelta):
    numLow = 0
    if wholeTrack:
        #
        # For whole track, loop on track, get one weight for each, make sure it is not 0,
        # and assign it to all projection points
        for i in range(1, ninView + 1):
            ind = av.indProjWgtList[i + av.ipStartWgtView[indv - 1] - 1 - 1]
            dev = (distRes[ninGroup + i - 1] - viewMedian) / (av.kfacRobust * rMADN)
            if dev <= 0.:
                trackWgt = 1.
            elif dev >= 1.:
                trackWgt = tooManyZeroDelta
            else:
                trackWgt = max(tooManyZeroDelta, pow(1. - dev * dev, 2))
            
            if trackWgt < av.smallWgtThreshold:
                numLow += 1
            for jj in range(av.irealStr[ind - 1], av.irealStr[ind]):
                if not(av.leavingOut and av.projLeftOut[jj - 1]):
                    av.weight[jj - 1] = trackWgt
    else:
        #
        # Otherwise, loop on projection points in the view-weight group
        for i in range(1, ninView + 1):
            ind = av.indProjWgtList[i + av.ipStartWgtView[indv - 1] - 1 - 1]
            dev = (distRes[ninGroup + i - 1] - viewMedian) / (av.kfacRobust * rMADN)
            if dev <= 0.:
                av.weight[ind - 1] = 1.
            elif dev >= 1.:
                av.weight[ind - 1] = 0.
            else:
                av.weight[ind - 1] = pow(1. - dev * dev, 2)
            
            if av.weight[ind - 1] < av.smallWgtThreshold:
                numLow += 1
    return numLow



#########################
#   robust statistics   #
#########################

valArray = list()
indexOffset = 0

def indexedFloatCompar(val):
  global valArray, indexOffset
  i = val - indexOffset
  return valArray[i]

'''
  Uses {sort} to sort indexes in [index] to the [n] floats in the array [x].
'''
def rsSortIndexedFloats(x, index):  # (x, index, n)
    global valArray, indexOffset
    valArray = x
    index.sort(key=indexedFloatCompar)
    indexOffset = 0

'''
  Selects item number [s] (numbered from 1) out of [num] items in the array
  [r], where items are considered in order from low to high.  [r] is partially
  rearranged while finding the item.  The algorithm runs in linear time and is faster
  than sorting the array first.  Returns 0 for [num] = 0.
'''
def percentileFloat(s, r, num):
    lo = 0
    up = num - 1
    if num <= 0:
        return 0.
    if num == 1:
        return r[0]
    s = max(1, min(num, s))
    s -= 1
    while up >= s and s >= lo:
        i = lo
        j = up
        temp = r[s]
        if s > num / 2:
            #
            # Operations for high percentiles
            r[s] = r[lo]
            r[lo] = temp
            while i < j:
                while r[j] > temp:
                    j -= 1
                r[i] = r[j]
                while i < j and r[i] <= temp:
                    i += 1
                r[j] = r[i]
            
            r[i] = temp
            if s < i:
                up = i - 1
            else:
                lo = i + 1

        else:
            #
            # Operations for low percentiles
            r[s] = r[up]
            r[up] = temp
            while i < j:
                while r[i] < temp:
                    i += 1
                r[j] = r[i]
                while i < j and r[j] >= temp:
                    j -= 1
                r[i] = r[j]
            
            r[j] = temp
            if s > j:
                lo = j + 1
            else:
                up = j - 1

    return r[s]

'''
  Computes the median of the [n] values in [x] in linear time.  The value is returned in 
  [median] and [x] is rearranged by one or two calls to @@percentileFloat@.  This 
  routine is faster than @rsMedian even for small [n], and much faster for large [n].
'''
def rsFastMedianInPlace(x, n):
    median = percentileFloat((n + 1) / 2, x, n)
    if n % 2 == 0:
        median = ( median + percentileFloat(n / 2 + 1, x, n) ) / 2.
    return median

'''
  Computes the median of the [n] values in [x] in linear time.  The value is returned in 
  [median], and [xjumble] is used for calling @@rsFastMedianInPlace@.
'''
def rsFastMedian(x, n, xjumble):
    xjumble = x[0:n]
    median = rsFastMedianInPlace(xjumble, n)
    return median

'''
  Computes the normalized median absolute deviation from the median for the
  [n] values in [x] in linear time, using the value already computed for the median in 
  [median].  The median absolute deviation is divided by .6745, which makes it 
  estimate sigma for a normal distribution (Wilcox, R.R., Robust Estimation and 
  Hypthesis Testing, 2nd edition, P. 78).  The result is returned in [MADN], and 
  [tmp] is used for storing absolute deviations.
'''
def rsFastMADN(x, n, median, tmp):
    for i in range(n):
        tmp[i] = abs(x[i] - median)
    MADN = rsFastMedianInPlace(tmp, n)
    MADN /= 0.6745
    return MADN



##################
#   regression   #
##################

'''
  Uses multiple linear regression to fit a polynomial of order
  [order] to [ndata] points whose (x,y) coordinates are in the arrays [x] and [y].
  It returns the coefficient of x to the i power in the array [slopes] and a
  constant term in [intcpt].  The equation fit is:  ^
  Y = intcpt + slopes\[0\] * X + slopes\[1\] * X**2 + ...  ^
  [work] is an array whose size must be at least ([order] + 1) * ([order] + 3 + [ndata]).
  The return value is the value returned by @@multRegress@.  Note that a Fortran 
  function polyfit in libhvem takes care of allocating [work] to the needed size and 
  calling the Fortran wrapper to this function.
'''
def polynomialFit(x, y, ndata, order, slopes, intcpt, work):
    wdim = order + 1
    xMean = work + wdim * ndata
    xSD = xMean + wdim
    mwork = xSD + wdim
    if not(order):
        return 1
    for i in range(ndata):
        for j in range(order):
            work[i + j * ndata] = pow(x[i], j + 1.)
        work[i + order * ndata] = y[i]
    
    return multRegress(work, ndata, 0, order, ndata, 1, 0, slopes, ndata, intcpt, xMean,
                      xSD, mwork)

'''
  Computes a multiple linear regression (least-squares fit) for the relationships between
  one or more dependent (output) variables and a set of independent (input) variables.  ^
  Input parameters:  ^
  [x]        - Input data matrix  ^
  [xSize]    - Size of the fastest-progressing dimension of [x]  ^
  [colFast]  - Nonzero if the column dimension is the fastest progressing one, i.e. if
  successive values in the array occur in successive columns  ^
  [numInpCol]  - Number of columns of data for input variables.  This is allowed to be 
  0 ^
  [numData]   - Number of rows of data; i.e., number of separate measurements ^
  [numOutCol] - Number of columns of data for output variables; i.e., number of 
  relationships to fit  ^
  [wgtCol]    - Column number with weighting factors if > 0, otherwise no weighting.  
  Columns are numbered from 0 when calling from C, or 1 calling from Fortran. ^
  [solSize]   - Size of the fastest-progressing dimension of [sol], the array/matrix to
  receive the solutions; must be at least [numInpCol]  ^
  [work]      - Any array for temporary use whose size must be at least 
  (numInpCol + numOutCol) * (numInpCol + numOutCol)  floating point elements  ^
  Outputs:  ^
  [sol]       - Matrix to receive the [numOutCol] sets of [numInpCol] coefficients.  Each
  set is placed in a column of [sol], where the row dimension is the fastest progressing 
  one ^
  [cons]      - Array to receive the [numOutCol] constant terms of the fits, or NULL
  to fit an equation with no constant term  ^
  [xMean]     - Array for the means of the [numInpCol] plus [numOutCol] data columns  ^
  [xSD]       - Array for the standard deviations of the [numInpCol] plus [numOutCol] 
  data
  columns (these will be far from correct if [cons] is NULL)  ^
  The input variables should be placed in the first [numInpCol] columns of [x] and the
  the output variables in the next [numOutCol] columns. The data element at row i,
  column j would be accessed as x\[i + j * xSize\] or x\[j\]\[i\] from C, or as x(i,j) 
  from Fortran if [colFast] is 0, or as x\[j + i * xSize\], x\[i\]\[j\], or x(j,i) if 
  [colFast] is nonzero.  ^
  The return value is 1 if [wgtCol] has an inappropriate value, or 3 if @gaussj returns
  with an error.
'''
def multRegress(x, xSize, colFast, numInpCol, numData, 
                numOutCol, wgtCol, sol, solSize, cons,
                xMean, xSD, work):
    colStride = 1 if colFast else xSize
    rowStride = xSize if colFast else 1
    mp = numInpCol + numOutCol
    fndata = numData
    if wgtCol > 0 and (wgtCol < mp or (not(colFast) and wgtCol >= xSize)):
        return 1
    
    # Get the unweighted means
    if wgtCol <= 0:
        for i in range(mp):
            dsum = 0.
            for k in range(numData):
                dsum += x[(k) * rowStride + (i) * colStride]
            xMean[i] = dsum / fndata
    # Or add up the weights and get the weighted means
    else:
        wsum = 0.
        for k in range(numData):
            wsum += x[(k) * rowStride + (wgtCol) * colStride]
        for i in range(mp):
            dsum = 0.
            for k in range(numData):
                dsum += x[(k) * rowStride + (i) * colStride] * \
                        x[(k) * rowStride + (wgtCol) * colStride]
            xMean[i] = dsum / wsum
    
    # Get the sums of squares and cross-products of deviations
    for i in range(mp):
        for j in range(i, mp):
            if i >= numInpCol and i != j:
                continue
            dsum = 0.
            if cons:
                if wgtCol <= 0:
                    for k in range(numData):
                        dsum += (x[(k) * rowStride + (i) * colStride] - xMean[i]) * \
                                (x[(k) * rowStride + (j) * colStride] - xMean[j])
                else:
                    for k in range(numData):
                        dsum += (x[(k) * rowStride + (i) * colStride] - xMean[i]) * \
                                (x[(k) * rowStride + (j) * colStride] - xMean[j]) * \
                                x[(k) * rowStride + (wgtCol) * colStride]
            else:
                if wgtCol <= 0:
                    for k in range(numData):
                        dsum += x[(k) * rowStride + (i) * colStride] * \
                                x[(k) * rowStride + (j) * colStride]
                else:
                    for k in range(numData):
                        dsum += x[(k) * rowStride + (i) * colStride] * \
                                x[(k) * rowStride + (j) * colStride] * \
                                x[(k) * rowStride + (wgtCol) * colStride]
            work[j * mp + i] = dsum
    
    # Get the SDs
    for i in range(mp):
        xSD[i] = (work[i * mp + i] / (fndata - 1.)) ** 0.5
    
    # If numInpCol = 0, then just return the c values as the weighted means
    if not(numInpCol):
        if cons:
            for j in range(numOutCol):
                cons[j] = xMean[j]
        return 0
    
    # Scale it by n - 1 to get covariance and by SD's to get the correlation matrix
    for i in range(numInpCol):
        for j in range(i, mp):
            den = xSD[i] * xSD[j]
            if den < 1.e-30:
                work[j * mp + i] = 1.
            else:
                work[j * mp + i] /= den * (fndata - 1.)
            if j < numInpCol:
                work[i * mp + j] = work[j * mp + i]
    
    # The matrix to be solved is now completely filled and symmetric so row/column
    # doesn't matter, but to call gaussj we need to transpose the final columns into sol
    for j in range(numOutCol):
        for i in range(numInpCol):
            sol[j + i * numOutCol] = work[(j + numInpCol) * mp + i]
    
    if gaussj(work, numInpCol, mp, sol, numOutCol, numOutCol):
        return 3
    
    # Scale the coefficients and transpose them back; get constant terms
    work = sol[0:numInpCol * numOutCol]
    for j in range(numOutCol):
        if cons:
            cons[j] = xMean[numInpCol + j]
        for i in range(numInpCol):
            if xSD[i] < 1.e-30:
                sol[i + solSize * j] = 0.
            else:
                sol[i + solSize * j] = work[j + i * numOutCol] * xSD[numInpCol + j] / xSD[i]
            if cons:
                cons[j] -= sol[i + solSize * j] * xMean[i]
    return 0



####################
#   Gauss-Jordan   #
####################

'''
  Solves the linear matrix equation A X = B by Gauss-Jordan elimination.  
  A is a square matrix of size [n] by [n] in array [a], dimensioned to [np] 
  columns.  B is a matrix with one row per row of A and [m] columns in array
  [b], dimensioned to [mp] columns.  The columns of [b] are replaced by
  the [m] solution vectors while [a] is reduced to a unit matrix.  It is 
  called the same from C and Fortran, but the matrices must be in row-major 
  order in both cases.  Specifically, in C, the matrices are indexed as 
  A\[row\]\[column\] and in Fortran they are indexed as A(column,row). ^
  The maximum value of [n] is 2000 but better and faster approaches should be 
  used long before reaching that point.  The routine returns -1 if [n] exceeds
  this value and 1 if the A matrix is singular.
'''
def gaussj(a, n, npp, b, m, mp): # np -> npp to avoid confliction with short for numpy
    MSIZ = 2000
    index = np.zeros((MSIZ, 2))
    pivot = np.zeros(MSIZ)
    ipivot = np.zeros(MSIZ)

    determ = 1.
    if n > MSIZ:
        return -1
    for j in range(n):
        ipivot[j]=0
    for i in range(n):
        amax=0.
        for j in range(n):
            if ipivot[j] != 1:
                for k in range(n):
                    if ipivot[k] == 0:
                        abstmp = a[j*npp+k]
                        if abstmp < 0:
                            abstmp = -abstmp
                        if amax < abstmp:
                            irow=j
                            icolum=k
                            amax=abstmp
                    elif ipivot[k] > 1:
                        return 1
        ipivot[icolum]=ipivot[icolum]+1
        if irow != icolum:
            determ = -determ
            for l in range(n):
                t=a[irow*npp+l]
                a[irow*npp+l]=a[icolum*npp+l]
                a[icolum*npp+l]=t
            for l in range(m):
                t=b[irow*mp+l]
                b[irow*mp+l]=b[icolum*mp+l]
                b[icolum*mp+l]=t
        index[i][0]=irow
        index[i][1]=icolum
        pivotmp=a[icolum*npp+icolum]

        pivot[i]=pivotmp
        determ *= pivotmp
        a[icolum*npp+icolum]=1.

        for l in range(n):
            a[icolum*npp+l]=a[icolum*npp+l]/pivotmp
        for l in range(m):
            b[icolum*mp+l]=b[icolum*mp+l]/pivotmp
        for l1 in range(n):
            t=a[l1*npp+icolum]
            if t != 0. and l1 != icolum:
                a[l1*npp+icolum]=0.
                for l in range(n):
                    a[l1*npp+l]=a[l1*npp+l]-a[icolum*npp+l]*t
                for l in range(m):
                    b[l1*mp+l]=b[l1*mp+l]-b[icolum*mp+l]*t
    for i in range(n):
        l=n-1-i
        if index[l][0] != index[l][1]:
            irow=index[l][0]
            icolum=index[l][1]
            for k in range(n):
                t=a[k*npp+irow]
                a[k*npp+irow]=a[k*npp+icolum]
                a[k*npp+icolum]=t
    return 0



####################
#   main process   #
####################

def findMedianResidual():
    xvfit = [0]*25
    slope = [0]*5
    tmpRes = np.zeros(0)
    tmpMedian = np.zeros(0)
    iorder = 2
    nFullFit = 15
    numIter = 3
    allZero = True
    iterCount = 0

    # Set up number of medians to find and smoothing iterations
    if av.patchTrackModel and av.robustByTrack:
        av.viewMedianRes[0] = 1.
        if av.numTrackGroups == 1:
            return
        numMedian = av.numTrackGroups
        numIter = 0
    else:
        numMedian = av.nview
    tmpMedian.resize(numMedian)

    # Find medians of either track groups or views
    for iv in range(1, numMedian + 1):
        ninViewSum = 0
        tmpRes.resize(0)
        for i in range(av.nrealPt):
            if ((av.realInTestSet and av.realInTestSet[i]) or 
                (av.leavingOut and av.realLeftOut[i])):
                continue
            #
            # Track group: compute mean residual of each track and save it
            if av.patchTrackModel and av.robustByTrack:
                pass
            #
            # View: get residual of each point
            else:
                for j in range(av.irealStr[i] - 1, av.irealStr[i + 1] - 1):
                    if av.isecView[j] == iv and not(av.leavingOut and av.projLeftOut[j]):
                        ninViewSum += 1
                        tmpRes = np.append(tmpRes, (pow(av.xresid[j], 2.) + pow(av.yresid[j], 2.))**0.5)
        av.viewMedianRes[iv - 1] = rsFastMedianInPlace(tmpRes, ninViewSum)
        if av.viewMedianRes[iv - 1] != 0:
            allZero = False
    iterCount += 1
    if allZero:
        for i in range(numMedian):
            av.viewMedianRes[i] = 1.
        return
    #
    # Make sure there are no zeros, just copy nearest value
    for iv in range(1, numMedian + 1):
        if av.viewMedianRes[iv - 1] == 0.:
            found = False
            # NEAR_VIEW_LOOP:
            for i in range(1, numMedian):
                if found:
                    break
                for iter in range(-1, 2, 2):
                    j = iv + iter * i
                    if j > 0 and j <= numMedian:
                        if av.viewMedianRes[j - 1] != 0.:
                            av.viewMedianRes[iv - 1] = av.viewMedianRes[j - 1]
                            found = True
                            break
    #
    # smooth with iterations
    tmpRes.resize((iorder + 1) * (iorder + 6 + numMedian))
    for iter in range(1, numIter + 1):
        tmpMedian = av.viewMedianRes[0:numMedian]
        for iv in range(1, numMedian + 1):
            ivst = max(1, iv - nFullFit / 2)
            ivnd = min(numMedian, iv + nFullFit / 2)
            nfit = ivnd + 1 - ivst
            if nfit < 5:
                iorder = 1
            for i in range(ivst, ivnd + 1):
                xvfit[i + 1 - ivst - 1] = i - iv
            polynomialFit(xvfit, tmpMedian[ivst - 1:], nfit, iorder, slope,
                          av.viewMedianRes[iv - 1:], tmpRes)
        iterCount += 1
    
    return

if __name__ == '__main__':
    mMaxCycles = 1000
    mRobustTotCycleFac = 3.
    mMinResRobust = 100
    maxWgtRings = 10

    robFailed = False
    
    maxTotCycles = abs(mMaxCycles) * mRobustTotCycleFac
    numTotCycles = 0
    numOneCycle = 0
    numBelowCrit = 0
    jpt = mMinResRobust
    index = maxWgtRings

    mNumWgtTotal = 0
    mNumWgtZero = 0
    mNumWgt1 = 0
    mNumWgt2 = 0
    mNumWgt5 = 0

    # @cz: Some vars are initialized in input_vars() in original C codes; 
    # here we initialize directly (maybe temporarily)
    mMinTiltView = 
    mNumProjPt = 
    mErrSave = 
    mVar = 
    mNvarSearch = 
    mIndAllReal = 
    mIndSave = 
    mJptSave = 
    mXyzErr = 

    fFinal = ?? # fFinal is the final error measure

    findMedianResidual()
    # @cz: We use setupWeightGroups instead of setupTrackWeightGroups
    ierr = setupWeightGroups(index, jpt, mMinTiltView)
    robTooFew = ierr != 0
    if robTooFew and av.leavingOut:
        print("WARNING: Too few data points to do robust fitting\n")
    av.robustWeights = 1
    for jpt in range(mNumProjPt):
        if not(av.leavingOut) and av.testSetFracStep <= 0:
            av.weight[jpt] = 1.
        mErrSave[jpt] = 1.
    fOriginal = fFinal
    fFinal *= 10.
    fLast = fFinal
    ifAveraged = 1
    metroRobust = 0
    mVarSave = mVar[0:mNvarSearch]
    while numTotCycles < maxTotCycles and not(robTooFew):
        fPrevious = fLast
        fLast = fFinal
        mWgtPrev = mErrSave[0:mNumProjPt]
        mErrSave = av.weight[0:mNumProjPt]
        computeWeights(mIndAllReal, mIndSave, mJptSave, mXyzErr)
        