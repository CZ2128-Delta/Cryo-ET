import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

'''
    class RobustFit: contains whole procedure of robust fitting.

    1. Compute median residual m(p) for each view p:
        1) Compute the median residual at each view;
        2) Smooth by three iterations of fitting a quadratic locally at each point.

    2. Set up weight groups.

    3. Iterative procedure:
        1) Compute MADN (normalized median absolute deviation from the median):
        
            MADN = median(|E(i) - median(E)|) / 0.6745      E(i) = e(i) / m(p)

                e(i): residual of i-th point
                m(p): median residual of p-th view where point i is located
            
        2) Compute weight w(i) for each point i:

            d(i) = (e(i) / m(p) - m_a) / (k * MADN)
            w(i) = (1 - d(i)^2)^2, for 0 <= d(i) <= 1,
                    1,             for d(i) < 0,
                    0,             for d(i) > 1.
            
                m_a: median of the values of e(i)/m(p)
                k  : a factor that determines how extreme a point must be to 
                     receive a weight of zero

    4. Output info and return weights.
'''
class RobustFit:
    '''
        xresid: residuals along x-axis of all points
        yresid: residuals along y-axis of all points
        robustFactor: k
    '''
    def __init__(self, xreal: dict, yreal: dict, xpred: dict, ypred: dict, 
                 views: np.ndarray, numFiducial: int, numView: int, numPoint: int, 
                 robustFactor: float=4.685, pixelSize: float=1, imageBinned: int=1) -> None:
        # Inputs
        self.xreal = xreal
        self.yreal = yreal
        self.xpred = xpred
        self.ypred = ypred
        self.views = views

        self.numFiducial = numFiducial
        self.numView = numView
        self.numPoint = numPoint

        self.k = robustFactor
        self.pixelSize = pixelSize
        self.imageBinned = imageBinned

        self.xresid = dict()
        self.yresid = dict()
        self.xrealList, self.yrealList = np.zeros(self.numPoint), np.zeros(self.numPoint)
        self.xpredList, self.ypredList = np.zeros(self.numPoint), np.zeros(self.numPoint)
        self.xresidList, self.yresidList = np.zeros(self.numPoint), np.zeros(self.numPoint)

        p = 0
        for i in range(self.numView):
            ikey = f"view_{i}"
            self.xresid[ikey] = dict()
            self.yresid[ikey] = dict()
            for jkey in self.xreal[ikey].keys():
                self.xresid[ikey][jkey] = xpred[ikey][jkey] - xreal[ikey][jkey]
                self.yresid[ikey][jkey] = ypred[ikey][jkey] - yreal[ikey][jkey]

                self.xrealList[p] = self.xreal[ikey][jkey]
                self.yrealList[p] = self.yreal[ikey][jkey]
                self.xpredList[p] = self.xpred[ikey][jkey]
                self.ypredList[p] = self.ypred[ikey][jkey]
                self.xresidList[p] = self.xresid[ikey][jkey]
                self.yresidList[p] = self.yresid[ikey][jkey]

                p += 1

        # Saves
        self.groups = np.zeros(self.numPoint)
        self.weights = np.ones(self.numPoint)
    
    '''
        Compute median residual m(p) for each view p.
        Return: viewMedianResid.
    '''
    def findMedianResidual(self) -> np.ndarray:
        viewMedianResid = np.zeros(self.numView)        

        # Find medians of views
        for i in range(self.numView):
            ikey = f"view_{i}"
            tmpResid = np.zeros(0)
            for jkey in self.xresid[ikey].keys():
                tmpResid = np.append(tmpResid, (self.xresid[ikey][jkey]**2 + 
                                                self.yresid[ikey][jkey]**2)**0.5)
            viewMedianResid[i] = np.median(tmpResid)

        # Smooth with iterations
        xvfit = [0]*25
        order = 2
        nFullFit = 15
        numIter = 3

        for _ in range(numIter):
            tmpMedian = viewMedianResid.copy()
            for iv in range(1, self.numView + 1):
                ivst = int(max(1, iv - nFullFit / 2))
                ivnd = int(min(self.numView, iv + nFullFit / 2))
                nfit = ivnd + 1 - ivst
                if nfit < 5:
                    order = 1
                for i in range(ivst, ivnd + 1):
                    xvfit[i + 1 - ivst - 1] = i - iv
                coeff = np.polyfit(xvfit[0:nfit], tmpMedian[ivst-1:ivnd], order)
                viewMedianResid[iv-1] = coeff[-1]

        return viewMedianResid
    
    '''
        Sets up weighting groups, consisting of sets of views divided into rings 
        if possible.
        TODO: implement it.
        Save: self.groups.
        Print: info.
        Return: tooFewPoints.
    '''
    def setupWeightGroups(self) -> bool:
        maxViewsForRings = [1, 10, 8, 6, 5, 5, 4, 4, 3, 3]
        irealRingList = np.zeros(self.numPoint)
        distReal = np.zeros(self.numPoint)

        tooFewPoints = 0
        return tooFewPoints

    '''
        Computes the normalized median absolute deviation from the median for the
        values in [arr]. The median absolute deviation is divided by .6745, which 
        makes it estimate sigma for a normal distribution (Wilcox, R.R., Robust 
        Estimation and Hypthesis Testing, 2nd edition, P. 78).
        Return: MADN.
    '''
    def computeMADN(self, arr: np.ndarray) -> float:
        tmp = np.zeros(len(arr))
        for i in range(len(arr)):
            tmp[i] = abs(arr[i] - np.median(arr))
        MADN = np.median(tmp) / 0.6745
        return MADN
    
    '''
        Computes the weights given the current set of residuals and the weighting groups.
        TODO: implement multi-group version. Currently one group as a whole.
        Save: self.weights.
    '''
    def computeWeights(self, viewMedianResid: np.ndarray) -> None:
        residuals = np.zeros(self.numPoint)

        # Loop on the views in the group and get the residuals divided by the smoothed 
        # median residual, namely e(i) / m(p).
        p = 0
        for i in range(self.numView):
            ikey = f"view_{i}"
            for jkey in self.xresid[ikey].keys():
                residuals[p] = (self.xresid[ikey][jkey]**2 + 
                                self.yresid[ikey][jkey]**2)**0.5 / viewMedianResid[i]
                p += 1
        
        # Get group overall median (m_a) and MADN, then compute weights.
        groupMedian = np.median(residuals)
        groupMADN = self.computeMADN(residuals)
        p = 0
        for i in range(self.numView):
            ikey = f"view_{i}"
            for jkey in self.xresid[ikey].keys():
                deviation = (residuals[p] - groupMedian) / (self.k * groupMADN)
                if deviation < 0:
                    self.weights[p] = 1
                elif deviation > 1:
                    self.weights[p] = 0
                else:
                    self.weights[p] = (1 - deviation ** 2) ** 2
                p += 1
        
        return
    
    '''
        Main procedure of robust fitting.
        Save: self.weights.
        Print: info.
    '''
    def process(self) -> None:
        # Initialization for computation
        maxCycles = 3000
        wgtPrev = np.ones(self.numPoint)
        errSave = np.ones(self.numPoint)
        deltaWgtMeanCrit, deltaWgtMaxCrit, maxDeltaWgtBelowCrit = 0.001, 0.01, 4

        ############### Test ###############
        yfit = self.ypredList
        yfit2 = self.ypredList

        # xfit = self.views.copy()
        # yfit = self.ypredList.copy()
        # xfitNoWgt = self.views.copy()
        # yfitNoWgt = self.ypredList.copy()
        ####################################

        fFinal = 0
        for i in range(self.numPoint):
            fFinal += (self.xresidList[i] ** 2 + 
                       self.yresidList[i] ** 2) * self.weights[i]

        # Initialization for output
        scaleXY = 0

        pixelUnits = "nm"
        pixelSize = self.pixelSize * self.imageBinned

        errSum, errSqSum, wgtErrSum, wgtSum = 0, 0, 0, 0

        # Compute m(p) and set up weight groups
        viewMedianResid = self.findMedianResidual()
        tooFewPoints = self.setupWeightGroups()
        if tooFewPoints:
            print("WARNING: Too few data points to do robust fitting.\n")
            return

        # Iterative procedure
        numCycles = 0
        numBelowCrit = 0
        fFinal *= 10
        fLast = fFinal
        ifAveraged = 1
        while numCycles < maxCycles and not tooFewPoints:
            fPrev = fLast
            fLast = fFinal
            wgtPrev = errSave.copy()
            errSave = self.weights.copy()
            # print(errSave)
            
            self.computeWeights(viewMedianResid)

            ############### Test ###############
            p = np.poly1d(np.polyfit(self.views, yfit, 1, w=self.weights))
            p2 = np.poly1d(np.polyfit(self.views, yfit2, 1))
            yfit = p(self.views)
            yfit2 = p2(self.views)
            self.yresidList = yfit - self.ypredList
            iw = 0
            for i in range(self.numView):
                ikey = f"view_{i}"
                for jkey in self.xresid[ikey].keys():
                    self.yresid[ikey][jkey] = yfit[iw] - self.ypredList[iw]
                    iw += 1

            # def fNoWgt(c, x, y):
            #     Ri = np.sqrt((x-c[0])**2 + (y-c[1])**2)
            #     return np.square(Ri - Ri.mean())
            
            # def f(c, x, y):
            #     Ri = np.sqrt((x-c[0])**2 + (y-c[1])**2)
            #     return np.square(Ri - Ri.mean()) * self.weights
            
            # xMean = np.mean(xfit)
            # yMean = np.mean(yfit)
            # cPred = xMean, yMean
            # c, _ = optimize.leastsq(f, cPred, args=(self.views, yfit))
            # r = (np.sqrt((self.views-c[0])**2 + (yfit-c[1])**2)).mean()
            # theta = np.linspace(-np.pi, np.pi, 18)
            # xfit = c[0] + r * np.cos(theta)
            # yfit = c[1] + r * np.sin(theta)
            # self.yresidList = yfit - self.ypredList
            # iw = 0
            # for i in range(self.numView):
            #     ikey = f"view_{i}"
            #     for jkey in self.xresid[ikey].keys():
            #         self.yresid[ikey][jkey] = yfit[iw] - self.ypredList[iw]
            #         iw += 1
            # # print(self.yresidList)
            # xMeanNoWgt = np.mean(xfitNoWgt)
            # yMeanNoWgt = np.mean(yfitNoWgt)
            # cPredNoWgt = xMeanNoWgt, yMeanNoWgt
            # cNoWgt, _ = optimize.leastsq(fNoWgt, cPredNoWgt, args=(self.views, yfitNoWgt))
            # rNoWgt = (np.sqrt((self.views-cNoWgt[0])**2 + (yfitNoWgt-cNoWgt[1])**2)).mean()
            # xfitNoWgt = cNoWgt[0] + rNoWgt * np.cos(theta)
            # yfitNoWgt = cNoWgt[1] + rNoWgt * np.sin(theta)
            ####################################

            for i in range(self.numPoint):
                fFinal += (self.xresidList[i] ** 2 + 
                           self.yresidList[i] ** 2) * self.weights[i]
            numCycles += 1

            # Count the weights below various levels and analyze change in weights
            nw0, nw1, nw2, nw5, nwTotal = 0, 0, 0, 0, 0
            errMean, errSd, prevMean, prevMax, allMean = 0, 0, 0, 0, 0
            # errMean: Mean change of weights < 0.5
            # errSd: Max change in all weights
            # prevMean: Mean difference between new weight and previous weights
            # prevMax: Max difference between new weight and previous weights
            # allMean: Mean change of all weights

            for iw in range(self.numPoint):
                nwTotal += 1
                if self.weights[iw] < 0.5:
                    nw5 += 1
                    errMean += abs(self.weights[iw] - errSave[iw])
                # print(self.weights)
                # print("\n#####################\n")
                # print(errSave)
                prevMean += abs(self.weights[iw] - wgtPrev[iw])
                allMean += abs(self.weights[iw] - errSave[iw])
                if self.weights[iw] == 0:
                    nw0 += 1
                if self.weights[iw] < 0.1:
                    nw1 += 1
                if self.weights[iw] < 0.2:
                    nw2 += 1
                errSd = max(errSd, abs(self.weights[iw] - errSave[iw]))
                prevMax = max(prevMax, abs(self.weights[iw] - wgtPrev[iw]))
                # print(errMean, errSd)

            if nw5 > 0:
                errMean /= nw5
            prevMean /= self.numPoint
            allMean /= self.numPoint

            if errMean < deltaWgtMeanCrit or errSd < deltaWgtMaxCrit:
                numBelowCrit += 1
            else:
                numBelowCrit = 0

            ############### Test ###############
            iw = 0
            wgtErrSum1, wgtSum1, errSum1 = 0, 0, 0
            for i in range(self.numView):
                ikey = f"view_{i}"
                for jkey in self.xresid[ikey].keys():
                    residErr1 = (self.xresid[ikey][jkey]**2 + self.yresid[ikey][jkey]**2)**0.5
                    
                    wgtErrSum1 += (self.weights[iw] ** 0.5) * residErr1
                    wgtSum1 += self.weights[iw] ** 0.5
                    errSum1 += residErr1

                    iw += 1
            print(errSum1/self.numPoint, wgtErrSum1/wgtSum1)
            # print(self.weights)
            ####################################
            
            # print(numCycles, ifAveraged, errMean, errSd, numBelowCrit)
            if (ifAveraged == 1 and errMean < deltaWgtMeanCrit and 
                errSd < deltaWgtMaxCrit) or numBelowCrit >= maxDeltaWgtBelowCrit:
                break

            if ifAveraged == 0 and ((prevMean < allMean and prevMax < errSd) or 
                                    (fFinal / fLast > 1.02 and fLast / fPrev > 1.01) or 
                                    fFinal / fLast > 1.05):
                for iw in range(self.numPoint):
                    self.weights[iw] = 0.5 * (errSave[iw] + wgtPrev[iw])
                    errSave[iw] = 1
                ifAveraged = 1
            else:
                ifAveraged = 0
        
        # Calculate output info
        for i in range(self.numPoint):
            scaleXY = max(scaleXY, abs(self.xrealList[i]))
            scaleXY = max(scaleXY, abs(self.yrealList[i]))

        rmsScale = scaleXY ** 2 / self.numPoint

        iw = 0
        for i in range(self.numView):
            ikey = f"view_{i}"
            for jkey in self.xresid[ikey].keys():
                self.xresid[ikey][jkey] *= scaleXY
                self.yresid[ikey][jkey] *= scaleXY
                residErr = (self.xresid[ikey][jkey]**2 + self.yresid[ikey][jkey]**2)**0.5
                
                wgtErrSum += (self.weights[iw] ** 0.5) * residErr
                wgtSum += self.weights[iw] ** 0.5
                errSum += residErr
                errSqSum += residErr ** 2

                iw += 1

        print(errSum/self.numPoint, wgtErrSum/wgtSum)

        # Output successful results
        print("Total cycles for robust fitting: %5d         Final F: %14.6f\n"
              %(numCycles, (fFinal * rmsScale) ** 0.5))
        print("Final mean and max weight change: %7.4f %7.4f\n"
              %(errMean, errSd))
        print("%6d weights: %4d are 0, %4d are < .1, %4d are < .2, %5d (%4.1f%%) are < .5\n"
              %(nwTotal, nw0, nw1, nw2, nw5, (100 * nw5) / nwTotal))
        
        errMeanNm = pixelSize * errSum / self.numPoint
        errSd = ((errSqSum - errSum ** 2 / self.numPoint) / (self.numPoint - 1)) ** 0.5
        wgtErrMean = pixelSize * wgtErrSum / wgtSum
        print("Residual error mean and sd:   %7.3f %7.3f %s\n"
              %(errMeanNm, pixelSize * errSd, pixelUnits))
        print("Residual error weighted mean: %7.3f         %s\n"
              %(wgtErrMean, pixelUnits))
        
        ############### Test ###############
        print(self.weights)

        plt.plot(self.views, self.ypredList, 'o', self.views, self.yrealList, '-')
        plt.plot(self.views, yfit, '-', color='green')
        plt.plot(self.views, yfit2, '-', color='red')
        plt.show()

        # plt.plot(self.views, self.ypredList, 'o', self.views, self.yrealList, '-')
        # plt.plot(xfit, yfit, '-', color='green')
        # plt.plot(xfitNoWgt, yfitNoWgt, '-', color='red')
        # plt.axis("equal")
        # plt.show()
        ####################################
        
        return

if __name__ == '__main__':
    numFiducial = 1
    numView = 12
    numPoint = 12

    viewList = np.linspace(start=1, stop=12, num=12)
    xrealList = viewList
    yrealList = np.linspace(start=1, stop=12, num=12)
    xpredList = viewList
    ypredList = np.linspace(start=1, stop=11, num=11) + np.random.normal(0, 0.1, 11)
    ypredList = np.append(ypredList, 1)

    # viewList = np.linspace(start=0, stop=11, num=12)
    # xrealList = viewList
    # yrealList = (viewList - 5.5) ** 2
    # np.random.seed(1)
    # xpredList = viewList
    # ypredList = (viewList - 5.5) ** 2 + np.random.normal(0, 4, 12)

    # xrealList = np.linspace(start=0, stop=11, num=18)
    # xpredList = np.linspace(start=0, stop=11, num=18)
    # theta = np.linspace(-np.pi, np.pi, 18)
    # viewList = np.cos(theta)
    # yrealList = np.sin(theta)
    # # np.random.seed(1)
    # ypredList = np.sin(theta) + np.random.normal(0, 0.2, 18)

    xreal, yreal = dict(), dict()
    xpred, ypred = dict(), dict()
    p = 0
    for i in range(numView):
        ikey = f"view_{i}"
        xreal[ikey], yreal[ikey] = dict(), dict()
        xpred[ikey], ypred[ikey] = dict(), dict()
        for j in range(numFiducial):
            jkey = f"marker_{j}"
            xreal[ikey][jkey] = xrealList[p]
            yreal[ikey][jkey] = yrealList[p]
            xpred[ikey][jkey] = xpredList[p]
            ypred[ikey][jkey] = ypredList[p]
            p += 1

    # print(yrealList, yreal)

    robustfit = RobustFit(xreal, yreal, xpred, ypred, viewList, numFiducial, numView, numPoint)
    robustfit.process()

    # plt.plot(viewList, ypredList, 'o', viewList, yrealList, '-')
    # plt.axis("equal")
    # plt.show()