# -*- coding:utf-8 -*-
#!/usr/bin/env python
""" Классы данных спектра и визуализации
Created on 22.12.2016, @author: dmagin
"""
import numpy as NP
from numpy import linalg as LA

import Geometry as Geo

import matplotlib.pyplot as plt #, mpld3
import matplotlib.tri as tri

def mR2G(mR):
    # Преобразование матрицы резистенсов в матрицу Грина.
    # Или, что тоже самое, - матрицы квадратов расстояний в матрицу корреляции
    # G = -(eL R eL)/2
    N = len(mR)
    eL = NP.eye(N) - NP.ones(N)/N
    return NP.dot(NP.dot(eL, mR), eL)/(-2.)

def Agregat(vArray, sFunction='Var', axis=0):
    # Считает разные агрегаты
    if sFunction == 'Var':
        return vArray.var(axis=axis)
    elif sFunction == 'Std':  
        return vArray.std(axis=axis)
    elif sFunction == 'Mean':  
        return vArray.mean(axis=axis)
    else:
        return None

def dicMarkStyle():
    # Словарь стилей маркера
    # http://matplotlib.org/api/markers_api.html#module-matplotlib.markers
    #return {v:k for k, v in Markers.MarkerStyle().markers.items()}
    dictMS = {'o': 'circle', 'h': 'hexagon', 'D': 'diamond', '*': 'star', '.': 'point', ',': 'pixel', 'None': 'nothing'}
    return {v:k for k, v in dictMS.items()}

def mColorMap():
    # перечень доступных цветовых карт
    # http://matplotlib.org/examples/color/colormaps_reference.html
    # ['hsv', 'gist_rainbow', 'jet', 'RdBu_r', 'gist_earth', 'rainbow', 'brg', 'winter', 'prism']
    cmaps = [#'== Perceptually Uniform Sequential ==',
                'viridis', 'inferno', 'plasma', 'magma',
             #'== Sequential ==',
                'Blues', 'BuGn', 'BuPu',
                'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd',
             #'== Sequential (2) ==',
                'afmhot', 'autumn', 'bone', 'cool',
                'copper', 'gist_heat', 'gray', 'hot',
                'pink', 'spring', 'summer', 'winter',
             #'== Diverging ==',
                'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'seismic',
             #'== Qualitative ==',
                 'Accent', 'Dark2', 'Paired', 'Pastel1',
                 'Pastel2', 'Set1', 'Set2', 'Set3',
             #'== Miscellaneous ==',
                'gist_earth', 'terrain', 'ocean', 'gist_stern',
                'brg', 'CMRmap', 'cubehelix',
                'gnuplot', 'gnuplot2', 'gist_ncar',
                'nipy_spectral', 'jet', 'rainbow',
                'gist_rainbow', 'hsv', 'flag', 'prism']
    
    return cmaps

class DataSpectr(object):
    '''Класс данных спектра'''
    
    def Spectr2Projections(self, sVal, mVector):
        '''Получаем вырожденные 2D-проекции спектра'''
        self.vIndex = [] # индексы вырожденных спектров
        self.vS = [] # собственные числа вырожденных спектров
        self.vDataX, self.vDataY = [], [] # собственные функции вырожденных спектров
        self.vDataZ = [] # собственные функции невырожденных спектров
         
        NumDims = len(sVal)
        self.numSp = NumDims - 1 # Всего уровней
        self.numSp2, self.numZ = 0, 0 # Вырожденных и невырожденных уровней
        bDegen = False # Признак вырожденного уровня 
        for i in range(0, NumDims-1):
            if bDegen:
                bDegen = False
                continue
            sV = sVal[i]
            if abs(sV) < self.error:
                continue # исключаем нулевой уровень Должен быть один
            if abs((sV - sVal[i+1])/sV) < self.error: # вырожденный уровень
                self.vIndex.append(i)
                self.vDataX.append(mVector[i])
                self.vDataY.append(mVector[i+1])
                self.vS.append(sV)
                self.numSp2 += 1
                bDegen = True
            else:
                bDegen = False
                self.vDataZ.append(mVector[i])
                self.numZ += 1
                if i == NumDims-2: # Последний уровень
                    self.vDataZ.append(mVector[i+1])
                    self.numZ += 1
   
    def SetFunction(self):
        # Добавление возмущения к исходным данным
        # F(R2) = w*Rd + 1/Rd, где w = dist^3/n^2, Rd = R2^degree
        # (в матричном виде ругается на деление на ноль при отрицательных степенях из-за нулевой диагонали)
        size = len(self.mR2)
        weight = self.distortion**3/(self.Size**2)
        mRD = NP.zeros((size, size))
        lr = range(size)
        for i in lr:
            for j in range(i+1, size):
                r2 = self.mR2[i][j]
                Rd = r2**self.degree
                mRD[i][j] = weight*Rd + 1/Rd
                mRD[j][i] = mRD[i][j] 
        self.mFunction = mRD 

    def SetSpectrData(self):
        self.SetFunction()
        mG = mR2G(self.mFunction)
        sG, vG = LA.eigh(mG)
        vG = NP.transpose(vG)
        self.Spectr2Projections(sG, vG)
        #print(sG)

    def CalcNumOfSpectrum(self):
        # Теоретические значения
        if self.BaseSet == 'Hex':
            NumSp = 3*self.Size*(self.Size-1)
            NumSpDeg = NumSp/3
        elif self.BaseSet == 'Square':
            NumSp = self.Size*self.Size
            NumSpDeg = (NumSp-NumSp%2)/4
            NumSp -= 1     
        else:
            NumSp = Size*Size
            NumSpDeg = int(NumSp/Symm)
        
        self.numSp, self.numSp2 = int(NumSp), int(NumSpDeg) 
        self.numZ = self.numSp - 2*self.numSp2 # Количество невырожденных спектров

    def BaseSets(self):
        return ['Hex', 'Square', '6-star', '5-star', '4-star', '6-border', '4-border', 'Circle']
        
    def SetPoints(self):
        '''Набор базовых точек'''
        if self.BaseSet == 'Hex': # Гексагональная сетка
            self.vSet = Geo.vSetHex(size=[self.Size, self.Size])
        
        elif self.BaseSet == 'Square':
            self.vSet = Geo.mSet2(self.Size, self.Size) # Квадратная сетка

        elif self.BaseSet == '6-star':
            self.vSet = Geo.mSetRegular(nPoints=6, mult=self.Size, div=self.Size-1, withZero=True) #, accur=10

        elif self.BaseSet == '5-star':
            self.vSet = Geo.mSetRegular(nPoints=5, mult=self.Size, div=self.Size-1, withZero=True) #, accur=10

        elif self.BaseSet == '4-star':
            self.vSet = Geo.mSetRegular(nPoints=4, mult=self.Size, div=self.Size-1, withZero=True) #, accur=10

        elif self.BaseSet == '6-border':
            self.vSet = Geo.mSetRegular(nPoints=6, r = self.Size, div=self.Size-2, withZero=False) #, accur=10

        elif self.BaseSet == '4-border':
            self.vSet = Geo.mSetRegular(nPoints=4, r = self.Size, div=self.Size-2, withZero=False) #, accur=10

        elif self.BaseSet == 'Circle':
            self.vSet = Geo.mSetRegular(nPoints=self.Size, r = self.Size)

        else:
            self.vSet = []
            #self.vSet = Geo.mSetRegular(nPoints=self.Symm, mult=self.Size, withZero=True, accur=7)
            #self.vSet = Geo.vSetTri(size=[self.Size, self.Size])

    def SetSpectr(self):
        self.SetPoints()
        self.mR2 = Geo.mSetToD2(self.vSet) # Получили матрицу квадратов расстояний
        #self.CalcNumOfSpectrum()
        self.SetSpectrData() 

    def __init__(self, BaseSet='Hex', Size=4, distortion=-0.04, degree=1, error=0.000001):
        self.BaseSet = BaseSet
        self.Size = Size 
        self.distortion = distortion
        self.degree = degree
        self.error = error
        self.SetSpectr()

class ShowSpectr(DataSpectr):
    '''Визуализация спектров'''
    def PlotTypes():
        return ['Points', 'Web', 'Mosaic', 'Contour']
        
    def SetPlots(self):
        self.fig.clear()
        self.ax = []
        nr, nc = 1, 1
        if self.NumPlots == 4: nr, nc = 2, 2 
        elif self.NumPlots == 9: nr, nc = 3, 3
        for i in range(self.NumPlots):
            ax = self.fig.add_subplot(nr, nc, i+1, frame_on=False) # aspect='equal' 2, 2, i+1 nrows=nr, ncols=nc, plot_number=i+1
            ax.set_aspect('equal') #, 'datalim'
            ax.axison = False
            self.ax.append(ax) 
        #self.ax.set_projection('lambert')
    
    def IniCanvas(self, ini=True):
        '''Инициализация холста'''
        if ini:
            self.fig = plt.figure(facecolor='white') #
            plt.subplots_adjust(left=0.0, right=1.0, bottom=0., top=0.97, wspace = 0.05, hspace = 0.05)
        self.SetPlots()
        
    def getMask(self, fmask, bmask=True):
        # Mask to hide triangles
        maxlim = 100
        lim = Geo.Val2Val(maxlim - self.iMask, [0, maxlim], [min(fmask), max(fmask)])
        return NP.where(fmask > lim, not bmask, bmask)

    def UpdateSpectr(self, ax, vXs, vYs, vSize, vColor, title=''):
        ax.autoscale_view(tight=None, scalex=True, scaley=True)
        cmap = self.cmap
        if self.inverseColor:
            cmap = cmap + "_r"

        vX, vY, cL = Geo.RotateData(vXs, vYs)
        
        if self.PlotType == 'Points':
            ax.scatter(vX, vY, s=vSize, c=vColor, marker=self.mark_style, cmap=cmap, alpha=self.alpha) #, norm=None, vmin=None, vmax=None, linewidths=None, verts=None, edgecolors=None
        else: # triangulation
            triang = tri.Triangulation(vX, vY)
            if self.useMask:
                xm = Agregat(vX[triang.triangles], self.MaskFunction, 1)
                ym = Agregat(vY[triang.triangles], self.MaskFunction, 1)
                triang.set_mask(self.getMask(xm*xm + ym*ym, self.inverseMask))

            if self.PlotType == 'Web':
                ax.scatter(vX, vY, s=vSize, c=vColor, marker=self.mark_style, cmap=cmap, alpha=self.alpha) #, norm=None, vmin=None, vmax=None, linewidths=None, verts=None, edgecolors=None
                ax.triplot(triang, marker='', linestyle='-', alpha=self.alpha) # можно выводить узлы одновременно с сеткой
                #ax.triplot(triang, marker=self.mark_style, ms=2**(self.msize), markerfacecolor='blue', linestyle='-', alpha=self.alpha) # color=vColor, cmap=self.cmap
            
            elif self.PlotType == 'Mosaic': 
                ax.tripcolor(triang, self.vColor, edgecolors='none', cmap=cmap, alpha=self.alpha) #shading='flat'

            elif self.PlotType == 'Contour':
                if self.useVectorColor: 
                    ax.tricontour(triang, self.vColor, cmap=cmap, alpha=self.alpha)
                else:
                    ax.tricontour(triang, self.vColor, colors='blue', alpha=self.alpha) #'k'

        if self.ShowTitle:
            ax.set_title(title, fontsize=11)
        
        ax.figure.canvas.draw()

    def ShowSp(self):
        vSize = 100*(2**(self.msize-4))    
        if self.useVectorMarksize:
            if self.inverseMarksize:
                vSize /= (700*self.vSize)
            else:
                vSize *= (self.vSize)
        
        vColor = None
        if self.useVectorColor:
            vColor = self.vColor
        
        if self.iSpectr >= self.numSp2 - self.NumPlots:
            self.iSpectr = self.numSp2 - self.NumPlots 
        if self.iSpectr < 0: self.iSpectr = 0

        for i in range(self.NumPlots):
            cIndex = self.iSpectr + i
            ax = self.ax[i]
            ax.clear()
            ax.axison = False
            if cIndex >= len(self.vIndex): continue
            spInd = self.vIndex[cIndex]
            if self.vS[cIndex] < 0: spInd += 1 # Поскольку нулевой уровень исключаем, то отрицательные индексы сдвигаются относительно положительных 
            title = str(spInd) + ':' + str(round(self.vS[cIndex], 3))
            self.UpdateSpectr(ax, self.vDataX[cIndex], self.vDataY[cIndex], vSize, vColor, title) 

    def ChangeIndex(self, val=0):
        self.iSpectr = val
        self.ShowSp()

    def SetVectorMS(self, index=-1):
        maxIndex = self.numZ-1
        if maxIndex < 0: return 
        if index >= 0: self.iMarksize = index
        if self.iMarksize < 0: self.iMarksize = 0
        elif self.iMarksize > maxIndex: self.iMarksize = maxIndex  
        
        vS = self.vDataZ[self.iMarksize]
        maxVal = NP.max(vS) #, minVal = NP.min(vS)
        self.vSize = (vS*vS)/(maxVal*maxVal) + 0.001

    def SetMarksize(self, val=0):
        self.SetVectorMS(int(float(val)))
        self.ShowSp()

    def ChangeMarksize(self, val=0):
        self.msize = val
        self.ShowSp()

    def SetColor(self, index=0):
        maxIndex = self.numZ-1
        if maxIndex < 0: return
        if index >= 0: self.iColor = index
        if self.iColor < 0: self.iColor = 0
        elif self.iColor > maxIndex: self.iColor = maxIndex
        self.vColor = self.vDataZ[self.iColor]

    def SetAlpha(self, val=0):
        self.alpha = val
        self.ShowSp()
        
    def SetData(self, show=False):
        self.SetSpectrData()
        if show: self.ShowSp()
        
    def ChangeDistortion(self, val=0):
        self.distortion = val
        self.SetData(True) 

    def ChangeDegree(self, val=0):
        self.degree = val
        self.SetData(True)

    def SetForm(self):
        self.SetSpectr()
        self.SetColor(-1)
        self.SetVectorMS(-1)
        self.ShowSp()
            
    def dictDefaultValues(self):
        return {'BaseSet': 'Hex',  #'Symm': 6, 
                'Size': 10, 'SizeRange': [2, 30], 'error': 0.000002,
                'distortion': -2., 'distortionRange': [-2., 2.], 'degree': 1, 'degreeRange': [0.01, 2.01],
                'iSpectr': 0, 'iColor': 0, 'iMarksize': 0,
                'NumPlots': 1, 'PlotType': 'Points', 'ShowTitle': False,
                'useMask': False, 'iMask': 0, 'inverseMask': False, 'MaskFunction': 'Std',
                'cmap': 'hsv', 'inverseColor': False, 'useVectorColor': False, 'alpha': 0.5,
                'msize': 2, 'mark_style': 'o', 'useVectorMarksize': False, 'inverseMarksize': False, 'MarksizeRange': [1, 10],
                }

    def Spectr2Dict(self):
        dictSp = self.dictDefaultValues() 
        for key in dictSp.keys():
            try: dictSp[key] = getattr(self, key) #dictSp[key] = self.__getattr__(key)
            except: print('Spectr2Dict(): invalid attribute name: ' + key)
        return dictSp
        
    def LoadFromDict(self, Dict=None, ini=False):
        dictSp = self.dictDefaultValues()
        if Dict is not None:
            dictSp.update(Dict)
        
        DataSpectr.__init__(self, BaseSet=dictSp['BaseSet'], Size=dictSp['Size'], distortion=dictSp['distortion'], degree=dictSp['degree'], error=dictSp['error'])
        for key, value in dictSp.items():
            try: self.__setattr__(key, value)
            except: print('LoadFromDict(): invalid attribute name: ' + key)

        self.IniCanvas(ini)
        self.SetColor(-1)
        self.SetVectorMS(-1)
        self.ShowSp()
            
    def __init__(self, dictAttr = None):
        self.LoadFromDict(dictAttr, ini=True)

