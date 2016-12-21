# -*- coding:utf-8 -*-
#!/usr/bin/env python
""" Визуализация спектров решеток.
Начало: ноябрь, 2016
Редактирован: 16.12.2016
"""
import numpy as NP
from numpy import linalg as LA

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt, mpld3
import matplotlib.tri as tri

from matplotlib.widgets import Slider, Button #, RadioButtons
import matplotlib.markers as Markers

import tkinter as tk
from tkinter import ttk
import tkinter.filedialog as filedialog 
 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import Geometry as Geo

def Val2Val(val, rangeFrom = [0, 100], rangeTo=[0, 1]):
    # Приведение числа из одного диапазона к другому - Normalisation
    # val - число из диапазона rangeFrom.
    # rangeTo - диапазон, к которому надо привести
    relval = (val-rangeFrom[0])/(rangeFrom[1]-rangeFrom[0])
    return rangeTo[0] + relval*(rangeTo[1] - rangeTo[0])  

def mR2G(mR):
    # Преобразование матрицы резистенсов в матрицу Грина.
    # Или, что тоже самое, - матрицы квадратов расстояний в матрицу корреляции
    # G = -(eL R eL)/2
    N = len(mR)
    eL = NP.eye(N) - NP.ones(N)/N
    return NP.dot(NP.dot(eL, mR), eL)/(-2.)

class DataSpectr(object):
    '''Класс данных спектра'''
    
    def Spectr2Projections(self, sVal, mVector):
        '''Получаем вырожденные 2D-проекции спектра'''
        self.vIndex = [] # индексы вырожденных спектров
        self.vS = [] # собственные числа вырожденных спектров
        self.vDataX, self.vDataY = [], [] # собственные функции вырожденных спектров
        self.vDataZ = [] # собственные функции невырожденных спектров
         
        curS = sVal[0],
        NumDims = len(sVal)
        for i in range(1, NumDims):
            sV = sVal[i]
            if abs(sV) < 2*self.error: continue # исключаем шум около ноля
            if abs(sV - curS) < self.error:
                self.vIndex.append(i-1)
                self.vS.append(sV)
                vX, vY, cL = Geo.RotateData(mVector[i], mVector[i-1])
                self.vDataX.append(vX)
                self.vDataY.append(vY)
                i += 1
            else:
                self.vDataZ.append(mVector[i-1]) #-1
            curS = sV
    
    def SetFunction(self):
        '''Добавление возмущения к исходным данным
        mRD = (1 - distortion)*mR2 + distortion * (mR2**degree) # В таком виде ругается на деление на ноль при отрицательных степенях из-за нулевой диагонали.
        for i in range(len(mR2)):
            mRD[i][i] = 0 
        '''
        size = len(self.mR2)
        #szcoeff = NP.sqrt(size)
        coeff = self.distortion #szcoeff*
        mRD = NP.zeros((size, size))
        lr = range(size)
        for i in lr:
            for j in range(i+1, size):
                r2 = self.mR2[i][j]
                Rd = r2**self.degree
                mRD[i][j] = coeff*Rd + 1/Rd #NP.sqrt(1/r2) #distortion/degv
                mRD[j][i] = mRD[i][j] 
        self.mFunction = mRD 

    def SetSpectrData(self):
        self.SetFunction()
        mG = mR2G(self.mFunction)
        sG, vG = LA.eigh(mG)
        vG = NP.transpose(vG)
        self.Spectr2Projections(sG, vG)

    def CalcNumOfSpectrum(self):
        if self.Symm == 4:
            NumSp = self.Size*self.Size
            NumSpDeg = (NumSp-NumSp%2)/4
            NumSp -= 1     
        elif self.Symm == 6:
            NumSp = 3*self.Size*(self.Size-1)
            NumSpDeg = NumSp/3
        else:
            NumSp = Size*Size
            NumSpDeg = int(NumSp/Symm)
        
        self.numSp, self.numSp2 = int(NumSp), int(NumSpDeg) 
        self.numZ = self.numSp - 2*self.numSp2 # Количество невырожденных спектров
        
    def SetPoints(self):
        '''Набор базовых точек'''
        #if self.Symm == 3:
            #self.vSet = Geo.vSetTri(size=[self.Size, self.Size])
        if self.Symm == 4:
            self.vSet = Geo.mSet2(self.Size, self.Size) # Квадратная сетка
        elif self.Symm == 6: # Гексагональная сетка
            self.vSet = Geo.vSetHex(size=[self.Size, self.Size])
        else:
            self.vSet = Geo.mSetRegular(nPoints=self.Symm, mult=self.Size, withZero=True, accur=7)
        #print(self.vSet)

    def SetSpectr(self):
        self.SetPoints()
        self.mR2 = Geo.mSetToD2(self.vSet) # Получили матрицу квадратов расстояний
        self.CalcNumOfSpectrum()
        self.SetSpectrData() 

    def __init__(self, Symm=4, Size=4, distortion=-0.04, degree=1, error=0.00001):
        self.Symm = Symm
        self.Size = Size 
        self.distortion = distortion
        self.degree = degree
        self.error = error
        self.SetSpectr()

class ShowSpectr(DataSpectr):
    '''Класс визуализации спектра'''
    def SetPlots(self):
        self.fig.clear()
        self.ax = []
        nr, nc = 1, 1
        if self.NumPlots > 1: nc = 2 
        if self.NumPlots > 2: nr = 2
        for i in range(self.NumPlots):
            ax = self.fig.add_subplot(nr, nc, i+1, frame_on=False) # aspect='equal' 2, 2, i+1 nrows=nr, ncols=nc, plot_number=i+1
            ax.set_aspect('equal', 'datalim')
            ax.axison = False
            self.ax.append(ax) 
    
    def IniCanvas(self):
        '''Инициализация холста'''
        self.fig = plt.figure(facecolor='white') #
        self.SetPlots()
        #fig, axes = plt.subplots(ncols=2, nrows=2) #,, sharex=sharex, sharey=sharey, subplot_kw=kw kw={'xticks': [], 'yticks': []} squeeze=True,   
        #self.ax.set_projection('lambert')
        #self.ax.set_autoscale_on(True)
        #self.ax.autoscale_view()
        #self.ax.autoscale_view(tight=None, scalex=True, scaley=True)
        plt.subplots_adjust(left=0.0, right=1.0, bottom=0, top=1.0) #, wspace = 0.1, hspace = 0.1

    def UpdateSpectr(self, ax, vX, vY, vSize, vColor):
        #self.fig.clear()
        ax.clear()
        ax.axison = False
        cmap = self.cmap
        if self.inverseColor:
            cmap = cmap + "_r"
        
        if self.PlotType == 'Points':
            ax.scatter(vX, vY, s=vSize, c=vColor, marker=self.mark_style, cmap=cmap, alpha=self.alpha) #, norm=None, vmin=None, vmax=None, linewidths=None, verts=None, edgecolors=None
        else: # triangulation
            triang = tri.Triangulation(vX, vY)
            if self.useVectorMarksize: # Mask off unwanted triangles.
                xmid = vX[triang.triangles].mean(axis=1)
                ymid = vY[triang.triangles].mean(axis=1)
                lim = ((self.iMarksize+5)/len(self.vDataZ)/4)**2 #self.vSize*self.vSize *(self.iMarksize+1)
                
                if self.inverseMarksize:
                    mask = NP.where(xmid*xmid + ymid*ymid <= lim, 1, 0)
                else:
                    mask = NP.where(xmid*xmid + ymid*ymid > lim, 1, 0)
                triang.set_mask(mask)

            if self.PlotType == 'Web':
                ax.triplot(triang, marker=self.mark_style, ms=2**(self.msize), markerfacecolor='blue', linestyle='-', alpha=self.alpha) # color=vColor, cmap=self.cmap,
            elif self.PlotType == 'Mosaic': # 'Web' and 'Mosaic' can be shown together 
                ax.tripcolor(triang, self.vColor, shading='flat', cmap=cmap, alpha=self.alpha) #shading='gouraud' 
        
        ax.figure.canvas.draw()

    def ShowSp(self):
        vSize = 100*(2**(self.msize))    
        if self.useVectorMarksize:
            if self.inverseMarksize:
                vSize /= (700*self.vSize)
            else:
                vSize *= (self.vSize)
        
        vColor = None
        if self.useVectorColor:
            vColor = self.vColor
        
        iSp = self.iSpectr
        if iSp < 0: return
        
        for i in range(self.NumPlots): 
            self.UpdateSpectr(self.ax[i], self.vDataX[iSp + i], self.vDataY[iSp + i], vSize, vColor) 

    def ChangeIndex(self, val=0):
        self.iSpectr = val
        self.ShowSp()

    def SetVectorMS(self, index=-1):
        maxIndex = len(self.vDataZ)-1
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
        maxIndex = len(self.vDataZ)-1
        if index >= 0: self.iColor = index
        if self.iColor < 0: self.iColor = 0
        elif self.iColor > maxIndex: self.iColor = maxIndex  
        self.vColor = self.vDataZ[self.iColor]

    def indUpdColor(self, val=0):
        self.SetColor(int(float(val)))
        self.ShowSp()

    def SetAlpha(self, val=0):
        self.alpha = val
        self.ShowSp()
        
    def SetData(self, show=False):
        self.SetSpectrData()
        if show: self.ShowSp()
        
    def SetForm(self, iniData=True):
        if iniData:
            self.SetSpectr()
        self.SetColor(-1)
        self.SetVectorMS(-1)
        self.ShowSp()
            
    def ChangeDistortion(self, val=0):
        self.distortion = val
        self.SetData(True) 

    def ChangeDegree(self, val=0):
        self.degree = val
        self.SetData(True)

    def dictDefaultValues(self):
        return {'Symm': 6, 'Size': 10, 'distortion': -0.04, 'degree': 1, 'error': 0.00001,
                'distortionRange': [-0.06, 0.06], 'degreeRange': [0.01, 2.01],
                'iSpectr': 0, 'iColor': 0, 'iMarksize': 0,
                'NumPlots': 1, 'PlotType': 'Points',
                'cmap': 'hsv', 'inverseColor': False, 'useVectorColor': True, 'alpha': 0.5,
                'msize': 2, 'mark_style': 'o', 'useVectorMarksize': False, 'inverseMarksize': False
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
        
        DataSpectr.__init__(self, Symm=dictSp['Symm'], Size=dictSp['Size'], distortion=dictSp['distortion'], degree=dictSp['degree'], error=dictSp['error'])
        for key, value in dictSp.items():
            try: self.__setattr__(key, value)
            except: print('LoadFromDict(): invalid attribute name: ' + key)

        if ini: 
            self.IniCanvas()
        else:
            self.SetPlots()
        
        self.SetForm(False)
            
    def __init__(self, dictAttr = None):
        self.LoadFromDict(dictAttr, ini=True)

###

def dicMarkStyle():
    # инвертированный словарь стилей маркера
    # http://matplotlib.org/api/markers_api.html#module-matplotlib.markers
    return {v:k for k, v in Markers.MarkerStyle().markers.items()}

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

def ValidCMap(smap):
    vColorMap = mColorMap()
    try:
        ind = vColorMap.index(smap)
    except:
        return '' # неправильный идентификатор цветовой карты
    if smap[0:2] == "==":
        return vColorMap[ind+1]
    return smap

###

def PlaceWidget(wd, row, col=0, colspan=1, stick='wn'):
    wd.grid(row=row, column=col, columnspan=colspan, sticky=stick)
    return row + col + (colspan - 1)

class ControlWidget(object):
    ''' Размещение виджетов '''
    def AddLabel(self, row, text = "", fg = "blue", colsp=2):
        frSpectr = ttk.Label(self.frControl, text=text, foreground=fg)
        return PlaceWidget(frSpectr, row, colspan=colsp)
    
    def SetBaseWidgets(self, rowControl):
        ## База (форма, размер, индекс) спектра
        def SetSizeLabel():
            self.lblSpSize.configure(text="Order: " + str(self.shSpectr.Size) + ', ' + str(self.shSpectr.numSp + 1))

        def SetIndexLabel():
            iSpectr = self.shSpectr.iSpectr
            self.lblSpIndex.configure(text="Index: " + str(iSpectr+1) + '/' + str(self.shSpectr.numSp2))

        def SetSpParameters():
            self.scIndex.configure(to=int(self.shSpectr.numSp2-self.shSpectr.NumPlots)) #, value=0
            self.scVMarksize.configure(to=int(self.shSpectr.numZ-1)) #, value=0
            self.scColor.configure(to=int(self.shSpectr.numZ-1)) #, value=0
            SetSizeLabel()
            SetIndexLabel()
        
        def ChangeSymm():
            self.shSpectr.Symm=Symm.get()
            self.shSpectr.SetForm()
            SetSpParameters()

        def ChangeSize(val):
            iSize = int(float(val))
            self.shSpectr.Size= iSize
            self.shSpectr.SetForm()
            SetSpParameters()
                        
        def ChangeIndex(val):
            ind = int(float(val))
            self.shSpectr.ChangeIndex(ind)
            SetIndexLabel()
    
        rowControl = self.AddLabel(rowControl, "Base")
        
        Symm = tk.IntVar()
        Symm.set(self.shSpectr.Symm)
        self.rbSymm6 = ttk.Radiobutton(self.frControl, text="Hex", variable=Symm, value=6, command=ChangeSymm) #.pack(anchor=W)
        self.rbSymm4 = ttk.Radiobutton(self.frControl, text="Square", variable=Symm, value=4, command=ChangeSymm) #.pack(anchor=W)
        rowControl = PlaceWidget(self.rbSymm6, rowControl, col=0, stick='ns')
        rowControl = PlaceWidget(self.rbSymm4, rowControl, col=1, stick='ns')
    
        self.lblSpSize = ttk.Label(self.frControl)
        PlaceWidget(self.lblSpSize, rowControl)
        SetSizeLabel()
        
        self.scSize = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=2, to=20, value=self.shSpectr.Size, command=ChangeSize)
        rowControl = PlaceWidget(self.scSize, rowControl, col=1, colspan=1, stick='ne')
        
        rowControl = self.AddLabel(rowControl, "Index", self.lblColor, 1)

        self.lblSpIndex = ttk.Label(self.frControl) #, text="Index: 0"
        PlaceWidget(self.lblSpIndex, rowControl)
        self.scIndex = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=self.shSpectr.numSp2-self.shSpectr.NumPlots, value=self.shSpectr.iSpectr, command=ChangeIndex) #tickinterval=5, , resolution=5
        rowControl = PlaceWidget(self.scIndex, rowControl, col=0, colspan=2, stick='ne')
        SetIndexLabel()
         
        return self.AddLabel(rowControl)

    def SetNumPlotsWidgets(self, rowControl):
        ## Количество рисунков (1, 2, 4)
        def SetNumPlots():
            self.shSpectr.NumPlots = NumPlots.get()
            self.scIndex.configure(to=int(self.shSpectr.numSp2-self.shSpectr.NumPlots)) #, value=0
            self.shSpectr.SetPlots()
            self.shSpectr.ShowSp()
            
        self.AddLabel(rowControl, "Num of plots")
    
        NumPlots = tk.IntVar(self.frControl)
        NumPlots.set(self.shSpectr.NumPlots)
        self.NumPlots1 = ttk.Radiobutton(self.frControl, text="1", variable=NumPlots, value=1, command=SetNumPlots)
        self.NumPlots2 = ttk.Radiobutton(self.frControl, text="2", variable=NumPlots, value=2, command=SetNumPlots)
        self.NumPlots4 = ttk.Radiobutton(self.frControl, text="4", variable=NumPlots, value=4, command=SetNumPlots)
    
        PlaceWidget(self.NumPlots1, rowControl, col=1, stick='w')
        PlaceWidget(self.NumPlots2, rowControl, col=1, stick='n')
        rowControl = PlaceWidget(self.NumPlots4, rowControl, col=1, stick='e')
    
        return self.AddLabel(rowControl)

    def SetPlotWidgets(self, rowControl):
        ## Тип рисунка
        def SetPlotType():
            self.shSpectr.PlotType = PlotType.get()
            self.shSpectr.ShowSp()
            
        rowControl = self.AddLabel(rowControl, "Plot type")
    
        PlotType = tk.StringVar(self.frControl)
        PlotType.set(self.shSpectr.PlotType)
        self.rbPoints = ttk.Radiobutton(self.frControl, text="Points", variable=PlotType, value='Points', command=SetPlotType)
        self.rbMosaic = ttk.Radiobutton(self.frControl, text="Mosaic", variable=PlotType, value='Mosaic', command=SetPlotType)
        self.rbWeb = ttk.Radiobutton(self.frControl, text="Web", variable=PlotType, value='Web', command=SetPlotType)
    
        PlaceWidget(self.rbPoints, rowControl, col=0, stick='w')
        PlaceWidget(self.rbMosaic, rowControl, col=1, stick='w')
        rowControl = PlaceWidget(self.rbWeb, rowControl, col=1, stick='e')
    
        return self.AddLabel(rowControl)

    def SetFunctionWidgets(self, rowControl):
        ## Функция
        def SetDistortionLabel():
            self.lblDistortion.configure(text="Distur: " + str(round(self.shSpectr.distortion, 3)))
            
        def SetDistortion(val=0):
            distortion = Val2Val(int(float(val)), [0, 100], [self.shSpectr.distortionRange[0], self.shSpectr.distortionRange[1]])  
            self.shSpectr.ChangeDistortion(distortion)
            SetDistortionLabel()
            
        def SetDegreeLabel():
            self.lblDegree.configure(text="Degree: " + str(round(self.shSpectr.degree, 2)))

        def SetDegree(val=0):
            degree = Val2Val(int(float(val)), [0, 100], [self.shSpectr.degreeRange[0], self.shSpectr.degreeRange[1]])  
            self.shSpectr.ChangeDegree(degree)
            SetDegreeLabel()

        rowControl = self.AddLabel(rowControl, "Function")
    
        #rowControl = self.AddLabel(rowControl, "Distortion", self.lblColor, 1)

        self.lblDistortion = ttk.Label(self.frControl)
        PlaceWidget(self.lblDistortion, rowControl)
        
        iniDist = Val2Val(self.shSpectr.distortion, [self.shSpectr.distortionRange[0], self.shSpectr.distortionRange[1]], [0, 100])
        self.scDistortion = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=100, value=iniDist, command=SetDistortion) 
        rowControl = PlaceWidget(self.scDistortion, rowControl, colspan=2, stick='ne')
        SetDistortionLabel() 
    
        rowControl = self.AddLabel(rowControl, "Degree", self.lblColor, 1)

        self.lblDegree = ttk.Label(self.frControl)
        PlaceWidget(self.lblDegree, rowControl)
        
        iniDegree = Val2Val(self.shSpectr.degree, [self.shSpectr.degreeRange[0], self.shSpectr.degreeRange[1]], [0, 100])
        self.scDegree = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=100, value=iniDegree, command=SetDegree) 
        rowControl = PlaceWidget(self.scDegree, rowControl, colspan=2, stick='ne')
        SetDegreeLabel() 
    
        return self.AddLabel(rowControl)

    def SetColorWidgets(self, rowControl):
        ## Цвет
        def CheckVectorColor(): 
            self.shSpectr.useVectorColor = useVectorColor.get()
            self.shSpectr.ShowSp()
        
        def ChangeColorMap(val=0):
            cmap = mColorMap()[int(float(val))]
            self.shSpectr.cmap = cmap
            self.lblCmap.configure(text=cmap)
            self.shSpectr.ShowSp()
            
        def SetVColor(val=0, ini=False):
            iColor = int(float(val))
            if not ini:
                self.shSpectr.indUpdColor(iColor)
            self.chkVectorColor.configure(text="Vector: " + str(iColor + 1))
            
        def InverseColor():
            self.shSpectr.inverseColor = inverseColor.get()
            self.shSpectr.ShowSp()
            
        def SetAlpha(val=0, ini=False):
            if not ini:
                self.shSpectr.SetAlpha(float(val))
            self.lblAlpha.configure(text="Alpha: " + str(round(float(val), 2)))
    
        rowControl = self.AddLabel(rowControl, "Color map", colsp=1)
        
        inverseColor = tk.BooleanVar()
        inverseColor.set(False)
        self.chkInverseColor = ttk.Checkbutton(self.frControl, text="Reverse", variable=inverseColor, onvalue=True, offvalue=False, command=InverseColor)
        rowControl = PlaceWidget(self.chkInverseColor, rowControl, col=1, stick='w') 
        
        cmap = self.shSpectr.cmap
        self.lblCmap = ttk.Label(self.frControl, text=cmap)
        PlaceWidget(self.lblCmap, rowControl)
        vColorMap = mColorMap()
        indCMap = vColorMap.index(cmap) 
        self.scColorMap = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=len(vColorMap)-1, value=indCMap, command=ChangeColorMap)
        rowControl = PlaceWidget(self.scColorMap, rowControl, col=1, stick='ne') 
    
        useVectorColor = tk.IntVar()
        useVectorColor.set(True)
        self.chkVectorColor = ttk.Checkbutton(self.frControl, variable=useVectorColor, onvalue=True, offvalue=False, command=CheckVectorColor)
        rowControl = PlaceWidget(self.chkVectorColor, rowControl, col=0, stick='ew') 
        
        self.scColor = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=self.shSpectr.numZ-1, value=self.shSpectr.iColor, command=SetVColor)
        rowControl = PlaceWidget(self.scColor, rowControl, col=1, stick='ne')
        SetVColor(val=self.shSpectr.iColor, ini=True) 
    
        self.lblAlpha = ttk.Label(self.frControl)
        PlaceWidget(self.lblAlpha, rowControl)
        #rowControl = self.AddLabel(rowControl, "Saturation", self.lblColor, 1)
        self.scAlpha = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0.1, to=1, value=self.shSpectr.alpha, command=SetAlpha)
        rowControl = PlaceWidget(self.scAlpha, rowControl, col=1, stick='ne')
        SetAlpha(self.scAlpha.get(), True) 
    
        return self.AddLabel(rowControl)

    def SetMarkerWidgets(self, rowControl):
        ## Markers
        def ChangeMarkform(val='o'):
            self.shSpectr.mark_style = dicMarkStyle()[val]
            self.shSpectr.ShowSp()

        def SetMarksizeLabel():
            self.lblMarksize.configure(text="Size: " + str(round(self.scMarksize.get()+2, 1)))

        def SetMarksize(val=0):
            self.shSpectr.ChangeMarksize(float(val))
            SetMarksizeLabel()
            
        def SetVmarksize(val=0, ini=False):
            ms = float(val)
            if not ini:
                self.shSpectr.SetMarksize(ms)
            self.chkMarkSize.configure(text="Vector: " + str(int(ms+1)))
    
        def CheckMarksize():
            self.shSpectr.useVectorMarksize = useVectorMarksize.get()
            self.shSpectr.ShowSp()

        def InverseMarksize():
            self.shSpectr.inverseMarksize = inverseMarksize.get()
            self.shSpectr.ShowSp()

        rowControl = self.AddLabel(rowControl, "Markers", colsp=1)
        
        vsMarkStyle = sorted(list(dicMarkStyle().keys()))
        sMarkStyle = tk.StringVar(self.frControl)
        self.omMarkStyle = ttk.OptionMenu(self.frControl, sMarkStyle, "circle", *vsMarkStyle, command=ChangeMarkform)  
        rowControl = PlaceWidget(self.omMarkStyle, rowControl, col=1, stick='w') 

        self.lblMarksize = ttk.Label(self.frControl)
        PlaceWidget(self.lblMarksize, rowControl)
        
        self.scMarksize = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=-1., to=8., value=self.shSpectr.msize, command=SetMarksize)
        rowControl = PlaceWidget(self.scMarksize, rowControl, col=1, stick='ne')
        SetMarksizeLabel() 
    
        useVectorMarksize = tk.BooleanVar()
        useVectorMarksize.set(False)
        self.chkMarkSize = ttk.Checkbutton(self.frControl, variable=useVectorMarksize, onvalue=True, offvalue=False, command=CheckMarksize)
        rowControl = PlaceWidget(self.chkMarkSize, rowControl, col=0, stick='ew') 
        
        self.scVMarksize = ttk.Scale(self.frControl, orient='horizontal', length=100, from_=0, to=self.shSpectr.numZ-1, value=self.shSpectr.iMarksize, command=SetVmarksize)
        rowControl = PlaceWidget(self.scVMarksize, rowControl, col=1, stick='ew')
        SetVmarksize(val=self.shSpectr.iMarksize, ini=True)

        inverseMarksize = tk.BooleanVar()
        inverseMarksize.set(False)
        self.chkInverseMarksize = ttk.Checkbutton(self.frControl, text="Reverse", variable=inverseMarksize, onvalue=True, offvalue=False, command=InverseMarksize)
        rowControl = PlaceWidget(self.chkInverseMarksize, rowControl, col=0, stick='wn') 
        
        return rowControl

    def UpdateControlWidgets(self):
        # Инициализация панели управления
        self.lblColor, self.colWidth2 = 'black', 140
        try: self.frControl.destroy()
        except: a=1
        self.frControl = ttk.Frame(self.root, height=2, borderwidth=10, relief='sunken') #bg='green', relief: flat, groove, raised, ridge, solid, sunken
        self.frControl.grid(row=0, column=1, sticky='nse')
        
        rowControl = self.SetBaseWidgets(0)
        rowControl = self.SetNumPlotsWidgets(rowControl)
        rowControl = self.SetPlotWidgets(rowControl)
        rowControl = self.SetFunctionWidgets(rowControl)
        rowControl = self.SetColorWidgets(rowControl)
        rowControl = self.SetMarkerWidgets(rowControl)
        
    def SetMenuWidget(self):
        def NewSpectr():
            self.shSpectr.LoadFromDict()
            self.UpdateControlWidgets()
        
        def OpenFromDict():
            options = {}
            options['filetypes'] = [('NPY', '.npy'), ('*', '.*')] #('json', '.json'), 
            options['defaultextension'] = '.npy'
            options['title'] = 'Choose file'
            filelist = filedialog.askopenfiles(**options)
            if filelist is None: return
            file = filelist[0] 
            dictSpectr = NP.load(file.name).item()
            self.shSpectr.LoadFromDict(dictSpectr)
            self.UpdateControlWidgets()
            file.close()

        def SaveDict():
            options = {}
            options['filetypes'] = [('NPY', '.npy')] #, ('json', '.json')
            options['defaultextension'] = '.npy'
            options['initialfile'] = 'mySpectr'
            options['title'] = 'Choose file name'
            file = filedialog.asksaveasfile(**options)
            if file is None: return
            NP.save(file.name, self.shSpectr.Spectr2Dict()) 
            file.close()

        def Image2File():
            options = {}
            options['filetypes'] = [('png', '.png'), ('pdf', '.pdf'), ('svg', '.svg')] #png, pdf, ps, eps, svg  # [('all files', '.*'), ('text files', '.txt')]
            options['defaultextension'] = '.png'
            #options['initialdir'] = 'C:\\'
            options['initialfile'] = 'mySpectr' + str(self.shSpectr.iSpectr) + '.png'
            #options['parent'] = root
            options['title'] = 'Choose file name'
            
            file = filedialog.asksaveasfile(**options)
            #filename = filedialog.asksaveasfilename(**options) ## defaultextension, filetypes, initialdir, initialfile, multiple, message, parent, title
            if file is None: return
            ftype = file.name[(file.name.index('.')+1):]
            fig = self.shSpectr.fig
            fig.savefig(file.name, format=ftype, transparent=True)
            file.close()

        def _quit():
            root.quit()     # stops mainloop
            root.destroy()  # this is necessary on Windows to prevent
        
        menuSpectr = tk.Menu(self.root) #создается объект Меню на главном окне
        self.root.config(menu=menuSpectr) #окно конфигурируется с указанием меню для него
 
        fm = tk.Menu(menuSpectr, tearoff=0) # стандартные операции
        fm.add_command(label="New", command=NewSpectr)
        fm.add_command(label="Open...", command=OpenFromDict)
        fm.add_command(label="Save...", command=SaveDict)
        fm.add_command(label="Save as picture...", command=Image2File)
        fm.add_separator()
        fm.add_command(label="Exit", command=_quit)
        menuSpectr.add_cascade(label="Spectr", menu=fm) #пункт располагается на основном меню
         
        hm = tk.Menu(menuSpectr, tearoff=0) # помощь
        hm.add_command(label="Help")
        hm.add_command(label="About") 
        menuSpectr.add_cascade(label="Help", menu=hm)

    def __init__(self, root, shSpectr):
        self.shSpectr = shSpectr
        self.root = root 
        
        ## Полотно для вывода
        self.canvas = FigureCanvasTkAgg(shSpectr.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky='nswe')
        self.canvas.show()

        self.SetMenuWidget()
            
        self.UpdateControlWidgets()

def main():
    root = tk.Tk()
    root.wm_title("Spectroscope 2D (v0.31)")
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    shSpectr = ShowSpectr()
    
    #print(shSpectr.Spectr2Dict())
    
    ctrWidgets = ControlWidget(root, shSpectr)

    ## Панель формы
    #btnQuit = ttk.Button(root, text='Quit', command=None)
    #btnQuit.grid(row=1, column=1, sticky='es')

    tk.mainloop()

main()


# TODO:
# * Сохранение спектра в файл (параметры генерации), запрос на сохранение при модификации, восстановление при открытии программы
# * Перенос в веб - mplD3
# * Подчистка файлов проекта
 
# * Добавить базы - 3, возможно, также 5 и 7, 8 - быстро не получилось. Можно попробовать использовать текущий спектр как базу
# * Анимация с сохранением в гифку
# * Учет/настройка цвета для режима Web
# * Сдвиг палитры. Если анимировать, то мозаики начнут играть цветом.
# * Управление доступностью элементов управления
# * Регулирование искажений внешней картинкой? - пока непонятно, как правильно искажать

# - done Добавление верхнего меню
# - done Сохранение картинки в файл
# - done Вывод одновременно 2-х или 4-х спектров
# - done Инвертирование палитры
# - done Выбор цветовой карты из дерева - сделано через слайдер

