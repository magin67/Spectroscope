# -*- coding:utf-8 -*-
#!/usr/bin/env python
""" Визуализация спектров. Управление
Created on 20.11.2016, @author: dmagin
"""
import numpy as NP

import SpectrData as SD
import Geometry as Geo

import tkinter as tk
from tkinter import ttk
import tkinter.filedialog as filedialog 
 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2TkAgg

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
        self.rbSymm6 = ttk.Radiobutton(self.frControl, text="Hex", variable=Symm, value=6, command=ChangeSymm)
        self.rbSymm4 = ttk.Radiobutton(self.frControl, text="Square", variable=Symm, value=4, command=ChangeSymm)
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
        ## Количество рисунков (1, 4, 9)
        def SetNumPlots():
            self.shSpectr.NumPlots = NumPlots.get()
            self.scIndex.configure(to=int(self.shSpectr.numSp2-self.shSpectr.NumPlots)) #, value=0
            self.shSpectr.SetPlots()
            self.shSpectr.ShowSp()
            
        self.AddLabel(rowControl, "Num of plots")
    
        NumPlots = tk.IntVar(self.frControl)
        NumPlots.set(self.shSpectr.NumPlots)
        self.NumPlots1 = ttk.Radiobutton(self.frControl, text="1", variable=NumPlots, value=1, command=SetNumPlots)
        self.NumPlots2 = ttk.Radiobutton(self.frControl, text="4", variable=NumPlots, value=4, command=SetNumPlots)
        self.NumPlots4 = ttk.Radiobutton(self.frControl, text="9", variable=NumPlots, value=9, command=SetNumPlots)
    
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

    def SetMaskWidgets(self, rowControl):
        ## Маска триангуляции
        def ChangeMaskFunction(val='Var'):
            self.shSpectr.MaskFunction = val
            self.shSpectr.ShowSp()

        def SetMask(val=0, ini=False):
            self.shSpectr.iMask = int(float(val))
            self.chkMask.configure(text='Use: ' + str(int(self.shSpectr.iMask)))
            if not ini:
                self.shSpectr.ShowSp()
            
        def CheckMask():
            self.shSpectr.useMask = useMask.get()
            self.shSpectr.ShowSp()
            
        def InverseMask():
            self.shSpectr.inverseMask = inverseMask.get()
            self.shSpectr.ShowSp()
            
        rowControl = self.AddLabel(rowControl, "Mask", colsp=1)

        vsMaskFunction = ['Std', 'Mean']
        sMaskFunction = tk.StringVar(self.frControl)
        self.omMaskFunction = ttk.OptionMenu(self.frControl, sMaskFunction, "Std", *vsMaskFunction, command=ChangeMaskFunction)  
        rowControl = PlaceWidget(self.omMaskFunction, rowControl, col=1, stick='w') 

        useMask = tk.IntVar()
        useMask.set(False)
        self.chkMask = ttk.Checkbutton(self.frControl, variable=useMask, onvalue=True, offvalue=False, command=CheckMask)
        rowControl = PlaceWidget(self.chkMask, rowControl, col=0, stick='ew') 
    
        self.scMask = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=100, value=self.shSpectr.iMask, command=SetMask)
        rowControl = PlaceWidget(self.scMask, rowControl, col=1, stick='ne')
        SetMask(val=self.shSpectr.iMask, ini=True) 

        inverseMask = tk.BooleanVar()
        inverseMask.set(self.shSpectr.inverseMask)
        self.chkInverseMask = ttk.Checkbutton(self.frControl, text="Inverse", variable=inverseMask, onvalue=True, offvalue=False, command=InverseMask)
        rowControl = PlaceWidget(self.chkInverseMask, rowControl, col=0, colspan=2, stick='wn') 

        return self.AddLabel(rowControl)

    def SetFunctionWidgets(self, rowControl):
        ## Функция
        def SetDistortionLabel():
            self.lblDistortion.configure(text="Disturbance: " + str(round(self.shSpectr.distortion, 3)))
            
        def SetDistortion(val=0):
            distortion = Geo.Val2Val(int(float(val)), [0, 100], [self.shSpectr.distortionRange[0], self.shSpectr.distortionRange[1]])  
            self.shSpectr.ChangeDistortion(distortion)
            SetDistortionLabel()
            
        def SetDegreeLabel():
            self.lblDegree.configure(text="Degree: " + str(round(self.shSpectr.degree, 2)))

        def SetDegree(val=0):
            degree = Geo.Val2Val(int(float(val)), [0, 100], [self.shSpectr.degreeRange[0], self.shSpectr.degreeRange[1]])  
            self.shSpectr.ChangeDegree(degree)
            SetDegreeLabel()

        rowControl = self.AddLabel(rowControl, "Function")
    
        #rowControl = self.AddLabel(rowControl, "Distortion", self.lblColor, 1)

        self.lblDistortion = ttk.Label(self.frControl)
        PlaceWidget(self.lblDistortion, rowControl)
        
        iniDist = Geo.Val2Val(self.shSpectr.distortion, [self.shSpectr.distortionRange[0], self.shSpectr.distortionRange[1]], [0, 100])
        self.scDistortion = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=100, value=iniDist, command=SetDistortion) 
        rowControl = PlaceWidget(self.scDistortion, rowControl, colspan=2, stick='ne')
        SetDistortionLabel() 
    
        rowControl = self.AddLabel(rowControl, "Degree", self.lblColor, 1)

        self.lblDegree = ttk.Label(self.frControl)
        PlaceWidget(self.lblDegree, rowControl)
        
        iniDegree = Geo.Val2Val(self.shSpectr.degree, [self.shSpectr.degreeRange[0], self.shSpectr.degreeRange[1]], [0, 100])
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
            self.lblCmap.configure(text='Map: ' + cmap)
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
    
        rowControl = self.AddLabel(rowControl, "Color", colsp=1)
        
        inverseColor = tk.BooleanVar()
        inverseColor.set(False)
        self.chkInverseColor = ttk.Checkbutton(self.frControl, text="Reverse", variable=inverseColor, onvalue=True, offvalue=False, command=InverseColor)
        rowControl = PlaceWidget(self.chkInverseColor, rowControl, col=1, stick='w') 
        
        cmap = self.shSpectr.cmap
        self.lblCmap = ttk.Label(self.frControl, text='Map: ' + cmap)
        PlaceWidget(self.lblCmap, rowControl)
        vColorMap = SD.mColorMap()
        indCMap = vColorMap.index(cmap) 
        self.scColorMap = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, from_=0, to=len(vColorMap)-1, value=indCMap, command=ChangeColorMap)
        rowControl = PlaceWidget(self.scColorMap, rowControl, col=1, stick='ne') 
    
        useVectorColor = tk.IntVar()
        useVectorColor.set(self.shSpectr.useVectorColor)
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
            self.lblMarksize.configure(text="Size: " + str(round(self.scMarksize.get(), 1)))

        def SetMarksize(val=0):
            self.shSpectr.ChangeMarksize(float(val))
            SetMarksizeLabel()
            
        def SetVmarksize(val=0, ini=False):
            ms = float(val)
            if not ini:
                self.shSpectr.SetMarksize(ms)
            self.chkMarkSize.configure(text="Vector: " + str(int(ms+3)))
    
        def CheckMarksize():
            self.shSpectr.useVectorMarksize = useVectorMarksize.get()
            self.shSpectr.ShowSp()

        def InverseMarksize():
            self.shSpectr.inverseMarksize = inverseMarksize.get()
            self.shSpectr.ShowSp()

        rowControl = self.AddLabel(rowControl, "Markers", colsp=1)
        
        vsMarkStyle = list(SD.dicMarkStyle().keys())
        sMarkStyle = tk.StringVar(self.frControl)
        self.omMarkStyle = ttk.OptionMenu(self.frControl, sMarkStyle, "circle", *vsMarkStyle, command=ChangeMarkform)  
        rowControl = PlaceWidget(self.omMarkStyle, rowControl, col=1, stick='wn') 

        self.lblMarksize = ttk.Label(self.frControl)
        PlaceWidget(self.lblMarksize, rowControl)
        
        self.scMarksize = ttk.Scale(self.frControl, orient='horizontal', length=self.colWidth2, 
                from_=self.shSpectr.MarksizeRange[0], to=self.shSpectr.MarksizeRange[1], value=self.shSpectr.msize, command=SetMarksize)
        rowControl = PlaceWidget(self.scMarksize, rowControl, col=1, stick='ne')
        SetMarksizeLabel() 
    
        useVectorMarksize = tk.BooleanVar()
        useVectorMarksize.set(self.shSpectr.useVectorMarksize)
        self.chkMarkSize = ttk.Checkbutton(self.frControl, variable=useVectorMarksize, onvalue=True, offvalue=False, command=CheckMarksize)
        rowControl = PlaceWidget(self.chkMarkSize, rowControl, col=0, stick='ew') 
        
        self.scVMarksize = ttk.Scale(self.frControl, orient='horizontal', length=100, from_=0, to=self.shSpectr.numZ-1, value=self.shSpectr.iMarksize, command=SetVmarksize)
        rowControl = PlaceWidget(self.scVMarksize, rowControl, col=1, stick='ew')
        SetVmarksize(val=self.shSpectr.iMarksize, ini=True)

        inverseMarksize = tk.BooleanVar()
        inverseMarksize.set(self.shSpectr.inverseMarksize)
        self.chkInverseMarksize = ttk.Checkbutton(self.frControl, text="Reverse", variable=inverseMarksize, onvalue=True, offvalue=False, command=InverseMarksize)
        rowControl = PlaceWidget(self.chkInverseMarksize, rowControl, col=0, stick='wn') 
        
        return rowControl

    def UpdateControlWidgets(self):
        # Инициализация панели управления
        self.lblColor, self.colWidth2 = 'black', 150
        try: self.frControl.destroy()
        except: a=1
        self.frControl = ttk.Frame(self.root, height=2, borderwidth=10, relief='sunken') # relief: flat, groove, raised, ridge, solid, sunken
        self.frControl.grid(row=0, column=1, sticky='nse')
        
        rowControl = self.SetBaseWidgets(0)
        rowControl = self.SetNumPlotsWidgets(rowControl)
        rowControl = self.SetPlotWidgets(rowControl)
        rowControl = self.SetFunctionWidgets(rowControl)
        rowControl = self.SetColorWidgets(rowControl)
        rowControl = self.SetMarkerWidgets(rowControl)
        rowControl = self.SetMaskWidgets(rowControl)
        
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
            root.destroy()  # this is necessary on Windows
        
        menuSpectr = tk.Menu(self.root)
        self.root.config(menu=menuSpectr)
 
        fm = tk.Menu(menuSpectr, tearoff=0) # стандартные операции
        fm.add_command(label="New", command=NewSpectr)
        fm.add_command(label="Open...", command=OpenFromDict)
        fm.add_command(label="Save...", command=SaveDict)
        fm.add_command(label="Save as picture...", command=Image2File)
        fm.add_separator()
        fm.add_command(label="Exit", command=_quit)
        menuSpectr.add_cascade(label="Spectr", menu=fm)
         
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
        self.shSpectr.ShowSp()

def main(title):
    # Основное тело
    root = tk.Tk()
    root.wm_title(title)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    shSpectr = SD.ShowSpectr()
    ctrWidgets = ControlWidget(root, shSpectr)

    tk.mainloop()

main("Spectroscope 2D (v0.33)")

# TODO:
# * Управление релевантностью элементов управления
# * Сохранение/восстановление контекста
# * Добавить базы - 3, 5, 7, 8. 
# * Использовать произвольную конфигурацию точек на плоскости, текущий спектр как базу.
# * Анимация с сохранением в гифку. Генерация случайного. 
# * Сдвиг палитры. Если анимировать, то мозаики начнут играть цветом.
