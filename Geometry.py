# -*- coding:utf-8 -*-
''' Created on 30.06.2015, @author: dmagin'''
import math
import numpy as NP

''' Общее '''

def Val2Val(val, rangeFrom = [0, 100], rangeTo=[0, 1]):
    # Приведение числа из одного диапазона к другому - Normalisation
    # val - число из диапазона rangeFrom.
    # rangeTo - диапазон, к которому надо привести
    relval = (val-rangeFrom[0])/(rangeFrom[1]-rangeFrom[0])
    return rangeTo[0] + relval*(rangeTo[1] - rangeTo[0])  

def RotateData(vDataX, vDataY):
    # Стабилизация данных относительно осей по первой точке
    indX = 0 
    cX, cY = vDataX[indX], vDataY[indX]
    cL = NP.sqrt(cX*cX + cY*cY)
    if cL == 0: return vDataX, vDataY, cL
    if cX < 0:
        cosFi, sinFi = -cX/cL, cY/cL
    else:
        cosFi, sinFi = cX/cL, -cY/cL
    mRotate = [(cosFi, -sinFi), (sinFi, cosFi)]
    mData = [vDataX, vDataY]
    mResult = NP.dot(mRotate, mData)
    
    return mResult[0], mResult[1], cL

def mAddSet(mSet1, mSet2):
    # Объединение двух наборов
    mSet = mSet1.copy()
    mSet3.extend(mSet2)
    return mSet

def r2(point1, point2):
    # Квадрат расстояния между точками
    return sum((x2 - x1)*(x2 - x1) for x1, x2 in zip(point1, point2))

''' Наборы координат узлов '''

def mSetToD2(mSet, coeff = 1.):
    # Формирование матрицы квадратов расстояний на основании координат узлов
    NumNodes = len(mSet)
    rangeIndex = range(NumNodes)
    mR = NP.zeros((NumNodes, NumNodes))
    for i in rangeIndex: 
        for j in rangeIndex:
            mR[i][j] = r2(mSet[i], mSet[j])*coeff 
    return mR

def mSetToD(mSet, coeff = 1., degree = 1.):
    # Формирование матрицы расстояний на основании координат узлов
    # Можно задавать степень - mSetToD(mSet, coeff = 1., degree = 3) вернет матрицу кубов расстояний
    NumNodes = len(mSet)
    rangeIndex = range(NumNodes)
    mR = NP.zeros((NumNodes, NumNodes))
    for i in rangeIndex: 
        for j in rangeIndex:
            mR[i][j] = (r2(mSet[i], mSet[j])**(degree/2))*coeff 
    return mR

def mDivideSet(mSet, N=2, mult=1, accur=8):
    # Мультипликация набора координат точек
    # Делит набор координат mSet на заданную кратность N и формирует новый набор. Повторяется mult раз
    mForDiv = mSet
    for m in range(mult):
        mSetResult = [] 
        for cSet1 in mForDiv:
            if not(cSet1 in mSetResult): mSetResult.append(cSet1)   
            for cSet2 in mForDiv:
                if cSet1 == cSet2: continue
                if not(cSet2 in mSetResult): mSetResult.append(cSet2)
                for i in range(1, N):
                    cSet = []
                    for k in range(len(cSet1)):
                        coord = round(cSet1[k] + i*(cSet2[k] - cSet1[k])/N, accur)
                        cSet.append(coord)
                    cSet = tuple(cSet)
                    if not(cSet in mSetResult): mSetResult.append(cSet)
        mForDiv = mSetResult 
    #print(mSetResult)
    return mSetResult

def mSet(Length):
    # Набор координат узлов (точек) квадратной решетки
    mSet = [] 
    for disp1 in range(Length):
        for disp2 in range(Length):
            mSet.append((disp1, disp2))
    return mSet

def mSet2(NumRow, NumCol = 1, aRow = 1, aCol = 1, dispRow = 0, dispCol = 0):
    # Набор координат узлов двумерной решетки: простая прямоугольная
    mSet = [] 
    for row in range(NumRow):
        for col in range(NumCol):
            mSet.append((col*aCol + dispCol*row, row*aRow + dispRow*col))
    return mSet

def vSet2vect(ev1, ev2, d1 = [-1, 1], d2 = [-1, 1], withLimit = False):
    # Набор координат узлов двумерной решетки: по базисным векторам
    mSet = [] 
    for n1 in range(d1[0], d1[1]+1):
        x1, y1 = ev1[0]*n1, ev1[1]*n1
        for n2 in range(d2[0], d2[1]+1):
            if withLimit and abs(n1-n2) > d1[1]: continue 
            x2, y2 = ev2[0]*n2, ev2[1]*n2
            mSet.append((x1 + x2, y1 + y2))
    return mSet

def vSetHex(r = 1, size = [1, 1]):
    # Гексагональная решетка: сота
    # r - шаг соты, минимальное расстояние между узлами
    # size[] - размеры сетки
    n1, n2 = size[0]-1, size[1]-1
    eSet = mSetRegular(3, r)
    vSet = vSet2vect(eSet[0], eSet[1], [-n1, n1], [-n2, n2], True)
    #print(vSet)
    return vSet

def mSet3(Length1, Length2 = 1, Length3 = 1):
    # Трехмерная прямоугольная решетка
    mSet = [] 
    for x in range(Length1):
        for y in range(Length2):
            for z in range(Length3):
                mSet.append((x, y, z))
    return mSet

def mSetRegular(nPoints = 3, r = 1., cc = [0, 0, 0], mult = 1, div = 0, withZero = False, accur=8):
    # Координаты правильных многоугольников
    #    nPoints - количество вершин,
    #    r - внешний радиус
    #    cc - координаты центра
    #    mult - количество размножений координат в направлении радиуса
    mSet = []
    xc, yc, zc = cc[0], cc[1], cc[2]
    fi = 2*NP.pi/float(nPoints)
    cTo = tuple(cc)
    cr = 0 
    for m in range(mult):
        cr += r
        for i in range(nPoints+1):
            x = round(xc + cr*NP.cos(fi*i), accur)
            y = round(yc + cr*NP.sin(fi*i), accur)
            cFrom = cTo
            cTo = (x, y, zc)
            if i != nPoints: 
                mSet.append(cTo)
            if div > 0 and i > 0: # Разбиваем линию между двумя точками
                #cFrom, cTo = mSet[len(mSet)-2], mSet[len(mSet)-1]
                dx, dy = (cTo[0] - cFrom[0])/(div+1), (cTo[1] - cFrom[1])/(div+1) 
                for j in range(div):
                    xd = round(cFrom[0] + (j+1)*dx, accur)   
                    yd = round(cFrom[1] + (j+1)*dy, accur)   
                    mSet.append((xd, yd, zc))

    if withZero: mSet.append(tuple(cc))
    return mSet

def mSetCross(N = 2, Length = 1., D = 2):
    # Координаты правильного креста
    mSet = []
    a = Length
    for i in range(N):
        mSet.append((i*a, 0, 0))
        if i == 0:
            continue
        mSet.append((-i*a, 0, 0))
        
        if D > 1: # двумерный
            mSet.append((0, i*a, 0))
            mSet.append((0, -i*a, 0))

        if D > 2: # трехмерный
            mSet.append((0, 0, i*a))
            mSet.append((0, 0, -i*a))

    return mSet

def mSetSquareBorder(Size = 1):
    # Набор координат узлов границы квадрата
    mSet = []
    for i in range(Size-1):
        mSet.append((0, i))
        mSet.append((i, Size-1))
        mSet.append((Size-1, Size-1-i))
        mSet.append((Size-1-i, 0))

    return mSet

def mSetAngle(N = 1, a = 1.):
    # Координаты прямого угла
    mSet = []
    for i in range(N):
        mSet.append((i*a, 0))
        if i == 0:
            continue

        mSet.append((0, i*a))
    return mSet

''' Координаты правильных многогранников, r - внешний радиус'''

def mD2Tetrahedron(r = 1.):
    # Тетраэдр - 4
    fi = 1./math.sqrt(2)
    coeff = r/math.sqrt(3./8.)/2
    mCoord = [(-1, 0, -fi), (1, 0, -fi),
            (0, -1, fi), (0, 1, fi)]
    return mSetToD2(mCoord, coeff)

def mD2Octahedron(r = 1.):
    # Октаэдр - 8
    coeff = r
    mCoord = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
    return mSetToD2(mCoord, coeff)

def mD2Cube(r = 1.):
    # Куб - 8
    coeff = r/math.sqrt(3.)
    mCoord = [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1),
              (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1)]
    return mSetToD2(mCoord, coeff)

def mD2Icosahedron(r = 1.):
    # Икосаэдр - 12
    fi = (1.+math.sqrt(5))/2
    coeff = r/math.sqrt(fi*math.sqrt(5))
    mCoord = [(0, -1, -fi), (0, -1, fi),(0, 1, -fi), (0, 1, fi),
            (-1, -fi, 0), (-1, fi, 0),(1, -fi, 0), (1, fi, 0),
            (-fi, 0, -1), (fi, 0, -1),(-fi, 0, 1), (fi, 0, 1)]
    return mSetToD2(mCoord, coeff)

def mD2Dodecahedron(r = 1.):
    # Додекаэдр - 20
    fi = (1.+math.sqrt(5))/2
    coeff = r/math.sqrt(3)
    mCoord = [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1), (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1),
              (0, -1/fi, -fi), (0, -1/fi, fi), (0, 1/fi, -fi), (0, 1/fi, fi),
              (-1/fi, -fi, 0), (-1/fi, fi, 0), (1/fi, -fi, 0), (1/fi, fi, 0),
              (-fi, 0, -1/fi), (-fi, 0, 1/fi), (fi, 0, -1/fi), (fi, 0, 1/fi)]
    return mSetToD2(mCoord, coeff)

'''Последовательности'''

def mSetGrid1(N, a, offset):
    # Одномерная сетка
    mSet = []
    for i in range(N):
        mSet.append((i*a+offset,))
    return mSet

def mSetFibo(N):
    # Фибоначчи
    mSet = []
    Fibo_2, Fibo_1 = 0, 1 
    mSet.append((1,)) 
    for i in range(1, N):
        Fibo = Fibo_1 + Fibo_2 
        mSet.append((Fibo,))
        Fibo_2 = Fibo_1
        Fibo_1 = Fibo 
    return mSet
