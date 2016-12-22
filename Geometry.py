# -*- coding:utf-8 -*-
''' Created on 30.06.2015, @author: ДМ '''
import math
import numpy as NP
#import GraphFlow as GF

def Val2Val(val, rangeFrom = [0, 100], rangeTo=[0, 1]):
    # Приведение числа из одного диапазона к другому - Normalisation
    # val - число из диапазона rangeFrom.
    # rangeTo - диапазон, к которому надо привести
    relval = (val-rangeFrom[0])/(rangeFrom[1]-rangeFrom[0])
    return rangeTo[0] + relval*(rangeTo[1] - rangeTo[0])  

def mSetToD2(mSet, coeff = 1.):
    ''' По набору координат узлов формируем матрицу квадратов расстояний'''
    NumNodes = len(mSet)
    mR = NP.zeros((NumNodes, NumNodes))
    i = 0
    for disp1 in mSet:
        j = 0
        for disp2 in mSet:
            if disp2 == disp1: 
                mR[i][j] = 0
            else:
                dr2 = 0
                for ik in range(len(disp1)):
                    dr = (disp2[ik] - disp1[ik])*coeff
                    dr2 += dr*dr 
                mR[i][j] = dr2
            j += 1
        i += 1 
    return mR

def mSetToD(mSet, coeff = 1., degree = 1.):
    ''' По набору координат узлов формируем матрицу расстояний.
    Можно управлять степенью - mSetToD(mSet, coeff = 1., degree = 3) вернет матрицу кубов расстояний
    '''
    NumNodes = len(mSet)
    mR = NP.zeros((NumNodes, NumNodes))
    i = 0
    for disp1 in mSet:
        j = 0
        for disp2 in mSet:
            if disp2 == disp1: 
                mR[i][j] = 0
            else:
                dr2 = 0
                for ik in range(len(disp1)):
                    d = (disp2[ik] - disp1[ik])
                    dr2 += d*d
                mR[i][j] = coeff*(dr2**(degree/2))
            j += 1
        i += 1 
    return mR

def mAddSet(mSet1, mSet2):
    '''Объединение двух наборов'''
    mSet = [] 
    for cSet in mSet1:
        mSet.append(cSet)
    for cSet in mSet2:
        mSet.append(cSet)
    return mSet

def mDivideSet(mSet, N=2, mult=1, accur=8):
    '''Делим заданный набор координат точек mSet на заданную кратность N и формируем новый набор. Повторяем mult раз'''
    #print(mSet)
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

''' Формирование набора узлов (точек) по квадратной решетке '''
def mSet(Length):
    mSet = [] 
    for disp1 in range(Length):
        for disp2 in range(Length):
            mSet.append((disp1, disp2))
    return mSet

def mSet2(NumRow, NumCol = 1, aRow = 1, aCol = 1, dispRow = 0, dispCol = 0):
    '''Двумерная: простая прямоугольная'''
    mSet = [] 
    for row in range(NumRow):
        for col in range(NumCol):
            mSet.append((col*aCol + dispCol*row, row*aRow + dispRow*col))
    return mSet

def vSet2vect(ev1, ev2, d1 = [-1, 1], d2 = [-1, 1], withLimit = False):
    '''Двумерная: по базисным векторам'''
    mSet = [] 
    for n1 in range(d1[0], d1[1]+1):
        x1, y1 = ev1[0]*n1, ev1[1]*n1
        for n2 in range(d2[0], d2[1]+1):
            if withLimit and abs(n1-n2) > d1[1]: continue 
            x2, y2 = ev2[0]*n2, ev2[1]*n2
            mSet.append((x1 + x2, y1 + y2))
    return mSet

def vSetHex(r = 1, size = [1, 1]):
    '''Гексагональная: сота
    r - шаг соты, минимальное расстояние между узлами
    size[] - размеры сетки
    '''
    n1, n2 = size[0]-1, size[1]-1
    eSet = mSetRegular(3, r)
    vSet = vSet2vect(eSet[0], eSet[1], [-n1, n1], [-n2, n2], True)
    #print(vSet)
    return vSet

def mSet3(Length1, Length2 = 1, Length3 = 1):
    ''' Трехмерная '''
    mSet = [] 
    for x in range(Length1):
        for y in range(Length2):
            for z in range(Length3):
                mSet.append((x, y, z))
    return mSet

def mSetRegular(nPoints = 3, r = 1., cc = [0, 0, 0], mult = 1, withZero = False, accur=8):
    ''' Координаты правильных многоугольников
    nPoints - количество вершин,
    r - внешний радиус
    cc - координаты центра
    mult - количество размножений координат в направлении радиуса'''
    mSet = []
    xc, yc, zc = cc[0], cc[1], cc[2]
    fi = 2*NP.pi/float(nPoints)
    for i in range(nPoints):
        cr = r
        for m in range(mult):
            x = round(xc + cr*NP.cos(fi*i), accur)
            y = round(yc + cr*NP.sin(fi*i), accur)
            mSet.append((x, y, zc))
            cr += r

    if withZero: mSet.append((0, 0, 0))
    return mSet

''' Координаты правильного креста'''
def mSetCross(N = 2, Length = 1., D = 2):
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

''' Координаты прямого угла'''
def mSetAngle(N = 1, a = 1.):
    mSet = []
    for i in range(N):
        mSet.append((i*a, 0))
        if i == 0:
            continue

        mSet.append((0, i*a))
    return mSet

''' Координаты правильных многогранников, r - внешний радиус'''

''' Тетраэдр - 4'''
def mD2Tetrahedron(r = 1.):
    fi = 1./math.sqrt(2)
    coeff = r/math.sqrt(3./8.)/2
    mCoord = [(-1, 0, -fi), (1, 0, -fi),
            (0, -1, fi), (0, 1, fi)]
    return mSetToD2(mCoord, coeff)

''' Октаэдр - 8'''
def mD2Octahedron(r = 1.):
    coeff = r
    mCoord = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)]
    return mSetToD2(mCoord, coeff)

''' Куб - 8'''
def mD2Cube(r = 1.):
    coeff = r/math.sqrt(3.)
    mCoord = [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1),
              (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1)]
    return mSetToD2(mCoord, coeff)

''' Икосаэдр - 12'''
def mD2Icosahedron(r = 1.):
    fi = (1.+math.sqrt(5))/2
    coeff = r/math.sqrt(fi*math.sqrt(5))
    mCoord = [(0, -1, -fi), (0, -1, fi),(0, 1, -fi), (0, 1, fi),
            (-1, -fi, 0), (-1, fi, 0),(1, -fi, 0), (1, fi, 0),
            (-fi, 0, -1), (fi, 0, -1),(-fi, 0, 1), (fi, 0, 1)]
    return mSetToD2(mCoord, coeff)

''' Додекаэдр - 20'''
def mD2Dodecahedron(r = 1.):
    fi = (1.+math.sqrt(5))/2
    coeff = r/math.sqrt(3)
    mCoord = [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1), (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1),
              (0, -1/fi, -fi), (0, -1/fi, fi), (0, 1/fi, -fi), (0, 1/fi, fi),
              (-1/fi, -fi, 0), (-1/fi, fi, 0), (1/fi, -fi, 0), (1/fi, fi, 0),
              (-fi, 0, -1/fi), (-fi, 0, 1/fi), (fi, 0, -1/fi), (fi, 0, 1/fi)]
    return mSetToD2(mCoord, coeff)

'''Последовательности'''

''' Одномерная сетка '''
def mSetGrid1(N, a, offset):
    mSet = []
    for i in range(N):
        mSet.append((i*a+offset,))
    return mSet

''' Фибоначчи '''
def mSetFibo(N):
    mSet = []
    Fibo_2 = 0 
    Fibo_1 = 1
    mSet.append((1,)) 
    for i in range(1, N):
        Fibo = Fibo_1 + Fibo_2 
        mSet.append((Fibo,))
        Fibo_2 = Fibo_1
        Fibo_1 = Fibo 
        
    return mSet

def RotateData(vDataX, vDataY):
    '''Стабилизирует данные относительно осей по первой точке'''
    #indX = vDataX.argmax()
    indX = 0 
    cX, cY = vDataX[indX], vDataY[indX]
    cL = NP.sqrt(cX*cX + cY*cY)
    if cL == 0: return vDataX, vDataY, cL
    if cX < 0:
        cosFi, sinFi = -cX/cL, cY/cL
    else:
        cosFi, sinFi = cX/cL, -cY/cL
    #print(cosFi, sinFi)
    mRotate = [(cosFi, -sinFi), (sinFi, cosFi)]
    mData = [vDataX, vDataY]
    mResult = NP.dot(mRotate, mData)
    
    return mResult[0], mResult[1], cL
