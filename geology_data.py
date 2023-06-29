# math libraries

import numpy as np
import sympy

from sympy import *
t = Symbol('t')

from scipy.interpolate import interp2d
from scipy.interpolate import griddata

import math

from skimage.measure import  points_in_poly

# data analysis library


import pandas as pd

from copernico import *
import utm

from Bezier import Bezier

# import plotting libraries

import plotly.offline as go_offline
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker

# imaging library

from PIL import Image

# Custom functions

def surface(grid,metodo): #metodo='linear', 'cubic' o 'nearest'

    tx=grid[0]
    ty=grid[1]
    tz=grid[2]

    X = np.linspace(min(tx), max(tx))
    Y = np.linspace(min(ty), max(ty))
    X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation
    Z = griddata((tx,ty),tz,(X,Y), 
                 method=metodo)
    
    return [tx,ty,tz,X,Y,Z]

def contact_tr(csv,elevation_gain):
    tr=pd.read_csv(csv)
    tr1=tr[['UTM_X','UTM_Y','elevation']].to_numpy()
    tx=[x[0] for x in tr1]
    ty=[x[1] for x in tr1]
    telv=[x[2]+elevation_gain for x in tr1]
    return [tx,ty,telv]

def contact_dat(x,y,elv,h,mode,nam,color,gr,t,ls):
    z=[x+h for x in elv]
    return go.Scatter3d(x=x, y=y, z=z,
            mode =mode,
            name=nam,
            legendgroup=gr,
            showlegend=t,
            line=dict(color=color,
                      dash=ls,
                      width=5)
                        )

def falls(i,falla,s,h,r,rr,dd):
    A=(falla[0][i],falla[1][i])
    c1z=falla[2][i]+h
    c2x=falla[0][i+1]
    c2y=falla[1][i+1]
    d=(c2x-A[0],c2y-A[1])
    nd=np.sqrt(d[0]**2+d[1]**2)
    dn=(d[0]/nd,d[1]/nd)
    A1=(A[0]+dn[0]*s,A[1]+dn[1]*s-dd)
    B=(A1[0]+dn[0]*250,A1[1]+dn[1]*250)
    C=(B[0]-r,B[1]+rr)
    return [[A1[0],B[0],C[0]],[A1[1],B[1],C[1]],[c1z,c1z,c1z]]

def eq(A,B): # A and B are points
    x1, y1 = A
    x2, y2 = B
    # calculate the distance between the two points
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    # calculate the angle between the two points
    angle = math.atan2(y2 - y1, x2 - x1)

    # calculate the coordinates of the third point (the first one)
    x3_1 = x1 + distance * math.cos(angle - (1 * math.pi / 3))
    y3_1 = y1 + distance * math.sin(angle - (1 * math.pi / 3))

    # calculate the coordinates of the third point (the second one)
    x3_2 = x1 + distance * math.cos(angle + (1 * math.pi / 3))
    y3_2 = y1 + distance * math.sin(angle + (1 * math.pi / 3))
    return[[x3_1,y3_1],[x3_2,y3_2]]

def tr(i,contacto,h,s,hh):
    A=(contacto[0][i],contacto[1][i])
    c1z=contacto[2][i]+h
    c2x=contacto[0][i+1]
    c2y=contacto[1][i+1]
    d=(c2x-A[0],c2y-A[1])
    nd=np.sqrt(d[0]**2+d[1]**2)
    dn=(d[0]/nd,d[1]/nd)
    B=(A[0]+dn[0]*150,A[1]+dn[1]*150)
    C=eq(A, B)[s]
    return [[A[0],B[0],C[0]],[A[1],B[1],C[1]],[c1z,c1z,c1z+hh]]

def fil_pol3(poly,grid):
    xypoly=[[poly[0][i],poly[1][i]] for i in range(len(poly[0]))]
    X1=grid[0]
    Y1=grid[1]
    Z1=np.array(grid[2])
    xy1=[[X1[i],Y1[i]] for i in range(len(X1))]

    mask1=points_in_poly(xy1,xypoly)

    Z1[np.where(~mask1)] = 0  # to avoid importing numpy.ma

    xyz1=[[X1[i],Y1[i],Z1[i]] for i in range(len(X1)) if Z1[i]!=0]

    xx1=[x[0] for x in xyz1]
    yy1=[x[1] for x in xyz1]
    zz1=[x[2] for x in xyz1]
    return [xx1,yy1,zz1]

def scatter(datos,nombre,opc,c,s,group,t):
    return go.Scatter3d(x=datos[0],
               y=datos[1],
               z=datos[2], 
               name=nombre,
               legendgroup=group,
               showlegend=t,
               mode='markers',
               marker=dict(
               size=s,
               color=datos[2],                # set color to an array/list of desired values
               colorscale=[[0, c], [1,c]],   # choose a colorscale
               opacity=opc))

def para_line(p1,p2):
    P1=np.array(p1)
    P2=np.array(p2)
    v = P2 - P1
    g = lambda t: [P1[0]+v[0]*t,P1[1]+v[1]*t,P1[2]+v[2]*t]
    return g


def get_csv(lista):
    csv,nombre,color=lista
    c1=pd.read_csv(DATADIR+csv+'.csv')
    c2=c1[['x','y','z']].to_numpy()
    cx=[x[0] for x in c2]
    cy=[x[1] for x in c2]
    cz=[x[2] for x in c2]
    return [[cx,cy,cz],nombre,color]

def plane_data(lista,grupo):
    return [go.Scatter3d(x=x[0][0], 
                         y=x[0][1], 
                         z=x[0][2],
                         mode ='lines',
                         #name=name,
                         legendgroup=grupo,
                         showlegend=False,
                         line=dict(color=x[2],width=3),
                         marker=dict(color=x[2],size=1)
                          ) for x in lista]

def plane_y(p1,p2,p3,x,z):
    P1=np.array(p1)
    P2=np.array(p2)
    P3=np.array(p3)
    v1 = P3 - P1
    v2 = P2 - P1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, P3)
    if b!=0:
        return [True,(d - a * x - c * z) / b ] 
    else:
        return [False, a*x+c*z-d]
    
def grid(A,B,r,s):
    A1=[A[0],A[1],300]
    x0=min(A[0],B[0])
    d1=abs(A[0]-B[0])//r
    x=[x0+i*d1 for i in range(r)]
    d2=(max(A[2],B[2]))//s
    z=[300+i*d2 for i in range(s)]
    y=[[plane_y(A,A1,B,x[i],z[j])[1] for i in range(r)] for j in range(s)]
    xyz=[]
    for i in range(r):
        for j in range(s):
            xyz=xyz+[(x[i],y[j][i],z[j])]
    cx=[xyz[i][0] for i in range(len(xyz))]
    cy=[xyz[i][1] for i in range(len(xyz))]
    cz=[xyz[i][2] for i in range(len(xyz))]
    return [cx,cy,cz]

def gridd(A,B,r,s):
    A1=[A[0],A[1],300]
    x0=min(A[0],B[0])
    d1=abs(A[0]-B[0])//r
    x=[x0+i*d1 for i in range(17,r+30,1)]
    d2=(max(A[2],B[2]))//s
    z=[300+i*d2 for i in range(s)]
    y=[[plane_y(A,A1,B,x[i],z[j])[1] for i in range(len(x))] for j in range(s)]
    xyz=[]
    for i in range(len(x)):
        for j in range(s):
            xyz=xyz+[(x[i],y[j][i],z[j])]
    cx=[xyz[i][0] for i in range(len(xyz))]
    cy=[xyz[i][1] for i in range(len(xyz))]
    cz=[xyz[i][2] for i in range(len(xyz))]
    return [cx,cy,cz]

def fil_pol(poly,grid):
    xzpoly=[[poly[0][i],poly[2][i]] for i in range(len(poly[0]))]
    X1=grid[0]
    Y1=np.array(grid[1])
    Z1=grid[2]
    xz1=[[X1[i],Z1[i]] for i in range(len(X1))]

    mask1=points_in_poly(xz1,xzpoly)

    Y1[np.where(~mask1)] = 0  # to avoid importing numpy.ma

    xyz1=[[X1[i],Y1[i],Z1[i]] for i in range(len(X1)) if Y1[i]!=0]

    xx1=[x[0] for x in xyz1]
    yy1=[x[1] for x in xyz1]
    zz1=[x[2] for x in xyz1]
    return [xx1,yy1,zz1]

t_points = np.arange(0, 1, 0.01) # Creates an iterable list from 0 to 1.

def bz(P, direc, L_int,L_heights,Q):
    if Q==[]:
        c=[P]+[[P[0]+direc[0]*L_int[i],
                P[1]+direc[1]*L_int[i],
                L_heights[i]] for i in range(len(L_int))]
    else:
        c=c=[P]+[[P[0]+direc[0]*L_int[i],
                P[1]+direc[1]*L_int[i],
                L_heights[i]] for i in range(len(L_int))]+[Q]
    xc=[x[0] for x in c]
    yc=[x[1] for x in c]
    zc=[x[2] for x in c]
    c_Bz=[[xc[i], zc[i]] for i in range(len(xc))]
    npc_Bz=np.array(c_Bz)
    c1_Bz= Bezier.Curve(t_points,npc_Bz) 
    x1=[x[0] for x in c1_Bz]
    z1=[x[1] for x in c1_Bz]
    y1=[plane_y(c[0],c[1],c[2],x1[i],z1[i])[1] for i in range(len(x1))]
    t1=[np.sqrt(x1[i]**2+y1[i]**2) for i in range(len(x1))]
    return [x1,y1,z1,npc_Bz,t1]  

# The equation of the plane determined by three points, z in terms of x and y

def plane(p1,p2,p3,x,y):
    P1=np.array(p1)
    P2=np.array(p2)
    P3=np.array(p3)
    v1 = P3 - P1
    v2 = P2 - P1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, P3)
    if c!=0:
        return [True,(d - a * x - b * y) / c ] 
    else:
        return [False, a*x+b*y-d]

# The equation of the plane determined by three points, y in terms of x and z

def plane_y(p1,p2,p3,x,z):
    P1=np.array(p1)
    P2=np.array(p2)
    P3=np.array(p3)
    v1 = P3 - P1
    v2 = P2 - P1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, P3)
    if b!=0:
        return [True,(d - a * x - c * z) / b ] 
    else:
        return [False, a*x+c*z-d]

# Implicit equation of the line determined by two points, y and z in terms of x 

def impl_line(p1,p2,x):
    P1=np.array(p1)
    P2=np.array(p2)
    v = P2 - P1
    y=((x-P1[0])/v[0])*v[1]+P1[1]
    z=((x-P1[0])/v[0])*v[2]+P1[2]
    return [x,y,z] 

#  parametric equation of the line determined by two points

def para_line(p1,p2):
    P1=np.array(p1)
    P2=np.array(p2)
    v = P2 - P1
    g = lambda t: [P1[0]+v[0]*t,P1[1]+v[1]*t,P1[2]+v[2]*t]
    return g

# the intersection of a line with a plane

def inter_line_plane(p1,p2,q1,q2,q3):
    r=para_line(p1,p2)
    g= lambda t:list(map(r,[t]))[0]
    x,y,z=list(map(g, [t]))[0]
    p=plane(q1,q2,q3,x,y)
    if p[0]==False:
        l=solve(p[1],t)
        sol=list(map(g,l))[0]
        return  [float(x) for x in sol]
    else:
        l=solve(z-p[1],t)
        sol=list(map(g,l))[0]
        return [float(x) for x in sol]

DATADIR='data/' # Directory with the data
FIGURESDIR='figures/' # Figures produced


# Topography data

topo=pd.read_csv(DATADIR+'topografy.csv')

#We split the data in topo to get three numpy arrays

topo1=topo[['UTM_X','UTM_Y','elevation']].to_numpy()

topox=[x[0] for x in topo1]
topoy=[x[1] for x in topo1]
topoelv=[x[2] for x in topo1]

topo_gr=[topox,topoy,topoelv]
topo_linear=surface(topo_gr,'linear')

# In the list cotasxy we save the dimensions of the model.

cotasxy=[min(topox)+200,max(topox)-2000,min(topoy)+100,max(topoy)-200]


# mechanic contacts

palomeque=contact_tr(DATADIR+'palomeque.csv',0)+['palomeque_thrust']
palomeque_ds=contact_tr(DATADIR+'palomeque_ds.csv',0)+['palomeque_desplazado']
calvillo=contact_tr(DATADIR+'calvillo.csv',0)+['calvillo_thrust']
calvillo_ds=contact_tr(DATADIR+'calvillo_ds.csv',0)+['calvillo_desplazado']

# faults

fault2=contact_tr(DATADIR+'fault2.csv',0)+['fault2']
fault1=contact_tr(DATADIR+'fault1.csv',0)+['fault1']
fault3=contact_tr(DATADIR+'fault3.csv',0)+['fault1']

# stratigraphic contacts

E12cal=contact_tr(DATADIR+'E1E2 calvillo.csv',0)+['E1-E2 calvillo']
E12pal=contact_tr(DATADIR+'E1E2 palomeque.csv',0)+['E1-E2 palomeque']
E23pal=contact_tr(DATADIR+'E2E3 palomeque.csv',0)+['E2-E3 palomeque']
E3Opal=contact_tr(DATADIR+'E3O palomeque.csv',0)+['E3-O palomeque']
PEpal=contact_tr(DATADIR+'PE palomeque.csv',0)+['P-E palomeque']
PEcal=contact_tr(DATADIR+'PE calvillo.csv',0)+['P-E calvillo']
E12cal_ds=contact_tr(DATADIR+'E1E2 calvillo_ds.csv',0)+['E1-E2 calvillo ds']
E12pal_ds=contact_tr(DATADIR+'E1E2 palomeque_ds.csv',0)+['E1-E2 palomeque ds']


cont=[palomeque_ds,palomeque,calvillo,calvillo_ds]
cab=[E12cal,E12pal,E23pal,E3Opal,PEpal,PEcal,E12cal_ds,E12pal_ds]
fal=[fault1,fault2,fault3]

contacts=[contact_dat(palomeque_ds[0],palomeque_ds[1],palomeque_ds[2],30,
                        "lines",'Thrusts','black','contacts',True,None)
          ]+[contact_dat(x[0],x[1],x[2],30,"lines",x[3],'black','contacts',False,None) for x in cont[1:]]

cabs=[contact_dat(E12cal[0],E12cal[1],E12cal[2],30,"lines",'Stratigraphic contacts and units','black',
                            'cabalgamientos',True,'dot')
              ]+[contact_dat(x[0],x[1],x[2],30,"lines",x[3],'black','cabalgamientos',False,'dot') for x in cab[1:]]

faults=[contact_dat(fault1[0],fault1[1],fault1[2],50,"lines",'Strike-slip faults','black','faults',True,None)
      ]+[contact_dat(x[0],x[1],x[2],30,"lines",x[3],'black','faults',False,None) for x in fal[1:]]

f1=falls(1,fault1,300,0,50,0,40)
f2=falls(0,fault2,-400,40,60,20,100)

namesx=[619200,620800]
namesy=[4196500,4196250]
namesz=[850,850]
namestxt=['Calvillo height','Palomeque height','A','A$^\prime$']

tx=[calvillo[0][1]+50,calvillo[0][3],calvillo[0][4],calvillo[0][5],calvillo[0][7]]
ty=[calvillo[1][1]-50,calvillo[1][3]-50,calvillo[1][4]-70,calvillo[1][5]-50,calvillo[1][7]-70]
tz=[calvillo[2][1]+40,calvillo[2][3]+30,calvillo[2][4]+50,calvillo[2][5]+15,calvillo[2][7]+20]

nnx=[618900,619350,620200,619350,619400,620120,620700,
     618850,619200,620900,621100,619700]
nny=[4196800,4196700,4196400,4196200,4196000,4196100,4196000,
     4195500,4195500,4195800,4195500,4195200]
nnz=[700,750,750,750,700,790,700,750,700,600,650,700]
nntxt=['M1','P','P','E1','E2','E1','E2','E1','E2','E3','O2','E1']

tcalvillo=[tr(1,calvillo,15,0,100),tr(3,calvillo,15,0,50),tr(5,calvillo,15,0,50),tr(7,calvillo,15,0,0)]

tcalvillo_dat=[go.Mesh3d(x=x[0],y=x[1],z=x[2],
                        alphahull=5, opacity=1, color='black',
                        i = np.array([0]),
                        j = np.array([1]),
                        k = np.array([2]),
                         legendgroup='contacts',
                         showlegend=False,
                       ) for x in tcalvillo]

tcalvillo_ds=[tr(1,calvillo_ds,25,1,30)]
tcalvillo_ds_dat=[go.Mesh3d(x=x[0],y=x[1],z=x[2],
                        alphahull=5, opacity=1, color='black',
                        i = np.array([0]),
                        j = np.array([1]),
                        k = np.array([2]),
                            legendgroup='contacts'
                       ) for x in tcalvillo_ds]

tpalomeque=[tr(1,palomeque,15,0,70),tr(4,palomeque,15,0,70),tr(5,palomeque,15,0,90),
            tr(6,palomeque,15,0,70),tr(7,palomeque,15,0,30)]

tpalomeque_dat=[go.Mesh3d(x=x[0],y=x[1],z=x[2],
                        alphahull=5, opacity=1, color='black',
                        i = np.array([0]),
                        j = np.array([1]),
                        k = np.array([2]),
                          legendgroup='contacts'
                       ) for x in tpalomeque]

tpalomeque_ds=[tr(1,palomeque_ds,5,1,0)]

tpalomeque_ds_dat=[go.Mesh3d(x=x[0],y=x[1],z=x[2],
                        alphahull=5, opacity=1, color='black',
                        i = np.array([0]),
                        j = np.array([1]),
                        k = np.array([2]), 
                        legendgroup='contacts',
                       ) for x in tpalomeque_ds]

tcont=tcalvillo_dat+tcalvillo_ds_dat+tpalomeque_dat+tpalomeque_ds_dat

# New topography

NT=pd.read_csv(DATADIR+'TN_el.csv')
NT1=NT[['x','y','elevation']].to_numpy()

NTx=[x[0] for x in NT1]
NTy=[x[1] for x in NT1]
NTelv=[x[2] for x in NT1]

gTN=[NTx,NTy,NTelv]

# Crossed sections

sections=pd.read_csv(DATADIR+'sections.csv')

sections1=sections[['UTM_X','UTM_Y','elevation','point']].to_numpy()
secx=[x[0] for x in sections1]
secy=[x[1] for x in sections1]
secelv=[x[2] for x in sections1]
secp=[x[3] for x in sections1]

# split the data
secx2=[(secx[i],secx[i+1]) for i in range(len(secx)-1)]
secy2=[(secy[i],secy[i+1]) for i in range(len(secx)-1)]
secelv2=[(secelv[i],secelv[i+1]) for i in range(len(secx)-1)]

points=[[(secx[i],secy[i],secelv[i]),secp[i]] for i in range(len(secx))]

# polygon at the upper face

pol1=[calvillo[0][:-5]+PEcal[0][::-1],
      calvillo[1][:-5]+PEcal[1][::-1],
      calvillo[2][:-5]+PEcal[2][::-1]]

pol2=[PEcal[0]+[calvillo[0][6]]+E12cal[0][::-1]+[fault1[0][1],PEcal[0][0]],
      PEcal[1]+[calvillo[1][6]]+E12cal[1][::-1]+[fault1[1][1],PEcal[1][0]]]

pol3=[E12cal[0]+calvillo[0][7:]+palomeque[0][:6][::-1]+[E12cal[0][0]],
      E12cal[1]+calvillo[1][7:]+palomeque[1][:6][::-1]+[E12cal[1][0]]]

pol4=[palomeque[0][:-2]+PEpal[0][::-1]+[palomeque[0][0]],
      palomeque[1][:-2]+PEpal[1][::-1]+[palomeque[1][0]]]

pol5=[PEpal[0]+[palomeque[0][-1]]+E12pal[0][::-1]+fault1[0][-4:-2][::-1]+[PEpal[0][0]],
      PEpal[1]+[palomeque[1][-1]]+E12pal[1][::-1]+fault1[1][-4:-2][::-1]+[PEpal[1][0]]]

pol6=[palomeque_ds[0]+E12pal_ds[0]+[palomeque_ds[0][0]],
      palomeque_ds[1]+E12pal_ds[1]+[palomeque_ds[1][0]]]

pol7=[E23pal[0]+E3Opal[0][::-1]+[620505,620029,E23pal[0][0]],
      E23pal[1]+E3Opal[1][::-1]+[4194758,4194758,E23pal[1][0]]]

pol8=[calvillo_ds[0][::-1]+[fault2[0][-2]]+E12cal_ds[0][::-1]+[calvillo_ds[0][-1]],
      calvillo_ds[1][::-1]+[fault2[1][-2]]+E12cal_ds[1][::-1]+[calvillo_ds[1][-1]]]

pol9=[[secx[0]-10,fault3[0][-1]]+calvillo_ds[0][::-1]+calvillo[0]+palomeque[0][6:]+[621420,secx[0]-10],
      [secy[2],fault3[1][-1]]+calvillo_ds[1][::-1]+calvillo[1]+palomeque[1][6:]+[secy[2],secy[2]]]

pol10=[fault3[0][-3:-1][::-1]+E12cal_ds[0]+fault1[0][1:6]+E12pal_ds[0][::-1]+palomeque_ds[0]+E12pal[0]+E23pal[0][::-1]+[620029,618332,fault3[0][-1]],
       fault3[1][-3:-1][::-1]+E12cal_ds[1]+fault1[1][1:6]+E12pal_ds[1][::-1]+palomeque_ds[1]+E12pal[1]+E23pal[1][::-1]+[4194758,4194758,fault3[1][-1]]]

pol11=[[620505]+E3Opal[0]+[cotasxy[1]],
       [4194758]+E3Opal[1]+[cotasxy[2]]]

pol_list=[pol1,pol2,pol3,pol4,pol5,pol6,pol7,pol8,pol9,pol10,pol11]
name_list=['P','E1','E2','E3','M1','O2','7','8','M1','10','Topography']
color_list=['#e6d7ff','#E97451','pink','#e6d7ff','#E97451','#C19A6B','#D1BEA8','#E97451','yellow','pink','#C19A6B']
filled=[fil_pol3(x,gTN) for x in pol_list]
boll=[False for i in range(len(pol_list)-1)]+[True]
p_dat=[scatter(filled[i],name_list[i],0.8,color_list[i],3,'ct',boll[i]) for i in range(len(pol_list))]


# The border

punto1=(min(topox)+200+1,min(topoy)+100+1,700)
punto2=(min(topox)+200+1,max(topoy)-200-1,700)
punto3=(max(topox)-2000-1,max(topoy)-201,700)
punto31=(max(topox)-2000-1,max(topoy)-150,700)
punto4=(max(topox)-2000-1,min(topoy)+100+1,700)

xp=[min(topox),min(topox),max(topox),max(topox)]
yp=[min(topoy),max(topoy),max(topoy),min(topoy)]
zp=[700,700,700,700]


# Straight line through two points

r1=para_line(punto1,punto2)
r2=para_line(punto1,punto4)
r3=para_line(punto2,punto3)
r4=para_line(punto4,punto3)
r41=para_line(punto4,punto31)

#100 points on a straight line

t_points = np.arange(0, 1, 0.01) # Creates an iterable list from 0 to 1.

l1=list(map(r1,t_points))
l2=list(map(r2,t_points))
l3=list(map(r3,t_points))
l4=list(map(r4,t_points))
l41=list(map(r41,t_points))

xl1=[x[0] for x in l1]
yl1=[x[1] for x in l1]
zl1=[x[2] for x in l1]

xl2=[x[0] for x in l2]
yl2=[x[1] for x in l2]
zl2=[x[2] for x in l2]

xl3=[x[0] for x in l3]
yl3=[x[1] for x in l3]
zl3=[x[2] for x in l3]

xl4=[x[0] for x in l4]
yl4=[x[1] for x in l4]
zl4=[x[2] for x in l4]

xl41=[x[0] for x in l41]
yl41=[x[1] for x in l41]
zl41=[x[2] for x in l41]


#into dataframes to comput latitudde and longitude

df1=pd.DataFrame(l1,columns = ['x','y','z'])
df2=pd.DataFrame(l2,columns = ['x','y','z'])
df3=pd.DataFrame(l3,columns = ['x','y','z'])
df4=pd.DataFrame(l4,columns = ['x','y','z'])

latlon1=[utm.to_latlon(xl1[i],yl1[i], 30, 'T') for i in range(len(xl1)) ]
c_lat1=[x[0] for x in latlon1]
c_lon1=[x[1] for x in latlon1]

latlon2=[utm.to_latlon(xl2[i],yl2[i], 30, 'T') for i in range(len(xl2)) ]
c_lat2=[x[0] for x in latlon2]
c_lon2=[x[1] for x in latlon2]

latlon3=[utm.to_latlon(xl3[i],yl3[i], 30, 'T') for i in range(len(xl3)) ]
c_lat3=[x[0] for x in latlon3]
c_lon3=[x[1] for x in latlon3]

latlon4=[utm.to_latlon(xl4[i],yl4[i], 30, 'T') for i in range(len(xl4)) ]
c_lat4=[x[0] for x in latlon4]
c_lon4=[x[1] for x in latlon4]

df1['lat']=c_lat1
df1['lon']=c_lon1

df2['lat']=c_lat2
df2['lon']=c_lon2

df3['lat']=c_lat3
df3['lon']=c_lon3

df4['lat']=c_lat4
df4['lon']=c_lon4

# get elevations

copernicus = CopernicusDEM(raster_paths=[DATADIR+'eu_dem_v11_E30N10.TIF'])

r1_el=copernicus.get_elevation(df1, lat_col='lat', lon_col='lon')
r2_el=copernicus.get_elevation(df2, lat_col='lat', lon_col='lon')
r3_el=copernicus.get_elevation(df3, lat_col='lat', lon_col='lon')
r4_el=copernicus.get_elevation(df4, lat_col='lat', lon_col='lon')

r1_el.to_csv(DATADIR+'r1_el.csv',index=False)
r2_el.to_csv(DATADIR+'r2_el.csv',index=False)
r3_el.to_csv(DATADIR+'r3_el.csv',index=False)
r4_el.to_csv(DATADIR+'r4_el.csv',index=False)


# with this we can calculate the topographic profile in the vertical plane determined by point1 and point2.

elvl1= r1_el['elevation'].tolist()
elvl2= r2_el['elevation'].tolist()
elvl3= r3_el['elevation'].tolist()
elvl4= r4_el['elevation'].tolist()


elvl41=elvl4

xl21=xl2+[xl4[-1]]
yl21=xl2+[yl4[-1]]
elvl21=elvl2+[elvl4[-1]]

topo_data=[go.Surface(x=topo_linear[3],y=topo_linear[4],z=topo_linear[5],
                    colorscale='YlGn',
                    opacity = 0.8,
                    name='Topography',
                    legendgroup='topo',
                    showlegend=True,
                    showscale=False),
           go.Surface(x=topo_linear[3],y=topo_linear[4],z=topo_linear[5],
                    colorscale='gray',
                    opacity = 0.3,
                    name='Topo',
                    legendgroup='topo',
                    showlegend=False,
                    showscale=False)]

text_data =[go.Scatter3d(
                   x=namesx, y=namesy, z=namesz,
                   text=namestxt,
                   name='Heights',
                   mode="text",
                   textfont=dict(color=["black","black"],size=13),
                   hoverinfo="skip"),
            go.Scatter3d(
                   x=nnx, y=nny, z=nnz,
                   text=nntxt,
                   name='Stratigraphic levels',
                   mode="text",
                   legendgroup='cabalgamientos',
                   showlegend=False,
                   textfont=dict(color='black',size=18),
                   hoverinfo="skip")]

# geological symbos for faults

geo_symbols=[go.Scatter3d(x=f1[0],
                        y=f1[1],
                        z=f1[2],
                       mode ='lines',
                       name='Sections',
                      legendgroup='faults',
                      showlegend=False,
                      line=dict(color='black',
                      width=5)
                       ),
             go.Scatter3d(x=f2[0],
                        y=f2[1],
                        z=f2[2],
                       mode ='lines',
                       name='Sections',
                      legendgroup='faults',
                      showlegend=False,
                      line=dict(color='black',
                      width=5)
                       )]

# outline the limits of the topography

outline_dat=[go.Scatter3d(x=xl1,
                        y=yl1,
                        z=elvl1,
                       mode ='lines',
                       name='Border',
                      legendgroup='border',
                      showlegend=True,
                      line=dict(color='black',
                      width=3)
                       ),
             go.Scatter3d(x=xl2,
                        y=yl2,
                        z=elvl2,
                       mode ='lines',
                       name='sur',
                      legendgroup='border',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ),
             go.Scatter3d(x=xl3,
                        y=yl3,
                        z=elvl3,
                       mode ='lines',
                       name='norte',
                      legendgroup='border',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ),
             go.Scatter3d(x=xl4,
                        y=yl4,
                        z=elvl4,
                       mode ='lines',
                       name='este',
                      legendgroup='border',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ),
             go.Scatter3d(x=[cotasxy[0],cotasxy[1],cotasxy[1],cotasxy[0],cotasxy[0]],
                        y=[cotasxy[2],cotasxy[2],cotasxy[3],cotasxy[3],cotasxy[2]],
                        z=[300,300,300,300,300],
                       mode ='lines',
                       name='East side',
                      legendgroup='upper',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ), 
             go.Scatter3d(x=[cotasxy[0],cotasxy[0]],
                        y=[cotasxy[2],cotasxy[2]],
                        z=[300,660],
                       mode ='lines',
                       name='East side',
                      legendgroup='upper',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ),
             go.Scatter3d(x=[cotasxy[0],cotasxy[0]],
                        y=[cotasxy[3],cotasxy[3]],
                        z=[300,686],
                       mode ='lines',
                       name='East side',
                      legendgroup='upper',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ), 
             go.Scatter3d(x=[cotasxy[1],cotasxy[1]],
                        y=[cotasxy[3],cotasxy[3]],
                        z=[300,562],
                       mode ='lines',
                       name='East side',
                      legendgroup='upper',
                      showlegend=False,
                      line=dict(color='black',
                      width=3)
                       ), 
             go.Scatter3d(x=[cotasxy[1],cotasxy[1]],
                        y=[cotasxy[2],cotasxy[2]],
                        z=[300,575],
                       mode ='lines',
                       name='Borders',
                      legendgroup='upper',
                      showlegend=True,
                      line=dict(color='black',
                      width=3)
                       )]


