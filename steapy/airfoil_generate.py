import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from pylab import *

'''https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil

   The NACA four-digit wing sections define the profile by
   For example, the NACA 2412 airfoil has a maximum camber of 2% located 40% (0.4 chords) 
   from the leading edge with a maximum thickness of 12% of the chord.

   The NACA 0015 airfoil is symmetrical, the 00 indicating that it has no camber. 
   The 15 indicates that the airfoil has a 15% thickness to chord length ratio: 
   it is 15% as thick as it is long.
   
   m = maximum camber
   p = loaction of maimum camber
   t = thickness of the airfoil
   c = chrod length
   
   NACA m p cc
   '''
#defining the camber line
def camber_line( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    m * (x / np.power(p,2)) * (2.0 * p - (x / c)),
                    m * ((c - x) / np.power(1-p,2)) * (1.0 + (x / c) - 2.0 * p ))

#derivative 
def dyc_over_dx( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    ((2.0 * m) / np.power(p,2)) * (p - x / c),
                    ((2.0 * m ) / np.power(1-p,2)) * (p - x / c ))

#defining the thickness
def thickness( x, t, c ):
    term1 =  0.2969 * (np.sqrt(x/c))
    term2 = -0.1260 * (x/c)
    term3 = -0.3516 * np.power(x/c,2)
    term4 =  0.2843 * np.power(x/c,3)
    term5 = -0.1015 * np.power(x/c,4)
    return 5 * t * c * (term1 + term2 + term3 + term4 + term5)

#defining the upper and lower coordinates of the airfoil
def naca_cambered(x, m, p, t, c=1):
    dyc_dx = dyc_over_dx(x, m, p, c)
    th = np.arctan(dyc_dx)
    yt = thickness(x, t, c)
    yc = camber_line(x, m, p, c)  
    y1,y2,y3,y4= x - yt*np.sin(th), yc + yt*np.cos(th),  x + yt*np.sin(th), yc - yt*np.cos(th)
    data1=hstack((y1,y3[::-1]))     # X Coordinates starting from 1 to 0 and then back to 1
    data2=hstack((y2,y4[::-1]))     # Y Coordinates starting from x=1 to x=0 and then back to x=1
    data=vstack((data1,data2)).T       
    savetxt('NACA'+'_'+str(int(m*100))+str(int(p*10))+str(int(t*100))+'.txt',data)
    return (y1,y2),(y3,y4)

def naca_symmetrical(x, t, c=1):
    yt = thickness(x, t, c)
    y1 = x
    y2 = yt
    y3 = x
    y4 = -yt
    data1=hstack((y1,y3[::-1]))     # X Coordinates starting from 1 to 0 and then back to 1
    data2=hstack((y2,y4[::-1]))     # Y Coordinates starting from x=1 to x=0 and then back to x=1
    data=vstack((data1,data2)).T       
    savetxt('NACA'+'_'+'00'+str(int(t*100))+'.txt',data)
    return (y1,y2),(y3,y4)
    
    
    

'''
#naca0012 
m = 0
p = 1
t = 0.12
c = 1.0
'''
    
##naca2412 
#m = 0.02
#p = 0.4
#t = 0.12
#c = 1.0

def coordinates_symmetrical(x, t, c):
    x = np.linspace(1, 0, 200)
    for item in naca_symmetrical(x, t, c):
        (item[0], item[1], 'b')
        

def coordinates_cambered(x,m,p,t,c):
    x = np.linspace(1,0,200)
    for item in naca_cambered(x, m, p, t, c):
        (item[0], item[1], 'b')
        
        

##plt.plot(x, camber_line(x, m, p, c), 'r')
#plt.axis('equal')
#plt.xlim((-0.05, 1.05))
#savefig('airfoil.pdf')
#show()