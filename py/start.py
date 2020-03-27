# Метод пристрелки
from numpy  import*
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm,os
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy.integrate import odeint
from scipy import linalg
import time
start = time.time()
c1 = 1.0
c2 = 0.8
c3 = 0.5
a =0.0
b = 1.0
nn =100
initial_state_0 =array( [a, c1, 0.0, 0.0])
initial_state_I =array( [a, 0.0, 1.0, 0.0])
initial_state_II =array( [a, 0.0, 0.0, 1.0])
to = a
tEnd =b
N = int(nn)
tau=(b-a)/N
def rungeKutta(f, to, yo, tEnd, tau):
         def increment(f, t, y, tau):
                  k1=tau*f(t,y)
                  k2=tau*f(t+(1/4)*tau,y+(1/4)*k1)
                  k3 =tau *f(t+(3/8)*tau,y+(3/32)*k1+(9/32)*k2)
                  k4=tau*f(t+(12/13)*tau,y+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3)
                  k5=tau*f(t+tau,y+(439/216)*k1-8*k2+(3680/513)*k3 -(845/4104)*k4)
                  k6=tau*f(t+(1/2)*tau,y-(8/27)*k1+2*k2-(3544/2565)*k3 +(1859/4104)*k4-(11/40)*k5)
                  return (16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6  
         t = []
         y= []
         t.append(to)
         y.append(yo)
         while to < tEnd:
                  tau = min(tau, tEnd - to)
                  yo = yo + increment(f, to, yo, tau)
                  to = to + tau
                  t.append(to)
                  y.append(yo)         
         return array(t), array(y)
def f(t, y):
         global theta
         f = zeros([4])
         f[0] = 1
         f[1] = -y [1]-y[2] +theta* sin(y[0])
         f[2] = -y[2]+y[3]
         f[3] = -y[2]
         return f
# Решение НЕОДНОРОДНОЙ системы -- theta = 1
theta = 1.0
yo =initial_state_0
t, y = rungeKutta(f, to, yo, tEnd, tau)
y2=[i[2] for i in y]
y3=[i[3] for i in y]
# Извлечение требуемых для решения задачи значений 
# Y20 = Y2(b), Y30 = Y3(b)
Y20 = y2[N-1]
Y30 = y3[N-1]
# Решение ОДНОРОДНОЙ системы -- theta = 0, задача I
theta = 0.0
yo= initial_state_I
t, y = rungeKutta(f, to, yo, tEnd, tau)
y2=[i[2] for i in y]
y3=[i[3] for i in y]
# Извлечение требуемых для решения задачи значений 
# Y21= Y2(b), Y31 = Y3(b)
Y21= y2[N-1]
Y31 = y3[N-1]
# Решение ОДНОРОДНОЙ системы -- theta = 0, задача II
theta = 0.0
yo =initial_state_II
t, y = rungeKutta(f, to, yo, tEnd, tau)
y2=[i[2] for i in y]
y3=[i[3] for i in y]
# Извлечение требуемых для решения задачи значений 
# Y211= Y2(b), Y311 = Y3(b)
Y211= y2[N-1]
Y311 = y3[N-1]
# Формирование системы линейных
# АЛГЕБРАИЧЕСКИХ уравнений для определния p2, p3
b1 = c2 - Y20
b2 = c3 - Y30
A = array([[Y21, Y211], [Y31, Y311]])
bb = array([[b1], [b2]])
 # Решение системы
p2, p3 = linalg.solve(A, bb)
# Окончательное решение краевой
# НЕОДНОРОДНОЙ задачи, theta = 1
theta = 1.0
yo = array([a, c1, p2, p3])
t, y = rungeKutta(f, to, yo, tEnd, tau)
y0=[i[0] for i in y]
y1=[i[1] for i in y]
y2=[i[2] for i in y]
y3=[i[3] for i in y]
# Проверка
print('y0[0]=', y0[0])
print('y1[0]=', y1[0])
print('y2[0]=', y2[0])
print('y3[0]=', y3[0])
print('y0[N-1]=', y0[N-1])
print('y1[N-1]=', y1[N-1])
print('y2[N-1]=', y2[N-1])
print('y3[N-1]=', y3[N-1])
j = N
xx = y0[:j]
yy1 = y1[:j]
yy2 = y2[:j]
yy3 = y3[:j]
stop = time.time()
print ("Время на модельную задачу: %f"%(stop-start))
plt.subplot(2, 1, 1)
plt.plot([a], [c1], 'ro')
plt.plot([b], [c2], 'go')
plt.plot([b], [c3], 'bo')
plt.plot(xx, yy1, color='r') #Построение графика
plt.plot(xx, yy2, color='g') #Построение графика
plt.plot(xx, yy3, color='b') #Построение графика
plt.xlabel(r'$x$') #Метка по оси x в формате TeX
plt.ylabel(r'$y_k(x)$') #Метка по оси y в формате TeX
plt.title(r'Метод пристрелки ', color='blue')
plt.grid(True) #Сетка
patch_y1 = mpatches.Patch(color='red',
                          label='$y_1$')
patch_y2 = mpatches.Patch(color='green',
                          label='$y_2$')
patch_y3 = mpatches.Patch(color='blue',
                           label='$y_3$')
plt.legend(handles=[patch_y1, patch_y2, patch_y3])
ymin, ymax = plt.ylim()
xmin, xmax = plt.xlim()
plt.subplot(2, 1, 2)
font = {'family': 'serif',
        'color': 'blue',
        'weight': 'normal',
        'size': 12,
        }
plt.text(0.2, 0.8, r'$\frac{dy_1}{dx}= - y_1 - y_2 + \sin(x),$',
         fontdict=font)
plt.text(0.2, 0.6,r'$\frac{dy_2}{dx}= - y_1 + y_3,$',
         fontdict=font)
plt.text(0.2, 0.4, r'$\frac{dy_3}{dx}= - y_2 - y_2,$',
         fontdict=font)
plt.text(0.2, 0.2, r'$y_1(a)=c_1, '
         r'\quad y_2(b)=c_2, \quad y_3(b)=c_3.$',
         fontdict=font)
plt.subplots_adjust(left=0.15)
plt.show()