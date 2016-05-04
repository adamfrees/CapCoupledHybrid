from numpy import *
from scipy.linalg import eigh
from scipy.integrate import ode
from time import time
from simulate_pulses import *
from numpy.random import normal


d1=0.64068
d2=0.58071
singleDipole = array([[1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]])
eps1=117.4+75./2.
Est1 = 52.
params1=(Est1,d1*Est1,d2*Est1)
eps2=101.7+75./2.
Est2 = 45.
params2=(Est2,d1*Est2,d2*Est2)
coupling=75.

ham=tensorProd(hybrid_Ham(eps1,params1),eye(3))+tensorProd(eye(3),hybrid_Ham(eps2,params2))+coupling/2.*tensorProd(singleDipole,singleDipole)
coupling = coupling/2.*tensorProd(singleDipole,singleDipole)
eigSolve = eigh(ham)
for x in range(eigSolve[1].shape[0]):
    eigSolve[1][:,x] = sign(eigSolve[1][argmax(abs(eigSolve[1][:,x])),x])*eigSolve[1][:,x]
initialState = (1./sqrt(2))*eigSolve[1][:,1].astype(complex)+(1./sqrt(2))*eigSolve[1][:,0].astype(complex)

freq1=((eigSolve[0][2]-eigSolve[0][0])+(eigSolve[0][3]-eigSolve[0][1]))/2.
freq2=((eigSolve[0][1]-eigSolve[0][0])+(eigSolve[0][3]-eigSolve[0][2]))/2.
freq1=(eigSolve[0][2]-eigSolve[0][0])
freq1=(eigSolve[0][3]-eigSolve[0][1])

#print eigSolve[0][1]-eigSolve[0][0],eigSolve[0][2]-eigSolve[0][0],eigSolve[0][3]-eigSolve[0][0]
#print dot(eigSolve[1][:,:4].T,dot(tensorProd(singleDipole,eye(3)),eigSolve[1][:,:4]))

#print freq1
#print freq2

toPlot = [[],[],[],[],[],[],[],[],[],[]]
startTime = 0.#endTime
endTime = 40.#50.#(float(int((startTime+gateTime)*freq/(2.*pi*hbar)))+1.)*(2.*pi*hbar)/freq

def perturb(t):
    mat = tensorProd(singleDipole,eye(3))
    freq = (freq1)/hbar
    return ham+sin(freq*t)*mat*.1
trueStates=[]
for stateNum in range(4):
    initialState=eigSolve[1][:,stateNum].astype(complex)
    trueState, currentTime, toPlot = simulate(initialState, ham,perturb,toPlot,startTime,endTime)
    trueStates += [trueState]
    for i in range(1,10):
        print i-1,toPlot[i][len(toPlot[i])-1]
for iters in range(101):
    d1 = normal(scale=2.5)
    eps1=117.4+75./2.+d1
    d2 = normal(scale=5.)
    eps2=101.7+75./2.+d2
    ham=tensorProd(hybrid_Ham(eps1,params1),eye(3))+tensorProd(eye(3),hybrid_Ham(eps2,params2))+coupling/2.*tensorProd(singleDipole,singleDipole)
    eigSolve = eigh(ham)
    for x in range(eigSolve[1].shape[0]):
        eigSolve[1][:,x] = sign(eigSolve[1][argmax(abs(eigSolve[1][:,x])),x])*eigSolve[1][:,x]
    fidelity=0.
    for stateNum in range(4):
        initialState=eigSolve[1][:,stateNum].astype(complex)
        state, currentTime, toPlot = simulate(initialState, ham,perturb,toPlot,startTime,endTime)
        for i in range(1,10):
            print i-1,toPlot[i][len(toPlot[i])-1]
        print d1,d2,abs(dot(trueStates[stateNum],state))
        fidelity+=abs(dot(trueStates[stateNum],state))
    fidelity=fidelity/4.
    print 'fidelity:',fidelity
'''
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

for i in range(1,5):
    #plt.plot(toPlot[0],map(lambda x: abs(x)**2,toPlot[i]),linewidth=2.)
    plt.plot(toPlot[0],map(lambda x: x,toPlot[i]),linewidth=2.)
newList = []
for i in range(len(toPlot[0])):
    leakage = 0.
    for j in range(5,10):
        leakage += abs(toPlot[j][i])**2
    newList += [leakage]
totalProb = 0.
for i in range(1,10):
    print i-1,toPlot[i][len(toPlot[i])-1]
    totalProb += toPlot[i][len(toPlot[i])-1]#abs(toPlot[i][len(toPlot[i])-1])**2
print 'Total Probability:',totalProb
plt.plot(toPlot[0],newList,linewidth=2.)
#plt.xlim(0.,0.15)
labels = ['00','01','10','11','Leakage']
plt.legend(labels,loc='center right')
plt.xlabel('Time (ns)')
plt.ylabel('State Population')
plt.title('Nothing')
plt.show()'''
