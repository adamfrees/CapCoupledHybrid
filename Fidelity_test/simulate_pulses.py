from numpy import *
from scipy.linalg import eigh
from scipy.integrate import ode
from time import time

hbar = 0.6582119 #ueV*ns

def hybrid_Ham(eps,params):
    '''
    4x4 matrix that describes a hybrid qubit
    '''
    (dER,delta1,delta2) = params
    Ham = array([[eps/2.,delta1,-delta2],[delta1,-eps/2.,0.],[-delta2,0.,-eps/2.+dER]])
    return Ham

def tensorProd(a,b):
    (a1size,a2size) = a.shape
    (b1size,b2size) = b.shape
    out = zeros((a1size*b1size,a2size*b2size))
    for i in range(a1size):
        for j in range(a2size):
            out[i*b1size:(i+1)*b1size,j*b2size:(j+1)*b2size] = a[i,j]*b
    return out

def rotatingFrame(toTransform,frame,time):
    eigSolve = eigh(frame)
    for x in range(eigSolve[1].shape[0]):
        eigSolve[1][:,x] = sign(eigSolve[1][argmax(abs(eigSolve[1][:,x])),x])*eigSolve[1][:,x]
    timeU = diag(map(lambda x: exp(1.j*(x)*time/hbar),eigSolve[0]))#eigSolve[0]))
    #return dot(timeU,dot(eigSolve[1].T,toTransform))
    return dot(timeU,dot(eigSolve[1].T,toTransform))

def schrod(t,psi,hamiltonian):
    try:
        a = hamiltonian(t)
    except TypeError:
        return -1.j/hbar*(dot(hamiltonian,psi))
    return -1.j/hbar*(dot(hamiltonian(t),psi))

def simulate(initial, staticHam,operation,toPlot,startTime,endTime):
    from progress.bar import Bar
    r = ode(schrod).set_integrator('zvode',method='bdf',atol=1.e-12,rtol=1.e-12)
    r.set_initial_value(initial,startTime)
    r.set_f_params(operation)

    dt = 1.e-4#(0.001*pi*hbar)#/freqtot#.001
    counter = 1
    bar = Bar('Processing', max=int((endTime-startTime)/dt))
    start = time()
    while r.successful() and r.t<endTime:
        r.integrate(r.t+dt)
        bar.next()
        if counter == 20:
            toPlot[0] += [r.t]
            newY = rotatingFrame(r.y,staticHam,r.t)
            for i,x in enumerate(newY):
                toPlot[i+1] += [abs(x)**2]#atan2(x.imag,x.real)]
            counter =0
        counter +=1
    bar.finish()
    print time()-start
    return r.y,r.t,toPlot
