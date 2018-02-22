#Script for calling the online MDCT implementation analysis and synthesis, for time measurement 
#and as simple implementation example.
#Gerald Schuller, November 2017

from pyrecplayfastLDFB import *
import sound
import os
import time

N=512 #Number of subbands and block size
#initialize filter bank memory:
initFB(N)

fb=np.loadtxt('fb2048t1023d512bbitcs.mat')
#N=1024:
#fb=np.loadtxt('fb4096t2047d1024bbitc.mat')

X, fs= sound.wavread('teststereo.wav')
print("X.shape:",X.shape) 
print("fs=", fs)
X=X*1.0/2**15
blocks=len(X)/N
y=np.zeros((blocks,2,N))
xrek=np.zeros((blocks*N,2))
startime=time.time()
#analysis filter bank:
for m in range(blocks):
  #print("X[m*N+np.arange(N).shape:", X[m*N+np.arange(N)].shape)
  #Analysis Filter Bank:
  #Left channel:
  y[m,0,:]=LDFB(X[m*N+np.arange(N),0],fb)
  #Right channel:
  y[m,1,:]=LDFB(X[m*N+np.arange(N),1],fb)
  #Synthesis Filter Bank:
  xrek[m*N+np.arange(N),0]=LDFBinv(y[m,0,:],fb);
  xrek[m*N+np.arange(N),1]=LDFBinv(y[m,1,:],fb);

endtime=time.time()
print("Duration analysis-synthesis: ", endtime-startime)

os.system('espeak -ven -s 120 '+'"The output of the synthesis MDCT"')
sound.sound(2**15*xrek,fs)



