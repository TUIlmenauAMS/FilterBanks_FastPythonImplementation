#Script for calling the online MDCT implementation analysis and synthesis, for time measurement 
#and as simple implementation example.
#Gerald Schuller, November 2017

from pyrecplayfastMDCT import *
import sound
import os
import time

N=1024 #Number of subbands and block size
#initialize filter bank memory:
initFB(N)

fb=np.sin(np.pi/(2*N)*(np.arange(0,int(1.5*N))+0.5))

X, fs= sound.wavread('test.wav')
print("X.shape:",X.shape) 
print("fs=", fs)
X=X*1.0/2**15
blocks=len(X)/N
y=np.zeros((blocks,N))
xrek=np.zeros(blocks*N)
startime=time.time()
#analysis filter bank:
for m in range(blocks):
  #print("X[m*N+np.arange(N).shape:", X[m*N+np.arange(N)].shape)
  #Analysis Filter Bank:
  y[m,:]=LDFB(X[m*N+np.arange(N)],fb)
  #Synthesis Filter Bank:
  xrek[m*N+np.arange(N)]=LDFBinv(y[m,:],fb);

endtime=time.time()
print("Duration analysis-synthesis: ", endtime-startime)

os.system('espeak -ven -s 120 '+'"The output of the synthesis MDCT"')
sound.sound(2**15*xrek,fs)



