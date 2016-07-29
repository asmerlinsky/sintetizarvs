
import matplotlib.pyplot as plt
import scipy.signal as signal
from scipy.io import wavfile
import numpy as np
from matplotlib import rcParams
import matplotlib
import sys
import sounddevice as sd
import time



 
 

###CARGO LOS ARCHIVOS PARA GRAFICAR
vs=np.loadtxt(sys.argv[1])
sintetizado=np.loadtxt('sintetizado.'+sys.argv[1]+'.dat')
envolvente=np.loadtxt('envolvente.'+sys.argv[1]+'.dat')
filtrado1=np.loadtxt('ptraqout.'+sys.argv[1]+'.dat')
pfinaldia=np.loadtxt('pfinal.'+sys.argv[1]+'.dat')
envcanto=np.loadtxt('envolvente.'+sys.argv[2]+'.dat')
canto=np.loadtxt(sys.argv[2])
pfinalnoche=np.loadtxt('pfinal.'+sys.argv[3]+'.dat')

##FILTRO BAJA CONTINUA
fn=44150      
a1, a2 = signal.butter(2,300/fn, btype='highpass', analog=False, output='ba')
# sintetizado=signal.filtfilt(a1,a2,sintetizado)
filtrado1=signal.filtfilt(a1,a2,filtrado1)
envcantofiltrada=signal.filtfilt(a1,a2,envcanto)
pfinaldia=signal.filtfilt(a1,a2,pfinaldia)
pfinalnoche=signal.filtfilt(a1,a2,pfinalnoche)

#normalizo
pfinalnoche/=np.max(np.abs(pfinalnoche))
sintetizado/=np.max(np.abs(sintetizado))
filtrado1/=np.max(np.abs(filtrado1))
pfinaldia/=np.max(np.abs(pfinaldia))


Fs=44150
t=np.linspace(0,len(sintetizado)/Fs,len(sintetizado))

#print("sintetizado max=",np.max(sintetizado))
#print("sintetizado min=",np.min(sintetizado))

if(np.max(sintetizado)==np.min(sintetizado)):
    sys.exit("MINIMO=MAXIMO")

    
f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20,10),dpi=80)
ax1.set_title('EMG')
ax1.plot(t[0:len(vs)],vs)
# ax1.set_ylim([np.min(sintetizado),np.max(sintetizado)])    
ax2.plot(t,envolvente)
ax2.set_title('envolvente')  

cmap = plt.get_cmap('Greys')

 
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True,figsize=(20,10),dpi=80)
ax1.set_title('envolvente canto')
ax1.plot(t,envcanto)
ax2.plot(t,envcantofiltrada)
ax3.specgram(envcantofiltrada,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)


f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(20,10),dpi=80)
ax1.plot(t,sintetizado)
ax1.set_title('vs')
ax2.specgram(sintetizado,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)

ax1.set_ylim([np.min(sintetizado),np.max(sintetizado)])






f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20,10))
ax1.plot(t[0:len(vs)],vs)
ax1.set_title('señal sintetizada')
ax2.set_title('sonograma')
ax2.specgram(sintetizado,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax2.set_ylim([0,10000])


f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20,10))
ax1.plot(t[0:len(filtrado1)],filtrado1)
ax1.set_title('P salida de la traquea')
ax2.set_title('sonograma')
ax2.specgram(filtrado1,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax2.set_ylim([0,10000])





f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(20,10))
ax1.set_title('pfinaldia')
ax1.plot(t[0:len(pfinaldia)],pfinaldia)
Pxx, freqs, bins, im = ax2.specgram(pfinaldia,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax2.set_ylim([0,10000])



sd.play(canto,44150)
time.sleep(1)

sd.play(pfinaldia, 44150)
time.sleep(1)

sd.play(pfinalnoche, 44150)

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True,sharey=True, figsize=(20,10))
plt.suptitle('SONOGRAMAS')
Pxx, freqs, bins, im = ax1.specgram(canto,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)
ax1.set_ylim([0,10000])
ax1.set_title('canto')
Pxx, freqs, bins, im = ax2.specgram(pfinaldia,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)

ax2.set_ylim([0,10000])
plt.xlabel("Tiempo (s)")
ax2.set_title('salida dia')
Pxx, freqs, bins, im = ax3.specgram(pfinalnoche,pad_to=600,cmap=cmap, NFFT=220, Fs=44150,noverlap=200)

ax3.set_ylim([0,10000])
ax3.set_title('salida noche')
wavfile.write('nochesint.wav',44150,np.int16(32767*pfinalnoche))
wavfile.write('diasint.wav',44150,np.int16(32767*pfinaldia))


plt.show(block=True)
