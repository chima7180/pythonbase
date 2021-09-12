#!/usr/bin/env python
# coding=utf-8

import numpy as np
import scipy.signal
import random
import itertools
import matplotlib
import matplotlib.pyplot as plt

F_ECH = 1e6
F_PORTEUSE = 50e3
F_MESSAGE = 5e3#2500#10000
T = 2

temps = np.arange(0, T, 1.0 / F_ECH)


def generer_porteuse(frequence=F_PORTEUSE,  phi=0):
    return np.sin(2 * np.pi * (temps * frequence) + phi)


def generer_porteuse_q(frequence=F_PORTEUSE,  phi=0):
    return np.cos(2 * np.pi * (temps * frequence) + phi)


def generer_NRZ(frequence=F_MESSAGE):
    resultat = np.zeros(temps.shape)
    for i, x in enumerate(temps):
        if (i % (F_ECH / frequence) == 0):
            etat = 1.0 * random.randint(0, 1)
        resultat[i] = etat
    return resultat


def modulation_ASK(source, k=1):
    porteuse = generer_porteuse()
    resultat = porteuse * (1 + k * source)
    return resultat


def modulation_FSK(message, ratio=1.5):
    porteuse_0 = generer_porteuse()
    porteuse_1 = generer_porteuse(frequence=ratio * F_PORTEUSE)
    result = porteuse_0 * (message == 0) + porteuse_1 * (message == 1)
    return result


def modulation_PSK(source):
    porteuse = generer_porteuse()
    resultat = porteuse * ((2 * source) - 1)
    return resultat


def calculer_tf(signal):
    tf = np.fft.rfft(signal * scipy.signal.blackmanharris(temps.shape[0]))
    frequences = np.fft.rfftfreq(2 * len(tf) - 1, 1.0 / F_ECH)
    return frequences, 20 * np.log10(np.abs(tf))


def filtre_passe_bas(signal_entree, frequence_coupure, ordre=4):
    b, a = scipy.signal.iirfilter(
        ordre, (2 * frequence_coupure / F_ECH), btype='lowpass')
    signal_sortie = scipy.signal.filtfilt(b, a, signal_entree)
    return signal_sortie


def demodulation_ASK(source, frequence_coupure):
    porteuse = generer_porteuse()
    produit = source * porteuse
    produit_filtre = filtre_passe_bas(produit, frequence_coupure)
    return (2.0 * produit_filtre)


def construire_sequence_bit(message, seuil=0.5):
    # Le message est échantilloné au milieu de chaque bit
    instants = []
    chaine = []
    for i in range(int(len(temps) / (F_ECH / F_MESSAGE))):
        pos = int(i * (F_ECH / F_MESSAGE) + ((F_ECH / F_MESSAGE) / 2))
        instants.append(temps[pos])
        chaine.append(message[pos] > seuil)
    return instants, chaine


def compute_BER(reference, message):
    error = 0
    for r, m in zip(reference, message):
        if r != m:
            error += 1
    return error / float(len(reference))


def add_noise(message, sigma=0.1):
    # FIXME : np.random.normal(0,sigma,message.shape[0])
    return message + np.random.normal(0, sigma, len(message))


def __group_bits(message, seuil=0.5, nb=2):
    instants, bits = construire_sequence_bit(message, seuil)
    _instants = np.array(instants)
    _bits = np.array(bits)
    if len(_bits) % nb != 0:
        _bits = np.append(_bits, [False] * (nb - (len(_bits) % nb)))
    ech = _instants[::nb] - (0.5 / F_MESSAGE)
    ech = np.append(ech, (ech[-1] + (nb * 1.0 / F_MESSAGE)))
    ech = ech[1:]
    grp = np.split(_bits, len(_bits) / nb)
    values = np.zeros(len(message))
    for t in range(len(grp)):
        mot = np.append([False] * (8 - nb), grp[t])
        # np.packbits(mot) ne fonctionne pas sur OSX comme sur Linux
        valeur = sum([2**(8 - idx - 1) for idx, v in enumerate(mot) if v])
        values[int(t * (nb * 1.0 / F_MESSAGE) * F_ECH) : int((t + 1) * (nb * 1.0 / F_MESSAGE) * F_ECH)] = valeur
    return values


def __modulateurIQ(signal_I, signal_Q):
    porteuse_sinus = generer_porteuse(frequence=F_PORTEUSE,  phi=0)
    porteuse_cosinus = generer_porteuse_q(frequence=F_PORTEUSE,  phi=0)

    result = (signal_I * porteuse_sinus) + (signal_Q * porteuse_cosinus)
    return result


def demodulateurIQ(signal):
    porteuse_sinus = generer_porteuse(frequence=F_PORTEUSE, phi=0)
    porteuse_cosinus = generer_porteuse_q(frequence=F_PORTEUSE, phi=0)
    signal_I_p = 2.0 * \
        filtre_passe_bas(signal * porteuse_sinus, 2.0 * F_MESSAGE)
    signal_Q_p = 2.0 * \
        filtre_passe_bas(signal * porteuse_cosinus, 2.0 * F_MESSAGE)
    return signal_I_p, signal_Q_p


def modulation_4_ASK(message, filtrage=True):
    values = (__group_bits(message, 0.5, 2) + 1) / 4.0
    signal_Q = np.zeros(len(message))
    if filtrage:
        values = filtre_passe_bas(values, F_MESSAGE/2.0)
    return __modulateurIQ(values, signal_Q)


def demodulation_4_ASK(signal_I_p, signal_Q_p):
    message = np.rint(signal_I_p * 4.0 - 1.0)
    message =  np.clip(message, 0, 3)
    return message.astype(np.uint8)


def modulation_4_PSK(message, filtrage=True):
    # Codage contient les positions de I et Q pour des valeurs allant de 0 à 3
    codage = []
    for i in range(4):
        pas = (2.0 * np.pi) / 4.0
        codage.append((np.sin(i * pas), np.cos(i * pas)))
    signal_I = np.zeros(len(message))
    signal_Q = np.zeros(len(message))
    values = __group_bits(message, 0.5, 2)
    for t, v in enumerate(values):
        (signal_I[t], signal_Q[t]) = codage[int(v)]
    if filtrage:
        signal_I = filtre_passe_bas(signal_I, F_MESSAGE/2.0)
        signal_Q = filtre_passe_bas(signal_Q, F_MESSAGE/2.0)
    return __modulateurIQ(signal_I, signal_Q)


def demodulation_4_PSK(signal_I_p, signal_Q_p):
    seuils = []
    pas = (2.0 * np.pi) / 4.0
    for i in range(4):
        seuils.append(((i + 0.5) * pas))
    #print(seuils)
    resultats = np.zeros(len(signal_I_p))
    # Calcul des angles et seuillage avec les angles
    angles = (np.arctan2(signal_I_p, signal_Q_p))
    # Conversion des angle -pi->+pi vers 0->2pi
    # http://stackoverflow.com/questions/37358016/numpy-converting-range-of-angles-from-pi-pi-to-0-2pi
    angles = (2 * np.pi + angles) * (angles < 0) + angles * (angles > 0)
    for idx, angle in enumerate(angles):
        for pos, value in enumerate(seuils):
            if angle < value:
                resultats[idx] = pos
                break
    return resultats


def modulation_16_PSK(message, filtrage = True):
    # Codage contient les positions de I et Q pour des valeurs allant de 0 à
    # 15 en code binaire... donc pas Gray
    codage = []
    for i in range(16):
        pas = (2.0 * np.pi) / 16.0
        codage.append((np.sin(i * pas), np.cos(i * pas)))
    signal_I = np.zeros(len(message))
    signal_Q = np.zeros(len(message))
    values = __group_bits(message, 0.5, 4)
    for t, v in enumerate(values):
        (signal_I[t], signal_Q[t]) = codage[int(v)]
    if filtrage:
        signal_I = filtre_passe_bas(signal_I, F_MESSAGE/4.0)
        signal_Q = filtre_passe_bas(signal_Q, F_MESSAGE/4.0)
    return __modulateurIQ(signal_I, signal_Q)


def demodulation_16_PSK(signal_I_p, signal_Q_p):
    seuils = []
    pas = (2.0 * np.pi) / 16.0
    for i in range(16):
        seuils.append(((i + 0.5) * pas))
    resultats = np.zeros(len(signal_I_p))
    # Calcul des angles et seuillage avec les angles
    angles = (np.arctan2(signal_I_p, signal_Q_p))
    # Conversion des angle -pi->+pi vers 0->2pi
    # http://stackoverflow.com/questions/37358016/numpy-converting-range-of-angles-from-pi-pi-to-0-2pi
    angles = (2 * np.pi + angles) * (angles < 0) + angles * (angles > 0)
    for idx, angle in enumerate(angles):
        for pos, value in enumerate(seuils):
            if angle < value:
                resultats[idx] = pos
                break
    return resultats


def modulation_16_QAM(message, filtrage = True):
    codage = ((-1, 1), (-1 / 3.0, 1), (1 / 3.0, 1), (1, 1),
              (-1, 1 / 3.0), (-1 / 3.0, 1 / 3.0), (1 / 3.0, 1 / 3.0), (1, 1 / 3.0),
              (-1, -1 / 3.0), (-1 / 3.0, -1 /
                               3.0), (1 / 3.0, -1 / 3.0), (1, -1 / 3.0),
              (-1, -1), (-1 / 3.0, -1), (1 / 3.0, -1), (1, -1))
    signal_I = np.zeros(len(message))
    signal_Q = np.zeros(len(message))
    values = __group_bits(message, 0.5, 4)
    for t, v in enumerate(values):
        (signal_I[t], signal_Q[t]) = codage[int(v)]
    if filtrage:
        signal_I = filtre_passe_bas(signal_I, F_MESSAGE/4.0)
        signal_Q = filtre_passe_bas(signal_Q, F_MESSAGE/4.0)
    return __modulateurIQ(signal_I, signal_Q)


def demodulation_16_QAM(signal_I_p, signal_Q_p):
    positions = ((-1, 1), (-1 / 3.0, 1), (1 / 3.0, 1), (1, 1),
                 (-1, 1 / 3.0), (-1 / 3.0, 1 /
                                 3.0), (1 / 3.0, 1 / 3.0), (1, 1 / 3.0),
                 (-1, -1 / 3.0), (-1 / 3.0, -1 /
                                  3.0), (1 / 3.0, -1 / 3.0), (1, -1 / 3.0),
                 (-1, -1), (-1 / 3.0, -1), (1 / 3.0, -1), (1, -1))
    # Conversion des position en coordonnées circulaires
    resultats = np.zeros(len(signal_I_p))
    for idx, iq in enumerate(zip(signal_I_p, signal_Q_p)):
        distances = []
        for pos, value in enumerate(positions):
            distances.append((value[0] - iq[0])**2 + (value[1] - iq[1])**2)
        val, n = min((val, n) for (n, val) in enumerate(distances))
        resultats[idx] = n
    return resultats


def decode_resultats(resultats, nb):
    # nb est la taille du groupe donc log2 de n-modulation
    message = []

    for i in range(int(len(temps) / (nb * (F_ECH / F_MESSAGE)))):
        pos = int(i * nb * (F_ECH / F_MESSAGE) +
                  (nb * (F_ECH / F_MESSAGE) / 2))
        v = bin(int(resultats[pos]))[2:]
        v = "0" * (nb - len(v)) + v
        message.extend(v)
    resultat = np.zeros(temps.shape)

    for i, x in enumerate(temps):
        if (i % (F_ECH / F_MESSAGE) == 0):
            if (i / (F_ECH / F_MESSAGE)) < len(message):
                value = message[int(i / (F_ECH / F_MESSAGE))]
            else:
                value = 0
        resultat[i] = value
    return resultat


def calculer_constellation(signal_I, signal_Q):
    plt.figure(figsize=(10,10))
    plt.hist2d(signal_I, signal_Q, bins=(100, 100), range=[[-1.1, 1.1], [-1.1, 1.1]] , cmap=plt.cm.reds)
    plt.colorbar()
    plt.grid()
    plt.show()



nrz = generer_NRZ()
nrz_fourier=calculer_tf(nrz)

plt.title('Tempo NRZ')
plt.plot(nrz)
plt.ylim(0,2)
plt.xlim(0,10000)
plt.show()
plt.title('transformer de fourrier')
plt.plot(nrz_fourier[0],nrz_fourier[1])
plt.ylim(-50,100)
plt.xlim(0,25000)
plt.show()


modulation_amplitude=modulation_ASK(nrz)
modulation_amplitude_fourier = calculer_tf(modulation_amplitude)

plt.title('Modulation Amplitude')
plt.plot(modulation_amplitude)
plt.ylim(-4,4)
plt.xlim(0,1000)
plt.show()
plt.title('transformer de fourrier')
plt.plot(modulation_amplitude_fourier[0],modulation_amplitude_fourier[1])
plt.ylim(-50,100)
plt.xlim(30000,70000)
plt.show()

nrz_filtre_passe_bas = filtre_passe_bas(nrz,F_MESSAGE)
plt.title('Temporel filtre passe bas')
plt.plot(nrz_filtre_passe_bas)
plt.plot(nrz)
plt.ylim(-0.5,1.5)
plt.xlim(0,10000)
plt.show()

#Démodulation

nrz_demoduler = demodulation_ASK(modulation_amplitude,F_MESSAGE) -1
plt.title('Temporel Demodulé')
plt.plot(nrz_demoduler)
plt.plot(nrz)
plt.ylim(-0.5,1.5)
plt.xlim(0,10000)
plt.show()

sequence_nrz_demoduler=construire_sequence_bit(nrz_demoduler)
sequence_nrz = construire_sequence_bit(nrz)
print("valeur du BER :", compute_BER(sequence_nrz_demoduler[1],sequence_nrz[1]))

#Sensibilité au bruit de la modulation ASK
k=[]
for sigma in range(5,250):
    nrz_bruit = add_noise(modulation_amplitude,sigma/100)
    nrz_bruit_demoduler = demodulation_ASK(nrz_bruit, F_MESSAGE)-1
    sequence_nrz_bruit = construire_sequence_bit(nrz_bruit)
    sequence_nrz_demoduler_2 = construire_sequence_bit(nrz_bruit_demoduler)
    k.append(compute_BER(sequence_nrz_demoduler_2[1],sequence_nrz_bruit[1]))

plt.title('Temporel Demodulé')
plt.plot(k)
plt.ylim(0.4,0.6)
plt.xlim(0,250)
plt.show()

#Modulation de fréquence
modulation_fréquence=modulation_FSK(nrz)
modulation_fréquence_fourier = calculer_tf(modulation_fréquence)
plt.title('Modulation fréquenciel')
plt.plot(modulation_fréquence)
plt.ylim(-4,4)
plt.xlim(0,1000)
plt.show()
plt.title('transformer de fourrier')
plt.plot(modulation_fréquence_fourier[0],modulation_fréquence_fourier[1])
plt.ylim(-50,100)
plt.xlim(30000,70000)
plt.show()

#modulation de phase
PSK = modulation_PSK(nrz_filtre_passe_bas)
for i in range(0,10,1):
    print(PSK[i])
PSK_fourier = calculer_tf(PSK)
plt.title('Modulation PSK')
plt.plot(PSK)
plt.ylim(-4,4)
plt.xlim(0,1000)
plt.show()
plt.title('transformer de fourrier')
plt.plot(PSK_fourier[0],PSK_fourier[1])
plt.ylim(-50,100)
plt.xlim(30000,70000)
plt.show()


v=[]
for sigma in range(5,250):
    nrz_bruit = add_noise(PSK,sigma/100)
    nrz_bruit_demoduler = demodulation_ASK(nrz_bruit, F_MESSAGE)-1
    sequence_nrz_bruit = construire_sequence_bit(nrz_bruit)
    sequence_nrz_demoduler_2 = construire_sequence_bit(nrz_bruit_demoduler)
    v.append(compute_BER(sequence_nrz_demoduler_2[1],sequence_nrz_bruit[1]))
plt.title('comparaison')
plt.plot(k)
plt.plot(v)
plt.ylim(10**-6,0)
plt.xlim(0,5)
plt.show()

for sigma in range(5,250):
    nrz_bruit = add_noise(PSK,sigma/100)
    nrz_bruit_demoduler = demodulation_ASK(nrz_bruit, F_MESSAGE)-1
    sequence_nrz_bruit = construire_sequence_bit(nrz_bruit)
    sequence_nrz_demoduler_2 = construire_sequence_bit(nrz_bruit_demoduler)
    v.append(compute_BER(sequence_nrz_demoduler_2[1],sequence_nrz_bruit[1]))
plt.title('comparaison')
plt.plot(k)
plt.plot(v)
plt.ylim(10**-6,0)
plt.xlim(0,5)
plt.show()