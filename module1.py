# Créé par theog, le 12/09/2021 en Python 3.7
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
