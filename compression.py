#!/usr/bin/env python
# coding=utf-8

from random import randint, choice


ch = "ctgtcgtcgcaatgatacgtacagtctcgcaatga"


def construire_chaine(nb=1000):
    ch = ""
    lettres = ["a", "c", "g", "t"]
    for i in range(nb):
        ch += choice(lettres)
    return ch


def construire_sequences(nb=100):
    ch = ""
    lettres = ["a", "c", "g", "t"]
    for i in range(nb):
        repetitions = randint(5, 25)
        ch += choice(lettres) * repetitions
    return ch


def compression_RLE(msg):
    dst = ""
    N = len(msg)
    tmp = [1, msg[0]]
    for i in range(1,N):
        if msg[i] == tmp[1] :
            tmp[0]= tmp[0] + 1

        else:
            dst = dst+str(tmp[0])+str(tmp[1])
            tmp = [1, msg[i]]
    dst = dst+str(tmp[0])+str(tmp[1])
    return dst

def compression_RLC(msg, marqueur='#'):
    dst = ""
    N = len(msg)
    tmp = [1, msg[0]]
    for i in range(1, N):
        if msg[i] == tmp[1]:
            tmp[0] = tmp[0] + 1

        else:
            if tmp[0] == 1:
                dst = dst + str(tmp[1])
            elif tmp[0] == 2:
                dst = dst + str(tmp[1]) + str(tmp[1])
            elif tmp[0] >= 3:
                dst = dst + "#" + str(tmp[0]) + str(tmp[1])
            tmp = [1, msg[i]]
    if tmp[0] == 1:
        dst = dst + str(tmp[1])
    elif tmp[0] == 2:
        dst = dst + str(tmp[1]) + str(tmp[1])
    elif tmp[0] >= 3:
        dst = dst + "#" + str(tmp[0]) + str(tmp[1])
    return dst


def compression_LZW(msg):
    pass


def calcul_taux_compression(src, dst):
    pass


if __name__ == "__main__":
    print(ch)
    print(compression_RLE(ch))
    print(compression_RLC(ch))
    # print(compression_LZW(ch))