import os
from time import time
from numpy import argmin, mean, require, empty
from numpy.random import random
from causality.CmCmGenerator import CmCmGenerator

__author__ = 'pheodor'


class Orchestrator():

    def __init__(self):
        self.generator = CmCmGenerator()

    def pi_from_sigma_k(self, parametr):
        # Создадим массив коэффициентов
        koeff = [parametr, 2.13, 0.1, 3.77]
        # Переберем силу связи в разумных пределах
        Nsigma = 100
        Nk = 12
        dsigma = 0.01
        dk = 0.01
        mPI_s_k = empty((Nk, Nsigma))
        mEself_s_k = empty((Nk, Nsigma))
        mEjoin_s_k = empty((Nk, Nsigma))
        # Создаём папку
        catalogdata = "./Data/Cm-to-Cm_0.1_3.77_II/"+"%4.2f" % parametr
        if not os.path.exists(catalogdata):
            os.makedirs(catalogdata)
        for bf in xrange(Nk):
            # Сгенерируем тестовые сигналы
            (radX, radY) = self.generator.generate(10000, 0.1, 0.15, koeff, bf*dk)
            # Зашумим сигналы
            for sigma in xrange(Nsigma):
                tmp = time()
                # Усредним занчения PI по сотне точек
                mPI_average = []; mEself_average = []; mEjoin_average = []
                for ll in xrange(100):
                    # Добавляем шум:
                    ssigma = (sigma+1)*dsigma
                    radXnoise = radX + random.normal(0, ssigma, len(radX))
                    radYnoise = radY + random.normal(0, ssigma, len(radY))
                    # Рассчитаем ошибки аппроксимации и PI:
                    fitS = fit_self(radXnoise, 3)[1]/(len(radX)-1)
                    fitJ = fit_join(radXnoise, radYnoise, 3)[1]/(len(radX)-1)
                    if len(fitJ)>0 and len(fitS)>0:
                        mPI_average.append(1 - fitJ[0]/fitS[0])
                        mEself_average.append(fitS[0])
                        mEjoin_average.append(fitJ[0])
                # Подсчитаем среднее значение PI для данного сигма:
                if mPI_average:
                    mPI_s_k[bf, sigma] = mean(require(mPI_average))
                else:
                    mPI_s_k[bf, sigma] = 0
                if mEself_average:
                    mEself_s_k[bf, sigma] = mean(require(mEself_average))
                else:
                    mEself_s_k[bf, sigma] = 0
                if mEjoin_average:
                    mEjoin_s_k[bf, sigma] = mean(require(mEjoin_average))
                else:
                    mEjoin_s_k[bf, sigma] = 0
                print "sigma=%4.2f,  k=%4.2f,  time=%6.3f"%(ssigma, bf*dk, time()-tmp)
            # Сохраним полученный массив в файл:
            f = open(os.path.join(catalogdata, "PI-from-sigma"+"%4.2f"%(bf*dk)+".txt"), 'w')
            for ii in xrange(Nsigma):
                f.write("%4.2f"%((ii+1)*dsigma) + ' ' + "%14.10f"%mPI_s_k[bf, ii] + '\n')
            f.close()
            # Сохраним полученный массив в файл:
            f2 = open(os.path.join(catalogdata, "Eself-from-sigma"+"%4.2f"%(bf*dk)+".txt"), 'w')
            for ii in xrange(Nsigma):
                f2.write("%4.2f"%((ii+1)*dsigma) + ' ' + "%14.10f"%mEself_s_k[bf, ii] + '\n')
            f2.close()
            # Сохраним полученный массив в файл:
            f3 = open(os.path.join(catalogdata, "Ejoin-from-sigma"+"%4.2f"%(bf*dk)+".txt"), 'w')
            for ii in xrange(Nsigma):
                f3.write("%4.2f"%((ii+1)*dsigma) + ' ' + "%14.10f"%mEjoin_s_k[bf, ii] + '\n')
            f3.close()
            # Записываем чистый ряд Y:
            with open(os.path.join(catalogdata, "radY_"+"%4.2f"%(bf)+".txt"), 'w') as fy:
                for ii in xrange(len(radY)-1):
                    fy.write("%12.8f %12.8f\n"%(radY[ii], radY[ii+1]))
            # Записываем чистый ряд X при не нулевой силе связи:
            with open(os.path.join(catalogdata, "radX_"+"%4.2f"%(bf)+".txt"), 'w') as fx:
                for ii in xrange(len(radX)-1):
                    fx.write("%12.8f %12.8f\n"%(radX[ii], radX[ii+1]))
        sk = argmin(mPI_s_k[:, :], axis=1)*dsigma
        fsk = open(os.path.join(catalogdata, "sk.txt"), 'w')
        for s in xrange(Nk):
            fsk.write("%12.8f %12.8f"%(dk*s, sk[s])+"\n")
        fsk.close()
        return True