# coding=utf-8
import os
from time import time
from numpy import argmin, mean, require, empty, random
from causality.CmCmGenerator import CmCmGenerator
from causality.CommonUtils import CommonUtils

__author__ = 'pheodor'


class Orchestrator():

    def __init__(self, defaultsKoeff, nSigma, nK):
        self.koeffArray = defaultsKoeff
        self.generator = CmCmGenerator(10000, 0.1, 0.15, self.koeffArray)
        self.nSigma = nSigma
        self.nk = nK
        self.dk = self.dSigma = 0.01
        self.pathToFolder = "./Data/Cm-to-Cm_" + self.koeffArray[1] + "_" + self.koeffArray[2] + "_II/"

    def pi_from_sigma_k(self, parameter):
        # Создадим массив коэффициентов
        mPI_s_k = mEself_s_k = mEjoin_s_k = empty((self.nk, self.nSigma))
        # Создаём папку
        if not os.path.exists(self.pathToFolder):
            os.makedirs(self.pathToFolder)
        for cf in xrange(self.nk):
            # Сгенерируем тестовые сигналы
            (radX, radY) = self.generator.generate(cf * self.dk)
            # Зашумим сигналы
            for sigma in xrange(self.nSigma):
                startTime = time()
                # Усредним занчения PI по сотне точек
                mPI_average = mEself_average = mEjoin_average = []
                for ll in xrange(100):
                    # Добавляем шум:
                    ssigma = (sigma + 1) * self.dSigma
                    radXNoise = radX + random.normal(0, ssigma, len(radX))
                    radYNoise = radY + random.normal(0, ssigma, len(radY))
                    # Рассчитаем ошибки аппроксимации и PI:
                    fitS = CommonUtils.fit_self(radXNoise, 3)[1]/(len(radX)-1)
                    fitJ = CommonUtils.fit_join(radXNoise, radYNoise, 3)[1]/(len(radX)-1)
                    if len(fitJ) > 0 and len(fitS) > 0:
                        mPI_average.append(1 - fitJ[0]/fitS[0])
                        mEself_average.append(fitS[0])
                        mEjoin_average.append(fitJ[0])
                # Подсчитаем среднее значение PI для данного сигма:
                if mPI_average:
                    mPI_s_k[cf, sigma] = mean(require(mPI_average))
                else:
                    mPI_s_k[cf, sigma] = 0
                if mEself_average:
                    mEself_s_k[cf, sigma] = mean(require(mEself_average))
                else:
                    mEself_s_k[cf, sigma] = 0
                if mEjoin_average:
                    mEjoin_s_k[cf, sigma] = mean(require(mEjoin_average))
                else:
                    mEjoin_s_k[cf, sigma] = 0
                print "sigma=%4.2f,  k=%4.2f,  time=%6.3f" % (sigma, cf*self.dk, time()-startTime)
            # Сохраним полученный массив в файл:
            f = open(os.path.join(self.pathToFolder, "PI-from-sigma"+"%4.2f" % (cf*self.dk)+".txt"), 'w')
            for ii in xrange(self.nSigma):
                f.write("%4.2f" % ((ii+1)*self.dSigma) + ' ' + "%14.10f" % mPI_s_k[cf, ii] + '\n')
            f.close()
            # Сохраним полученный массив в файл:
            f2 = open(os.path.join(self.pathToFolder, "Eself-from-sigma"+"%4.2f" % (cf*self.dk)+".txt"), 'w')
            for ii in xrange(self.nSigma):
                f2.write("%4.2f" % ((ii+1)*self.dSigma) + ' ' + "%14.10f" % mEself_s_k[cf, ii] + '\n')
            f2.close()
            # Сохраним полученный массив в файл:
            f3 = open(os.path.join(self.pathToFolder, "Ejoin-from-sigma"+"%4.2f" % (cf*self.dk)+".txt"), 'w')
            for ii in xrange(self.nSigma):
                f3.write("%4.2f" % ((ii+1)*self.dSigma) + ' ' + "%14.10f" % mEjoin_s_k[cf, ii] + '\n')
            f3.close()
            # Записываем чистый ряд Y:
            with open(os.path.join(self.pathToFolder, "radY_"+"%4.2f" % (cf)+".txt"), 'w') as fy:
                for ii in xrange(len(radY)-1):
                    fy.write("%12.8f %12.8f\n" % (radY[ii], radY[ii+1]))
            # Записываем чистый ряд X при не нулевой силе связи:
            with open(os.path.join(self.pathToFolder, "radX_"+"%4.2f" % (cf)+".txt"), 'w') as fx:
                for ii in xrange(len(radX)-1):
                    fx.write("%12.8f %12.8f\n" % (radX[ii], radX[ii+1]))
        sk = argmin(mPI_s_k[:, :], axis=1) * self.dSigma
        fsk = open(os.path.join(self.pathToFolder, "sk.txt"), 'w')
        for s in xrange(self.nk):
            fsk.write("%12.8f %12.8f" % (self.dk*s, sk[s])+"\n")
        fsk.close()
        return True