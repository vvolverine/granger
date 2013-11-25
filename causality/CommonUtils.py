# coding=utf-8
from numpy import require, empty
from numpy.linalg import linalg

__author__ = 'pheodor'


class CommonUtils():

    def _polyn(self, P, d):
        """ Для всех базисных функций многомерного полинома кроме нулевой
            рассчитывает индексы родительской функции и элемента вектора состояния,
            на который её нужно домножить. Обходит все функции с помощью дерева.
            P - порядок полинома, d - размерность.
        """
        # Списки, содержащие номера родительских функций и номера домножаемых компонент:
        coel = []
        coej = []

        def polrep(ii, k, q):
            """ Рекурсивная функция, добавляющая новые индексы.
                k - минимальная величина индекса,
                ii - текущий номер родительской функции,
                q - полиномиальный порядок родительской функции.
            """
            for j0 in xrange(k, d):
                coel.append(ii)
                coej.append(j0)
                i = len(coel)
                if q < P:
                    polrep(i, j0, q+1)
        polrep(0, 0, 1)
        # Превращаем списки в массивы и возвращаем:
        return require(coel), require(coej)

    def fit_self(self, radX, rate):
        N = len(radX)
        v = empty((N-2, 2))
        v[:, 0] = radX[1:-1]
        v[:, 1] = radX[:-2]
        # Выясним номера и коэффициенты для полинома
        (c1, c2) = self._polyn(rate, v.shape[1])
        f = empty((N-2, len(c1)+1))
        f[:, 0] = 1
        for i in xrange(len(c1)):
            f[:, i+1] = f[:, c1[i]]*v[:, c2[i]]
        return linalg.lstsq(f, radX[2:])

    def fit_join(self, radX, radY, rate):
        N = len(radX)
        v = empty((N-2, 3))
        v[:, 0] = radX[1:-1]
        v[:, 1] = radX[:-2]
        v[:, 2] = radY[1:-1]
        # Выясним номера и коэффициенты для полинома
        (c1, c2) = self._polyn(rate, v.shape[1])
        f = empty((N-2, len(c1)+1))
        f[:, 0] = 1
        for i in xrange(len(c1)):
            f[:, i+1] = f[:, c1[i]]*v[:, c2[i]]
        return linalg.lstsq(f, radX[2:])