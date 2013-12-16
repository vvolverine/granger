from multiprocessing import Pool, cpu_count
import argparse
from numpy import arange
from causality.Orchestrator import Orchestrator

__author__ = 'pheodor'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='describe orchestrator configuration')
    parser.add_argument("--nSigma", help='define maximum sigma value(in points)', default=100, type=int)
    parser.add_argument("--nK", help='define maximum coupling force(in points)', default=12, type=int)
    parser.add_argument("-d", help='define delimiter of points value', default=.01, type=float)
    parser.add_argument("--folder", help='set base dir for text results')
    parser.add_argument("--length", help='set the time range length', default=10000, type=int)
    parser.add_argument("--starts", help='set initial value of X and Y time ranges', nargs='*', type=float)
    parser.add_argument("--weights", help='define default value in coefficient array', nargs='*', type=float)
    args = parser.parse_args()
    orchestrator = Orchestrator(args.weights, args.starts, args.length, args.nK, args.nSigma, args.d, args.folder)
    pool = Pool(processes=cpu_count())
    pool.map(orchestrator.pi_from_sigma_k, arange(0.1, 0.3, 0.05))