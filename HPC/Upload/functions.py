from classes import *

def gen_data(name):
    p = np.linspace(0,0.005,10)
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:
        ok, errr = avg15_coin("not steane", 3, noise=r, qec=True,k=1)
        y.append(ok), err.append(errr)
        # ok1, errr1 = avg15_coin(3, 15, noise=r, err=True, k=1)
        # y_qec.append(ok1), err_qec.append(errr1)

    #data = np.array((p, y, y_qec, err, err_qec))
    data = np.array((p, y, err))
    np.savetxt("RotSurf_test_more_qec{}.txt".format(name), data, delimiter=",")