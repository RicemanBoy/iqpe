from classes import *

def gen_data(name):
    p = np.linspace(0,0.001,10)
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:
        ok, errr = avg15("rotsurf9", 3, noise=r, qec=False,k=1)
        y.append(ok), err.append(errr)
        # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
        # y_qec.append(ok1), err_qec.append(errr1)
        # y = avg14_coin_success("steane", 3, r, True, k=1)

    #data = np.array((p, y, y_qec, err, err_qec))
    data = np.array((p, y, err))
    np.savetxt("nFT_RotSurf_PS{}.txt".format(name), data, delimiter=",")