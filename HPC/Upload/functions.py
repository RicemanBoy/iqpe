from classes import *

def gen_data(name):
    p = np.linspace(0,0.007,10)
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:
        ok, errr = avg15_repcode(7, 3, r, qec = False, k = 1, bias = -0.9)
        ok1, errr1 = avg15_repcode(7, 3, r, qec = True, k = 1, bias = -0.9)
        # ok, errr = avg15("steane", 3, noise=r, qec=True,k=1)
        y.append(ok), err.append(errr)
        # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
        y_qec.append(ok1), err_qec.append(errr1)
        # y = avg14_coin_success("steane", 3, r, True, k=1)

    data = np.array((p, y, y_qec, err, err_qec))
    #data = np.array((p, y, err))
    np.savetxt("Repcode_d7_biasm0_9{}.txt".format(name), data, delimiter=",")