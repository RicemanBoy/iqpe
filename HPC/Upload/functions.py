from classes import *

def gen_data(name):
    p = [np.linspace(0,0.001,10)[5]]
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:
        ok, errr = avg15_coin("steane", 3, noise=r, qec=True,k=1)
        y.append(ok), err.append(errr)
        # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
        # y_qec.append(ok1), err_qec.append(errr1)
        # y = avg14_coin_success("steane", 3, r, True, k=1)

    #data = np.array((p, y, y_qec, err, err_qec))
    data = np.array((p, y, err))
    np.savetxt("Steane_noPS_qec4_noflags_5{}.txt".format(name), data, delimiter=",")