from classes import *

def gen_data(name):
    p = np.linspace(0,0.001,10)
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:
        ok, errr = avg15("steane", 3, noise=r, qec=False,k=1)
        y.append(ok), err.append(errr)
        # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
        # y_qec.append(ok1), err_qec.append(errr1)
        # y = avg14_coin_success("steane", 3, r, True, k=1)

    #data = np.array((p, y, y_qec, err, err_qec))
    data = np.array((p, y, err))
    np.savetxt("nFT_steane_PS{}.txt".format(name), data, delimiter=",")
    # dustin = RotSurf16q(2)
    # dustin.h(0)
    # dustin.x(1)

    # dustin.t(0)
    # dustin.t(1)
    # dustin.cnot(0,1)
    # dustin.tdg(1)
    # dustin.cnot(0,1)
    # dustin.sdg(0)
    # dustin.h(0)

    # dustin.readout(0,1000,0)

    # #data = np.array((p, y, y_qec, err, err_qec))
    # data = np.array((dustin.zeros, dustin.ones))
    # np.savetxt("rotsurf16_testt.txt".format(name), data, delimiter=",")