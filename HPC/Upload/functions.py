# import pyfiles.steane_ftec as s
# import pyfiles.rotsurf_ftec as rsf
# import pyfiles.bigsteane as bst
# import pyfiles.bigsteane_ftec as bstf
import pyfiles.classes as c
import numpy as np

def gen_data(name):                           #code OG
    # p = [0.5*1e-3]
    # p = [10**(-i-1) for i in range(7)][1:]
    p = np.linspace(0.00,0.001,10)
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:  
        y_list = c.avg7_RZ_ramsey_htoff(3, 3, r, qec = False, k = 1, bias = -1e99, short=True)        
        y.append(np.mean(y_list)), err.append(np.std(y_list))
        y1_list = c.avg7_RZ_ramsey_htoff(3, 3, r, qec = True, k = 1, bias = -1e99, short=True)      
        y_qec.append(np.mean(y1_list)), err_qec.append(np.std(y1_list))

    # data = np.array((p, y, y_qec, err, err_qec))
    data = np.array((p, y, err))
    np.savetxt("d3_rz_short_shorqec{}.txt".format(name), data, delimiter=",")
