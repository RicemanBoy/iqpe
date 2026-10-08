# import pyfiles.steane_ftec as s
# import pyfiles.rotsurf_ftec as rsf
# import pyfiles.bigsteane as bst
# import pyfiles.bigsteane_ftec as bstf
import pyfiles.classes as c
import numpy as np

def gen_data(name):                           #code OG
    p = [np.linspace(0.00,0.001,10)[1]]
    # p = s.np.linspace(0.0015, 0.003, 3)
    # p = [np.linspace(0,0.005,6)[2]]
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:  
        y_list = c.avg7_ramsey_htoff(5, 3, r, qec = True, k = 1, bias = -1e99)        
        y.append(np.mean(y_list)), err.append(np.std(y_list))
        # y1_list = c.avg7_ramsey_htoff(3, 3, r, qec = True, k = 1, bias = -1e99)      
        # y_qec.append(np.mean(y1_list)), err_qec.append(np.std(y1_list))

    # data = np.array((p, y, y_qec, err, err_qec))
    data = np.array((p, y, err))
    np.savetxt("d5_htoffcxbeam_infbias_qec_0{}.txt".format(name), data, delimiter=",")
