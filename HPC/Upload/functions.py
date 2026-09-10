import pyfiles.steane_ftec as s
import pyfiles.rotsurf_ftec as rsf
import pyfiles.bigsteane as bst
import pyfiles.bigsteane_ftec as bstf

def gen_data(name):                           #code OG
    p = s.np.linspace(0.00,0.005,10)
    # p = s.np.linspace(0.0015, 0.003, 3)
    # p = [np.linspace(0,0.005,6)[2]]
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:  
        y_list = bstf.avg7_ramsey("bigsteane", 3, r, qec = False, k = 1)    
        y.append(s.np.mean(y_list)), err.append(s.np.std(y_list))
        y1_list = bstf.avg7_ramsey("bigsteane", 3, r, qec = True, k = 1, post = False) 
        y_qec.append(s.np.mean(y1_list)), err_qec.append(s.np.std(y1_list))

    data = s.np.array((p, y, y_qec, err, err_qec))
    #data = np.array((p, y, err))
    s.np.savetxt("bigsteane_ftqec_lots{}.txt".format(name), data, delimiter=",")
