import steane_ftec as s

def gen_data(name):                           #code OG
    p = s.np.linspace(0.00,0.005,20)
    # p = [np.linspace(0,0.005,6)[2]]
    y, y_qec = [],[]
    err, err_qec = [], []

    for r in p:  
        y_list = s.avg7_ramsey("steane", 3, r, qec = False, k = 1)    
        y.append(s.np.mean(y_list)), err.append(s.np.std(y_list))
        y1_list = s.avg7_ramsey("steane", 3, r, qec = True, k = 1, post = True) 
        y_qec.append(s.np.mean(y1_list)), err_qec.append(s.np.std(y1_list))

    data = s.np.array((p, y, y_qec, err, err_qec))
    #data = np.array((p, y, err))
    s.np.savetxt("steane_OG_ftqec+PS{}.txt".format(name), data, delimiter=",")
