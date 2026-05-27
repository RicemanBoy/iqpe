from classes import *

def gen_data(name):
    p = np.linspace(0,0.007,10)
    y, y_qec, ys = [], [], []
    err, err_qec, y1s = [], [], []

    for r in p:
        ok, errr, y_list = avg15_repcode(3, 3, r, qec = False, k = 1, bias = 0.9)
        ok1, errr1, y1_list = avg15_repcode(3, 3, r, qec = True, k = 1, bias = 0.9)
        # ok, errr = avg15("steane", 3, noise=r, qec=True,k=1)
        y.append(ok), err.append(errr)
        # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
        y_qec.append(ok1), err_qec.append(errr1)
        # y = avg14_coin_success("steane", 3, r, True, k=1)
        ys, y1s = ys + y_list, y1s + y1_list

    data = np.array((p, y, y_qec, err, err_qec))
    #data = np.array((p, y, err))

    max_len = max(
        data.shape[1],
        len(ys),
        len(y1s)
    )

    def pad_row(row, length):
        padded = np.full(length, np.nan)
        padded[:len(row)] = row
        return padded

    padded_old = np.array([
        pad_row(row, max_len) for row in data
    ])

    row1 = pad_row(ys, max_len)
    row2 = pad_row(y1s, max_len)

    final_data = np.vstack([padded_old, row1, row2])

    np.savetxt("Repcode_d3_bias0_9{}.txt".format(name), final_data, delimiter=",")



# def gen_data(name):
#     p = np.linspace(0,0.007,10)
#     y, y_qec = [],[]
#     err, err_qec = [], []

#     for r in p:
#         ok, errr = avg15_repcode(7, 3, r, qec = False, k = 1, bias = -0.9)
#         ok1, errr1 = avg15_repcode(7, 3, r, qec = True, k = 1, bias = -0.9)
#         # ok, errr = avg15("steane", 3, noise=r, qec=True,k=1)
#         y.append(ok), err.append(errr)
#         # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
#         y_qec.append(ok1), err_qec.append(errr1)
#         # y = avg14_coin_success("steane", 3, r, True, k=1)

#     data = np.array((p, y, y_qec, err, err_qec))
#     #data = np.array((p, y, err))
#     np.savetxt("Repcode_d7_biasm0_9{}.txt".format(name), data, delimiter=",")
