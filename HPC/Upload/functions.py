from classes import *


def gen_data(name):                   #code mit npz
    p = np.linspace(0.00,0.005,6)
    bias = [-1e4, -100, -10, 0, 10, 100, 1e4]

    y = np.zeros((len(bias), len(p), 7))           #7 Winkel !!!
    y_qec = np.zeros((len(bias), len(p), 7))       

    for i, b in enumerate(bias):
        for j, r in enumerate(p):
            _, _, y_list = avg7_repcode("z", 3, 3, r, qec = False, post= False, k = 1, bias = b)
            _, _, y_qec_list = avg7_repcode("z", 3, 3, r, qec = True, post= False, k = 1, bias = b)

            y[i, j, :] = y_list
            y_qec[i, j, :] = y_qec_list

    np.savez(
    f"Phasecode_d3_exactangles_nftqec{name}.npz",
    p=p,
    bias=bias,
    y=y,
    y_qec=y_qec
    )

# def gen_data(name):                       #code1
#     p = np.linspace(0,0.007,10)
#     bias = np.linspace(-10, 10, 10)
#     y, y_qec = [], []

#     for ö in bias:
#         for r in p:
#             _, _, y_list = avg15_repcode(5, 3, r, qec = False, k = 1, bias = ö)
#             _, _, y_qec_list = avg15_repcode(5, 3, r, qec = True, k = 1, bias = ö)
#             y, y_qec = y + y_list, y_qec + y_qec_list

#     data = np.array((p, bias))
#     #data = np.array((p, y, err))

#     max_len = max(
#         data.shape[1],
#         len(y),
#         len(y_qec)
#     )

#     def pad_row(row, length):
#         padded = np.full(length, np.nan)
#         padded[:len(row)] = row
#         return padded

#     padded_old = np.array([
#         pad_row(row, max_len) for row in data
#     ])

#     row1 = pad_row(y, max_len)
#     row2 = pad_row(y_qec, max_len)

#     final_data = np.vstack([padded_old, row1, row2])

#     np.savetxt("Repcode_d5{}.txt".format(name), final_data, delimiter=",")

# def gen_data(name):                           #code 0 
#     p = np.linspace(0,0.007,10)
#     y, y_qec, ys = [], [], []
#     err, err_qec, y1s = [], [], []

#     for r in p:
#         ok, errr, y_list = avg15_repcode(3, 3, r, qec = False, k = 1, bias = 0.9)
#         ok1, errr1, y1_list = avg15_repcode(3, 3, r, qec = True, k = 1, bias = 0.9)
#         # ok, errr = avg15("steane", 3, noise=r, qec=True,k=1)
#         y.append(ok), err.append(errr)
#         # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
#         y_qec.append(ok1), err_qec.append(errr1)
#         # y = avg14_coin_success("steane", 3, r, True, k=1)
#         ys, y1s = ys + y_list, y1s + y1_list

#     data = np.array((p, y, y_qec, err, err_qec))
#     #data = np.array((p, y, err))

#     max_len = max(
#         data.shape[1],
#         len(ys),
#         len(y1s)
#     )

#     def pad_row(row, length):
#         padded = np.full(length, np.nan)
#         padded[:len(row)] = row
#         return padded

#     padded_old = np.array([
#         pad_row(row, max_len) for row in data
#     ])

#     row1 = pad_row(ys, max_len)
#     row2 = pad_row(y1s, max_len)

#     final_data = np.vstack([padded_old, row1, row2])

#     np.savetxt("Repcode_d3_bias0_9{}.txt".format(name), final_data, delimiter=",")



# def gen_data(name):                           #code OG
#     p = np.linspace(0,0.001,10)[4:8]
#     y, y_qec = [],[]
#     err, err_qec = [], []

#     for r in p:
#         ok, errr = avg15_coin("steane", 3, r, qec = True, k = 1)
#         # ok, errr = avg15("steane", 3, noise=r, qec=True,k=1)
#         y.append(ok), err.append(errr)
#         # ok1, errr1 = avg15_coin("notsteane", 3, noise=r, qec=True, k=1)
#         # y = avg14_coin_success("steane", 3, r, True, k=1)

#     data = np.array((p, y, err))
#     #data = np.array((p, y, err))
#     np.savetxt("FTSteane_qec_4_7_neww{}.txt".format(name), data, delimiter=",")
