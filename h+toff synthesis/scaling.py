"""Gate count vs precision for rysynth. Run: python3 scaling.py"""
import math
import time
import statistics as st
import rysynth as R

ANGLES = [0.4321, 1.0, 2.2, 3.9, 5.31]          # generic angles, no symmetry
EXPS = list(range(1, 11))                        # eps = 10^-1 .. 10^-12

print("  eps      k     total    CCX      H      X    err/eps   sec")
print("  " + "-" * 62)
rows = []
for e in EXPS:
    eps = 10.0 ** (-e)
    ks, tot, ccx, hs, xs, ratios = [], [], [], [], [], []
    t0 = time.time()
    for th in ANGLES:
        r = R.synthesize_ry(th, eps)
        c = r['counts']
        err = R.clean_ancilla_error(r['gates'], th)
        assert err <= eps * 1.001, (th, eps, err)
        ks.append(r['k'])
        tot.append(c['total'])
        ccx.append(c.get('CCX', 0))
        hs.append(c.get('H', 0))
        xs.append(c.get('X', 0))
        ratios.append(err / eps)
    dt = time.time() - t0
    rows.append((e, st.mean(ks), st.mean(tot), st.mean(ccx), st.mean(hs),
                 st.mean(xs), max(ratios), dt))
    print("  1e-%-2d  %5.1f  %7.1f %6.1f %6.1f %6.1f   %6.3f  %5.2f"
          % (e, st.mean(ks), st.mean(tot), st.mean(ccx), st.mean(hs),
             st.mean(xs), max(ratios), dt))


def slope(xs, ys):
    mx, my = st.mean(xs), st.mean(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / \
           sum((x - mx) ** 2 for x in xs)


L = [e * math.log2(10) for e, *_ in rows]        # log2(1/eps)
print("\n  least-squares fits against log2(1/eps):")
print("    k        = %.3f * log2(1/eps) %+.2f   (Ross-Selinger optimum: 3.0)"
      % (slope(L, [r[1] for r in rows]),
         st.mean([r[1] for r in rows]) - slope(L, [r[1] for r in rows]) * st.mean(L)))
for name, idx in (("total", 2), ("CCX", 3), ("H", 4)):
    s = slope(L, [r[idx] for r in rows])
    b = st.mean([r[idx] for r in rows]) - s * st.mean(L)
    print("    %-8s = %.3f * log2(1/eps) %+.2f" % (name, s, b))

print("\n  Solovay-Kitaev comparison (c=3.97, same asymptotic constant "
      "convention):")
for e in (4, 8, 12):
    L1 = e * math.log2(10)
    print("    eps=1e-%-2d   this work ~%5.0f gates   SK ~%.1e gates"
          % (e, slope(L, [r[2] for r in rows]) * L1, L1 ** 3.97))
