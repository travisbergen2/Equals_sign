#!/usr/bin/env python3
"""
rh_calculator_twin.py -- Python twin of Room IX "The Equals Sign" (the RH calculator), v0.2.
Computes, with committed algorithms, every reference value the page self-checks against; verifies the page's
double-precision algorithms (Euler-Maclaurin zeta, Lanczos log-gamma) against mpmath; ports Paper 14's committed
A(m) = POLE(m) - PRIME(m) + ARCH(m) construction verbatim from the E-HECKE-1 census code (ehecke_census.py,
2026-08-26) and reproduces the paper's N-14-1 anchors; writes rh_calculator_embed.json =
{zeros: first 10,000 ordinates (Odlyzko, 9 dp), ref: {...}}.
Inputs: zeros1.txt (Odlyzko first 100,000 zeros; sha256 3436c916a7878261ac183fd7b9448c9a4736b8bbccf1356874a6ce1788541632)
Nothing here bears on the truth of RH: every value is a finite instance of an equivalent form.
"""
import math, cmath, json, os, time
import numpy as np
from mpmath import mp, mpf, mpc as MPC, zeta as mpzeta, loggamma as mploggamma, siegeltheta, psi as mpsi, pi as MPPI, log as mplog, catalan

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
ZPATH = next(p for p in [os.path.join(HERE, 'zeros1.txt'), 'zeros1.txt', '/agent/workspace/zeros1.txt'] if os.path.exists(p))
zeros = np.loadtxt(ZPATH); assert zeros.shape == (100000,)
NZ = 10000
Z = zeros[:NZ]; T = float(Z[-1]); T100k = float(zeros[-1])
ref = {'NZ': NZ, 'T': T}
TWO_PI = 2*math.pi
GAMMA_E = 0.57721566490153286061

# ---------- 1. a zero's dictionary ----------
def rim_theta(g): return 2*math.atan(1/(2*g))
ref['rim_theta_deg_zero1'] = math.degrees(rim_theta(Z[0]))
def tail_inv_gamma2(Tt): return (math.log(Tt/TWO_PI) + 1)/(TWO_PI*Tt)
def li_lambda(zs, Tt, n):
    th = 2*np.arctan(1/(2*zs)); return float(np.sum(4*np.sin(n*th/2)**2) + n*n*tail_inv_gamma2(Tt))
ref['li_lambda_n1_10'] = [li_lambda(Z, T, n) for n in range(1, 11)]
ref['li_lambda_100k_n1_10'] = [li_lambda(zeros, T100k, n) for n in range(1, 11)]
lam1_exact = 1 + GAMMA_E/2 - 0.5*math.log(4*math.pi); ref['li_lambda1_exact'] = lam1_exact
assert abs(ref['li_lambda_100k_n1_10'][0] - lam1_exact) < 2e-8
assert abs(ref['li_lambda_n1_10'][0] - lam1_exact) < 1e-7, ref['li_lambda_n1_10'][0]
def moment_tail(n, Tt):   # int_T^inf gamma^{-2n} dN(gamma), Riemann-von Mangoldt density
    return Tt**(1-2*n)/(TWO_PI*(2*n-1)) * (math.log(Tt/TWO_PI) + 1.0/(2*n-1))
ref['p2_plus_tail'] = float(np.sum(Z**-2.0) + moment_tail(1, T))
ref['p2_exact'] = 0.023104993115419
assert abs(ref['p2_plus_tail'] - ref['p2_exact']) < 1e-8, ref['p2_plus_tail']
ref['count_zeros_le_100'] = int(np.sum(Z <= 100)); assert ref['count_zeros_le_100'] == 29
ref['count_zeros_le_1000'] = int(np.sum(Z <= 1000)); assert ref['count_zeros_le_1000'] == 649
def N_main(Tt): return (Tt/TWO_PI)*math.log(Tt/(TWO_PI*math.e)) + 7/8
ref['N_main_100'] = N_main(100.0)
ref['unfolded_gap_1_2'] = float((Z[1]-Z[0])*math.log(Z[0]/TWO_PI)/TWO_PI)
def ef_term(x, beta, g):
    rho = complex(beta, g); return float(-2*(x**rho/rho).real)
ref['ef_term_zero1_x1000'] = ef_term(1000.0, 0.5, Z[0])

# ---------- 2. Davenport-Heilbronn traitor zero ----------
DH_B, DH_G = 0.808517182456637, 85.6993484853776
rho = complex(DH_B, DH_G); mirror = 1 - rho.conjugate()
ref['dh_beta'] = DH_B; ref['dh_gamma'] = DH_G
ref['dh_recip_circle_defect'] = (1-2*DH_B)/abs(rho)**2
ref['dh_pair_forced_sign'] = (2*DH_B-1)*(1/abs(mirror)**2 - 1/abs(rho)**2); assert ref['dh_pair_forced_sign'] > 0
quad = [rho, rho.conjugate(), mirror, mirror.conjugate()]
def li_quad(n): return sum(1 - (1 - 1/r)**n for r in quad).real
first_neg = next(n for n in range(1, 3000) if li_quad(n) < 0); ref['dh_li_first_negative_n'] = first_neg; assert first_neg == 537
ref['dh_li_contrib_n80'] = li_quad(80)
ref['dh_disc_radius'] = abs(1 - 1/rho)
ref['dh_weil_leg_deg_a0p1'] = math.degrees(0.01*DH_G*(2*DH_B-1))
ref['dh_visibility_a_min'] = math.sqrt(math.pi/(4*(DH_B-0.5)*DH_G))
ref['dh_N_below'] = int(np.sum(Z < DH_G))   # true zeros below the traitor's height (23)

# ---------- 3. integers ----------
def factor(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0)+1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0)+1
    return f
def sopfr(n): return sum(p*e for p, e in factor(n).items())
def sigma(n):
    s = 1
    for p, e in factor(n).items(): s *= (p**(e+1)-1)//(p-1)
    return s
def gaussian(p):
    if p == 2: return (1, 1, 45.0, -1.0)
    if p % 4 == 3: return None
    for b in range(1, int(math.isqrt(p))+1):
        a2 = p - b*b; a = math.isqrt(a2)
        if a*a == a2 and a >= b: return (a, b, math.degrees(math.atan2(b, a)), (complex(a, b)**4).real/p**2)
def eisenstein(p):
    if p == 3: return (2, 1, 30.0, -1.0)
    if p % 3 == 2: return None
    w = complex(-0.5, math.sqrt(3)/2)
    for b in range(1, int(math.isqrt(p))+2):
        for a in range(b, int(math.isqrt(p))+2):
            if a*a - a*b + b*b == p:
                z = a + b*w; th = math.degrees(cmath.phase(z)) % 60; th = min(th, 60-th)
                return (a, b, th, (z**6).real/p**3)
g13 = gaussian(13); assert (g13[0], g13[1]) == (3, 2) and abs(g13[3] + 119/169) < 1e-12
e7 = eisenstein(7); assert abs(e7[3] + 143/343) < 1e-9
e13 = eisenstein(13); assert abs(e13[3] - 253/2197) < 1e-9
ref['gauss_13'] = list(g13); ref['eis_7'] = list(e7); ref['eis_13'] = list(e13)
ref['sopfr_fixed_points_le_100'] = [n for n in range(2, 101) if sopfr(n) == n]
def harmonic(n): return sum(1.0/k for k in range(1, n+1))
def lagarias(n):
    H = harmonic(n); rhs = H + math.exp(H)*math.log(H); return {'sigma': sigma(n), 'H': H, 'rhs': rhs, 'margin': rhs - sigma(n)}
ref['lagarias_5040'] = lagarias(5040); assert ref['lagarias_5040']['margin'] > 0
def robin_margin(n): return math.exp(GAMMA_E)*n*math.log(math.log(n)) - sigma(n)
ref['robin_margin_10080'] = robin_margin(10080); assert ref['robin_margin_10080'] > 0

# ---------- 4. sieve to 2e6 ----------
NS = 2_000_000
spf = np.zeros(NS+1, dtype=np.int32)
for i in range(2, NS+1):
    if spf[i] == 0: spf[i:NS+1:i][spf[i:NS+1:i] == 0] = i
n_arr = np.arange(NS+1); is_prime = (spf == n_arr) & (n_arr >= 2)
Lam = np.zeros(NS+1); mu = np.zeros(NS+1, dtype=np.int8); Om = np.zeros(NS+1, dtype=np.int16); pp = np.zeros(NS+1, dtype=np.uint8); mu[1] = 1
for n in range(2, NS+1):
    p = spf[n]; m = n // p
    Om[n] = Om[m] + 1; mu[n] = 0 if spf[m] == p else -mu[m]
    pp[n] = 1 if (m == 1 or (pp[m] and spf[m] == p)) else 0
    if pp[n]: Lam[n] = math.log(p)
theta_c = np.cumsum(np.where(is_prime, np.log(np.maximum(n_arr, 1)), 0.0)); psi_c = np.cumsum(Lam); pi_c = np.cumsum(is_prime.astype(np.int64))
M_c = np.cumsum(mu.astype(np.int64)); L_c = np.cumsum(np.where(n_arr >= 1, (-1)**Om, 0).astype(np.int64))
chi4 = np.where(is_prime & (n_arr % 4 == 1), 1, 0) - np.where(is_prime & (n_arr % 4 == 3), 1, 0); race_c = np.cumsum(-chi4)
def li(x):
    if x <= 1: return float('-inf')
    lx = math.log(x); s = 0.0; fact = 1.0; inner = 0.0
    for n in range(1, 200):
        fact *= n
        if (n-1) % 2 == 0: inner += 1.0/n
        term = ((-1)**(n-1)) * lx**n / (fact * 2**(n-1)) * inner; s += term
        if abs(term) < 1e-17*abs(s) and n > 20: break
    return GAMMA_E + math.log(lx) + math.sqrt(x)*s
tab = {}
for k in range(2, 7):
    x = 10**k
    tab[str(x)] = dict(pi=int(pi_c[x]), theta=float(theta_c[x]), psi=float(psi_c[x]), M=int(M_c[x]), L=int(L_c[x]), race=int(race_c[x]), li=li(x))
ref['sieve_table'] = tab; ref['sieve_cap'] = NS
assert tab['1000000']['pi'] == 78498 and tab['1000000']['M'] == 212 and tab['1000000']['L'] == -530 and tab['1000000']['race'] == 147

# ---------- 5. explicit formula with the bundled zeros ----------
def psi_explicit(x, zs):
    s = x - math.log(TWO_PI) - 0.5*math.log(1 - x**-2.0)
    return float(s - np.sum(2*np.sqrt(x)*(0.5*np.cos(zs*math.log(x)) + zs*np.sin(zs*math.log(x)))/(0.25 + zs*zs)))
for x in (1000.5, 10000.5, 100000.5):
    ref[f'psi_explicit_x{int(x)}'] = psi_explicit(x, Z); ref[f'psi_true_x{int(x)}'] = float(psi_c[int(x)])
    assert abs(ref[f'psi_explicit_x{int(x)}'] - psi_c[int(x)]) < 8.0

# ---------- 6. function-field world ----------
def hasse(p, a, b):
    assert (4*a**3 + 27*b*b) % p != 0
    s = 0
    for x in range(p):
        r = (x*x*x + a*x + b) % p
        if r == 0: continue
        s += 1 if pow(r, (p-1)//2, p) == 1 else -1
    a_p = -s; return dict(a_p=a_p, points=p+1-a_p, abs_alpha_over_sqrt_p=abs(complex(a_p, math.sqrt(4*p-a_p*a_p))/2)/math.sqrt(p))
ref['hasse_1009'] = hasse(1009, 1, 1); ref['hasse_10007'] = hasse(10007, 1, 1)

# ---------- 7. a point s ----------
B2K = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330, 854513/138, -236364091/2730]
def zeta_em(s, K=12):
    N = max(30, int(abs(s.imag)) + 20)
    tot = sum(complex(n)**(-s) for n in range(1, N)) + N**(1-s)/(s-1) + 0.5*N**(-s)
    rising = s
    for k in range(1, K+1):
        tot += B2K[k-1]/math.factorial(2*k) * rising * N**(-s-2*k+1); rising *= (s+2*k-1)*(s+2*k)
    return tot
LANCZOS_C = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
def lgamma_c(z):
    if z.real < 0.5: return cmath.log(cmath.pi/cmath.sin(cmath.pi*z)) - lgamma_c(1-z)
    z -= 1; x = LANCZOS_C[0]
    for i in range(1, 9): x += LANCZOS_C[i]/(z+i)
    t = z + 7.5
    return 0.5*math.log(2*math.pi) + (z+0.5)*cmath.log(t) - t + cmath.log(x)
def xi_c(s): return s*(s-1)/2 * cmath.exp(-s/2*math.log(math.pi) + lgamma_c(s/2)) * zeta_em(s)
def rs_theta(t): return (lgamma_c(complex(0.25, t/2))).imag - t/2*math.log(math.pi)
mp.dps = 30
checks = {}
for name, s in [('s2', complex(2, 0)), ('rho1', complex(0.5, float(Z[0]))), ('s_0p75_50', complex(0.75, 50)), ('s_3_4', complex(3, 4)), ('s_half_100p5', complex(0.5, 100.5)), ('s_half_1000p5', complex(0.5, 1000.5))]:
    ze = zeta_em(s); zm = complex(mpzeta(MPC(s.real, s.imag)))
    checks[name] = dict(s=[s.real, s.imag], zeta_em=[ze.real, ze.imag], abs_err_vs_mpmath=abs(ze-zm)); assert abs(ze - zm) < 1e-10*max(1, abs(zm))
ref['zeta_checks'] = checks
for z in [complex(0.25, float(Z[0])/2), complex(1.5, 0), complex(0.25, 50)]:
    assert abs(lgamma_c(z) - complex(mploggamma(MPC(z.real, z.imag)))) < 1e-11*max(1, abs(lgamma_c(z)))
ref['xi_half'] = xi_c(complex(0.5, 0)).real; assert abs(2*ref['xi_half'] - 0.99424155637662) < 1e-11
ref['rs_theta_at_gamma1'] = rs_theta(float(Z[0])); assert abs(ref['rs_theta_at_gamma1'] - float(siegeltheta(mpf(Z[0])))) < 1e-9
ref['rs_theta_100'] = rs_theta(100.0)
Zt = cmath.exp(1j*rs_theta(100.5))*zeta_em(complex(0.5, 100.5)); ref['Z_100p5'] = Zt.real
def lagarias_sum(s, zs, Tt):
    tot = 0j
    for g in zs: tot += 1/(s - complex(0.5, g)) + 1/(s - complex(0.5, -g))
    return tot + 2*(s-0.5)*tail_inv_gamma2(Tt)
ls = lagarias_sum(complex(0.75, 50), Z, T); ls100k = lagarias_sum(complex(0.75, 50), zeros, T100k)
ref['lagarias_sum_0p75_50'] = [ls.real, ls.imag]; ref['lagarias_sum_0p75_50_100k'] = [ls100k.real, ls100k.imag]
assert abs(ls100k - complex(2.323231, -1.441188)) < 2e-5 and abs(ls - ls100k) < 1e-3
ref['lagarias_angle_deg'] = math.degrees(cmath.phase(ls))

# ---------- 8. LIVE HANKEL (TSF Theorem 4.1): q_m = sum gamma^{-2(m+1)} + tail, double precision ----------
def q_seq(zs, Tt, M, plant=None):
    q = []
    for m in range(M):
        n = m+1
        v = float(np.sum(zs**(-2.0*n))) + moment_tail(n, Tt)
        if plant is not None:   # plant = (gamma0, delta): replace the on-line zero gamma0 by the pair z = gamma0 -/+ i delta
            g0, dl = plant; v += -g0**(-2.0*n) + 2*(complex(g0, -dl)**(-2*n)).real
        q.append(v)
    return q
def hankel_lmin(q, k): return float(np.linalg.eigvalsh(np.array([[q[i+j] for j in range(k)] for i in range(k)]))[0])
KH = 10
qH = q_seq(Z, T, 2*KH-1)
ref['hankel_q'] = qH
ref['hankel_lmin_live'] = [hankel_lmin(qH, k) for k in range(1, KH+1)]
HANKEL_EXACT = [0.023104993115419, 8.43684e-8, 1.11461e-13, 5.41937e-20, 1.66617e-26, 4.12689e-33, 5.13634e-40]   # Newton route, 2026-09-27
ref['hankel_lmin_exact_k1_7'] = HANKEL_EXACT
for k in range(1, 8): assert abs(ref['hankel_lmin_live'][k-1]/HANKEL_EXACT[k-1] - 1) < 2e-5, (k, ref['hankel_lmin_live'][k-1], HANKEL_EXACT[k-1])
ref['hankel_noise_floor'] = 2.0**-52 * qH[0]
def corr_lmin(q, k):
    H = np.array([[q[i+j] for j in range(k)] for i in range(k)]); d = np.sqrt(np.diag(H)); return float(np.linalg.eigvalsh(H/np.outer(d, d))[0])
ref['hankel_corr_lmin'] = [corr_lmin(qH, k) for k in range(1, 8)]
qP = q_seq(Z, T, 2*KH-1, plant=(float(Z[0]), 0.1))
lp = [hankel_lmin(qP, k) for k in range(1, KH+1)]
ref['hankel_plant_zero1_beta0p6_lmin'] = lp
kstar = next(k for k in range(1, KH+1) if lp[k-1] < 0); ref['hankel_plant_zero1_beta0p6_kstar'] = kstar; assert kstar == 5, lp
ref['hankel_kstar_law_traitor'] = 1.83*ref['dh_N_below'] + 5.4

# ---------- 9. PAPER 14 A(m) -- ported verbatim from ehecke_census.py (E-HECKE-1, 2026-08-26) ----------
DELTA = 0.05
LOGPI = float(mplog(MPPI))
PSI_QUARTER = float(mpsi(0, mpf(1)/4)); PSI1_QUARTER = float(mpsi(1, mpf(1)/4))
assert abs(PSI_QUARTER - (-GAMMA_E - math.pi/2 - 3*math.log(2))) < 1e-14 and abs(PSI1_QUARTER - (math.pi**2 + 8*float(catalan))) < 1e-12
ref['psi_quarter'] = PSI_QUARTER; ref['psi1_quarter'] = PSI1_QUARTER; ref['log_pi'] = LOGPI
def sieve_bytes(n):
    s = bytearray([1])*(n+1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return s
def prime_powers(limit, s):
    out = []
    for p in range(2, limit+1):
        if s[p]:
            n = p; e = 1
            while n <= limit: out.append((n, p, e)); n *= p; e += 1
    out.sort(); return out
def tent_hits(n, delta, M):
    ln = math.log(n); m0 = int(ln/delta)
    for m in (m0-1, m0, m0+1, m0+2):
        if 0 <= m <= M:
            w = 1.0 - abs(ln - m*delta)/delta
            if w > 0: yield m, w
def prime_array_zeta(delta, M, cutoff, s):
    arr = np.zeros(M+1)
    for n, p, e in prime_powers(cutoff, s):
        c = math.log(p)/math.sqrt(n)
        for m, w in tent_hits(n, delta, M): arr[m] += c*w
    return arr
def S_series(a0, step, t, tol=1e-22):
    if t == 0: return PSI1_QUARTER/(step*step)      # psi'(a0/step)/step^2 with a0/step = 1/4
    s = 0.0; j = 0
    while True:
        a = a0 + step*j; term = math.exp(-a*t)/(a*a); s += term
        if term < tol and a*t > 1: break
        j += 1
        if j > 500000: break
    return s
def arch_array_zeta(delta, M):
    a0, step, const = 0.5, 2.0, PSI_QUARTER - LOGPI
    S = [S_series(a0, step, m*delta) for m in range(M+2)]
    arr = np.empty(M+1); arr[0] = const + (2.0/delta)*(S[0]-S[1])
    for m in range(1, M+1): arr[m] = -(1.0/delta)*(S[m-1] - 2*S[m] + S[m+1])
    return arr
def pole_array(delta, M):
    C = (8.0/delta)*(math.cosh(delta/2)-1.0); return np.array([2*C*math.cosh(m*delta/2) for m in range(M+1)])
def toe(arr, Mi):
    idx = np.abs(np.subtract.outer(np.arange(Mi), np.arange(Mi))); return arr[idx]
def mu_min(arr, Mi): return float(np.linalg.eigvalsh(toe(arr, Mi))[0])
def wall_scan(arr, delta, Wcap, undec=1e-8):
    for Mi in range(1, int(round(Wcap/delta))+1):
        mu = mu_min(arr, Mi)
        if mu < -undec: return Mi*delta, mu
    return None, None
MBIG = 300; CUT_BIG = 3445000
sbig = sieve_bytes(CUT_BIG)
POLE = pole_array(DELTA, MBIG); ARCH = arch_array_zeta(DELTA, MBIG); PRIME = prime_array_zeta(DELTA, MBIG, CUT_BIG, sbig)
A = POLE - PRIME + ARCH
ref['A_delta'] = DELTA; ref['A_cutoff'] = CUT_BIG
ref['A_first_12'] = [float(v) for v in A[:12]]; ref['POLE_first_6'] = [float(v) for v in POLE[:6]]; ref['ARCH_first_6'] = [float(v) for v in ARCH[:6]]; ref['PRIME_first_6'] = [float(v) for v in PRIME[:6]]
ref['A_lag300'] = dict(POLE=float(POLE[300]), PRIME=float(PRIME[300]), ARCH=float(ARCH[300]), A=float(A[300]))
assert abs(POLE[300] - 90.41) < 0.02 and abs(PRIME[300] - 90.59) < 0.02, (POLE[300], PRIME[300])   # Paper 14 Sec 4: +90.41 / -90.59 at lag 300
anchors = {}
targets = {2.5: 1.34e-3, 5.0: 1.85e-4, 7.5: 6.4e-5, 10.0: 4.4e-5, 12.5: 3.1e-5, 15.0: 2.2e-5}
census_mine = {2.5: 0.001335335532445962, 5.0: 0.0001848120699362628, 7.5: 6.382446636081962e-05, 10.0: 4.4007386876924974e-05, 12.5: 3.12447556013296e-05, 15.0: 2.211753658283158e-05}
for W, tgt in targets.items():
    mu = mu_min(A, int(round(W/DELTA))); anchors[str(W)] = mu
    assert abs(mu/census_mine[W] - 1) < 1e-9, (W, mu, census_mine[W])          # identical to the E-HECKE-1 recomputation
    assert abs(mu - tgt) <= 2*10.0**(math.floor(math.log10(abs(tgt)))-1)         # Paper 14 N-14-1 anchors, census tolerance rule
ref['A_mu_min_by_window'] = anchors
# zero side: A_zeta(m) = sum_{gamma>0} 2 B(gamma) cos(m delta gamma), B(r) = 2(1 - cos(delta r))/(delta r^2)
def B_tent(r, delta): return 2*(1 - np.cos(delta*r))/(delta*r*r)
Bz = B_tent(Z, DELTA)
Azeta = np.array([float(np.sum(2*Bz*np.cos(m*DELTA*Z))) for m in range(MBIG+1)])
ref['A_zeta_first_12'] = [float(v) for v in Azeta[:12]]
ref['A_vs_Azeta_max_abs_301'] = float(np.max(np.abs(A - Azeta))); ref['A_vs_Azeta_rms_301'] = float(np.sqrt(np.mean((A - Azeta)**2)))
ref['A_vs_Azeta_first120zeros_max'] = float(np.max(np.abs(A - np.array([float(np.sum(2*Bz[:120]*np.cos(m*DELTA*Z[:120]))) for m in range(MBIG+1)]))))
assert abs(ref['A_vs_Azeta_first120zeros_max'] - 0.23) < 0.03, ref['A_vs_Azeta_first120zeros_max']    # Paper 14: max |A - A_zeta| = 0.23 with 120 zeros
# ablations (E-WALL / E-HECKE-1 E1a language, zeta scaffold)
M40 = 60
ref['wall_pole_only'] = wall_scan(pole_array(DELTA, M40), DELTA, 3.0)[0]; assert abs(ref['wall_pole_only'] - 2*DELTA) < 1e-9
ref['wall_arch_zeta_only'] = wall_scan(arch_array_zeta(DELTA, M40), DELTA, 3.0)[0]
ref['wall_pole_plus_arch_zeta'] = wall_scan(pole_array(DELTA, M40) + arch_array_zeta(DELTA, M40), DELTA, 3.0)[0]; assert abs(ref['wall_pole_plus_arch_zeta'] - 0.80) < 1e-9   # E-WALL N=1 arm
s21 = sieve_bytes(21); PR21 = prime_array_zeta(DELTA, M40, 21, s21)
tbl = {}
for Mi in (14, 16, 30, 60):
    W = Mi*DELTA
    tbl[str(round(W, 2))] = dict(no_primes=mu_min(POLE[:M40+1] + ARCH[:M40+1], Mi), primes_le_21=mu_min(POLE[:M40+1] - PR21 + ARCH[:M40+1], Mi), all_primes=mu_min(A, Mi))
ref['A_prime_bond_table'] = tbl
assert tbl['0.7']['no_primes'] > 0 and tbl['0.8']['no_primes'] < 0 and all(v['primes_le_21'] > 0 for v in tbl.values())
# per-prime restoration: first negative window as prime powers enter one at a time
restore = []
run = np.zeros(M40+1)
for n, p, e in prime_powers(21, s21):
    c = math.log(p)/math.sqrt(n)
    for m, w in tent_hits(n, DELTA, M40): run[m] += c*w
    restore.append(dict(n=n, wall=wall_scan(POLE[:M40+1] - run + ARCH[:M40+1], DELTA, 3.0)[0]))
ref['A_restore_walls'] = restore
print(f"[A(m)] anchors {anchors}  walls pole-only {ref['wall_pole_only']} arch-only {ref['wall_arch_zeta_only']} scaffold {ref['wall_pole_plus_arch_zeta']}  ({time.time()-T0:.1f}s)")

# ---------- write embed ----------
embed = dict(zeros=[float(f"{z:.9f}") for z in Z], T=T, ref=ref,
             source="Odlyzko zeros1 (first 10,000 of 100,000), sha256 3436c916a7878261ac183fd7b9448c9a4736b8bbccf1356874a6ce1788541632")
out = os.path.join(HERE, 'rh_calculator_embed.json')
json.dump(embed, open(out, 'w'), indent=None, separators=(',', ':'))
print("ALL TWIN ASSERTS PASSED ->", out, f"({os.path.getsize(out)} bytes, {time.time()-T0:.1f}s)")
for k in ['li_lambda_n1_10', 'p2_plus_tail', 'hankel_lmin_live', 'hankel_corr_lmin', 'hankel_plant_zero1_beta0p6_lmin', 'hankel_kstar_law_traitor', 'dh_N_below',
          'A_first_12', 'A_zeta_first_12', 'A_lag300', 'A_vs_Azeta_max_abs_301', 'A_vs_Azeta_rms_301', 'A_vs_Azeta_first120zeros_max', 'A_prime_bond_table', 'A_restore_walls',
          'psi_explicit_x1000', 'psi_true_x1000', 'psi_explicit_x100000', 'psi_true_x100000', 'lagarias_sum_0p75_50', 'count_zeros_le_1000']:
    print(f"  {k}: {ref[k]}")
