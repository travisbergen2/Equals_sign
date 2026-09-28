# Room XI · The Equals Sign — an RH calculator (v0.2)

Type a zero (`#7` or an ordinate), an integer, a height `x`, a point `σ+it`, or a curve `p=1009 a=1 b=1`
over F_p. Press `=`.

- **The dial** — the Bergen Disc, z = 1 − 1/ρ: pole at the center, critical line on the rim, s = ∞ at the
  vanishing point; all 10,000 tabulated zeros as rim ticks (they crowd within θ₁ = 4.05° of z = 1), the current
  object and its mirror orbit plotted (radial offset exaggerated by a stated factor when off the rim), and an
  unrolled rim on a log scale showing the 1/γ crowding.
- **The tower** — the live Hankel test of Two Spectral Fingerprints, Theorem 4.1: q_m = Σγ^{−2(m+1)} over the
  10,000 zeros plus the Riemann–von Mangoldt tail, H_k = [q_{i+j}], λ_min(H_k) for k = 1…10 by cyclic Jacobi,
  next to the exact Newton-route values (k ≤ 7) and the certification floor κ(D^{−½}H_kD^{−½})·ε
  (Demmel–Veselić). Displace a zero off the line (β field, or the preset) and the tower reports the first
  certified negative eigenvalue k*; the traitor zero shows why its Hankel signature sits beyond double precision
  while the Li row detects it at n = 537.
- **The windows** — Paper 14's A(m) = POLE(m) − PRIME(m) + ARCH(m) at δ = 0.05 with prime powers to 3,445,000,
  ported verbatim from the E-HECKE-1 census implementation: μ_min(T_M) at windows 2.5 … 15 against the paper's
  N-14-1 anchors, the explicit formula live (arithmetic A(m) against the zero side Σ2B(γ)cos(mδγ) over 10,000
  zeros), the lag-300 cancellation, the E-WALL ablations (POLE-only wall 2δ, ARCH-only, POLE+ARCH scaffold 0.80),
  the prime–pole bond table and the one-prime-at-a-time restoration walls, and the P-14-1 conditional certificate.
- **Dictionary / board** — every representation of the object, and what each RH-equivalent says about it, tagged
  `FREE` (identity; detects nothing), `EQUIVALENT` (theorem "RH ⟺ …", can fail), `MEASURED`, `DISPROVEN`
  (Mertens, Pólya), `PROVEN WORLD` (curves over F_p — Hasse).

## Files

- `index.html` — the self-contained page (one external request: the Google Fonts stylesheet; no analytics);
  embeds the first 10,000 ordinates and the twin's reference values; self-checks 89 items on load.
- `receipts/rh_math.js` — the numeric core shared by the page and Node.
- `receipts/check_page_math.js` — `node check_page_math.js` runs the page's own self-check headlessly.
- `receipts/rh_calculator_twin.py` — the Python twin: reference values, mpmath verification of the
  double-precision algorithms, the ported A(m) construction with its census asserts. Needs `zeros1.txt`
  (Odlyzko's first 100,000 zeros, https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1 ,
  sha256 `3436c916a7878261ac183fd7b9448c9a4736b8bbccf1356874a6ce1788541632`).
- `receipts/build_page.py` — inlines the embed JSON and `rh_math.js` into `receipts/index.template.html`.
- `receipts/rh_calculator_embed.json` — first 10,000 ordinates (9 dp) + reference values.

Run: `python3 receipts/rh_calculator_twin.py && node receipts/check_page_math.js && python3 receipts/build_page.py`.

## Algorithms (browser and twin share them)

ζ(s): Euler–Maclaurin, N = ⌊|t|⌋ + 20 direct terms, 12 Bernoulli terms (≤ 2e-13 vs mpmath on the test points).
log Γ: Lanczos (g = 7) with reflection. li(x): Ramanujan's series. Sieve: smallest-prime-factor table to 2·10⁶
(counting functions) and a byte sieve to 3,445,000 (Paper 14's prime term). Eigenvalues: cyclic two-sided Jacobi.
Li: λₙ = Σ_pairs 4 sin²(nθ/2) plus the tail n²(ln(T/2π)+1)/(2πT). Hankel moments: Σγ^{−2n} plus the tail
T^{1−2n}(ln(T/2π) + 1/(2n−1))/(2π(2n−1)). A(m): POLE(m) = 2C_δ cosh(mδ/2); PRIME(m) = Σ Λ(n)n^{−½} tent_δ(ln n − mδ);
ARCH(m) by the exact-series reduction with S(t) = Σ_j e^{−(2j+½)t}/(2j+½)², S(0) = ψ′(¼)/4, constant ψ(¼) − ln π.

## Fence

Nothing here proves or disproves RH. Every row is a finite instance of an equivalent statement at its grade.
Citations flagged from memory in the twin's comments: Lagarias 1999, Robin 1984, Hasse 1933, Schoenfeld 1976,
Bays–Hudson 2000, Rubinstein–Sarnak 1994, Barnet-Lamb–Geraghty–Harris–Taylor 2011, Demmel–Veselić 1992.
