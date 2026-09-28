/* rh_math.js — shared numeric core of Room XI "The Equals Sign" (the RH calculator), v0.2.
   Runs in the browser (inlined into index.html by build_page.py) and in Node (receipts/check_page_math.js).
   Every routine mirrors rh_calculator_twin.py; Paper 14's A(m) is ported from ehecke_census.py (E-HECKE-1). */
(function (root) {
'use strict';
const TWO_PI = 2*Math.PI, GAMMA_E = 0.57721566490153286061;
const CATALAN = 0.915965594177219015054603514932;
const PSI_QUARTER = -GAMMA_E - Math.PI/2 - 3*Math.log(2);        // psi(1/4)
const PSI1_QUARTER = Math.PI*Math.PI + 8*CATALAN;                 // psi'(1/4)
const LOGPI = Math.log(Math.PI);
const EPS = 2.0**-52;
let ZEROS = [], T = 0;
function setZeros(z) { ZEROS = z; T = z[z.length-1]; }

// ---------------- complex ----------------
const C = (re, im=0) => ({re, im});
const cadd=(a,b)=>C(a.re+b.re,a.im+b.im), csub=(a,b)=>C(a.re-b.re,a.im-b.im);
const cmul=(a,b)=>C(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re);
const cdiv=(a,b)=>{const d=b.re*b.re+b.im*b.im; return C((a.re*b.re+a.im*b.im)/d,(a.im*b.re-a.re*b.im)/d);};
const cabs=a=>Math.hypot(a.re,a.im), carg=a=>Math.atan2(a.im,a.re), conj=a=>C(a.re,-a.im);
const cexp=a=>{const e=Math.exp(a.re); return C(e*Math.cos(a.im), e*Math.sin(a.im));};
const clog=a=>C(Math.log(cabs(a)), carg(a));
const cscale=(a,k)=>C(a.re*k,a.im*k);
const cpowr=(x,s)=>cexp(cscale(s,Math.log(x)));
const cinv=a=>cdiv(C(1),a);
const csin=z=>C(Math.sin(z.re)*Math.cosh(z.im), Math.cos(z.re)*Math.sinh(z.im));
function cpowint(a,n){let r=C(1),b=a;while(n>0){if(n&1)r=cmul(r,b);b=cmul(b,b);n=Math.floor(n/2);}return r;}

// ---------------- special functions ----------------
const LG=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
function clgamma(z){
  if(z.re<0.5) return csub(clog(cdiv(C(Math.PI), csin(cscale(z,Math.PI)))), clgamma(csub(C(1),z)));
  z=csub(z,C(1)); let x=C(LG[0]);
  for(let i=1;i<9;i++) x=cadd(x, cdiv(C(LG[i]), cadd(z,C(i))));
  const t=cadd(z,C(7.5));
  return cadd(cadd(C(0.5*Math.log(2*Math.PI)), cmul(cadd(z,C(0.5)), clog(t))), cadd(cscale(t,-1), clog(x)));
}
const B2K=[1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510,43867/798,-174611/330,854513/138,-236364091/2730];
function fact(n){let f=1;for(let i=2;i<=n;i++)f*=i;return f;}
function zetaEM(s){
  const N=Math.max(30, Math.trunc(Math.abs(s.im))+20);
  let tot=C(0);
  for(let n=1;n<N;n++) tot=cadd(tot, cpowr(n, cscale(s,-1)));
  tot=cadd(tot, cdiv(cpowr(N, csub(C(1),s)), csub(s,C(1))));
  tot=cadd(tot, cscale(cpowr(N, cscale(s,-1)), 0.5));
  let rising=s;
  for(let k=1;k<=12;k++){
    const e=csub(cscale(s,-1), C(2*k-1));
    tot=cadd(tot, cscale(cmul(rising, cpowr(N,e)), B2K[k-1]/fact(2*k)));
    rising=cmul(rising, cmul(cadd(s,C(2*k-1)), cadd(s,C(2*k))));
  }
  return tot;
}
function xi(s){
  const half=cscale(s,0.5);
  const pref=cexp(cadd(cscale(half,-LOGPI), clgamma(half)));
  return cmul(cmul(cscale(cmul(s,csub(s,C(1))),0.5), pref), zetaEM(s));
}
const rsTheta=t=>clgamma(C(0.25,t/2)).im - t/2*LOGPI;
function li(x){
  if(x<=1) return -Infinity;
  const lx=Math.log(x); let s=0,f=1,inner=0;
  for(let n=1;n<200;n++){ f*=n; if((n-1)%2===0) inner+=1/n;
    const term=((n-1)%2===0?1:-1)*Math.pow(lx,n)/(f*Math.pow(2,n-1))*inner; s+=term;
    if(Math.abs(term)<1e-17*Math.abs(s)&&n>20) break; }
  return GAMMA_E+Math.log(lx)+Math.sqrt(x)*s;
}

// ---------------- zero side ----------------
const tailInvG2=Tt=>(Math.log(Tt/TWO_PI)+1)/(TWO_PI*Tt);
const momentTail=(n,Tt)=>Math.pow(Tt,1-2*n)/(TWO_PI*(2*n-1))*(Math.log(Tt/TWO_PI)+1/(2*n-1));
const rimTheta=g=>2*Math.atan(1/(2*g));
function liLambda(n){let s=0;for(const g of ZEROS){const v=Math.sin(n*rimTheta(g)/2);s+=4*v*v;}return s+n*n*tailInvG2(T);}
const NMain=Tt=>(Tt/TWO_PI)*Math.log(Tt/(TWO_PI*Math.E))+7/8;
function efTerm(x,beta,g){const rho=C(beta,g);return -2*cdiv(cpowr(x,rho),rho).re;}
function psiExplicit(x){let s=x-Math.log(TWO_PI)-0.5*Math.log(1-Math.pow(x,-2));const lx=Math.log(x),sx=Math.sqrt(x);
  for(const g of ZEROS) s-=2*sx*(0.5*Math.cos(g*lx)+g*Math.sin(g*lx))/(0.25+g*g); return s;}
function lagariasSum(s){let tot=C(0);for(const g of ZEROS){tot=cadd(tot,cadd(cinv(csub(s,C(0.5,g))),cinv(csub(s,C(0.5,-g)))));}
  return cadd(tot, cscale(csub(s,C(0.5)), 2*tailInvG2(T)));}
function liQuadContrib(rho,n){const m=csub(C(1),conj(rho)); let s=0; for(const r of [rho,conj(rho),m,conj(m)]) s+=1-cpowint(csub(C(1),cinv(r)),n).re; return s;}
function firstNegativeLi(rho,max=6000){for(let n=1;n<=max;n++) if(liQuadContrib(rho,n)<0) return n; return null;}
function countZerosLE(x){let c=0;for(const g of ZEROS){if(g<=x)c++;else break;}return c;}

// ---------------- integer side ----------------
let SV=null; const SIEVE_CAP=2000000;
function buildSieve(){ if(SV) return SV; const N=SIEVE_CAP; const spf=new Int32Array(N+1);
  for(let i=2;i<=N;i++){ if(spf[i]===0){ for(let j=i;j<=N;j+=i) if(spf[j]===0) spf[j]=i; } }
  const mu=new Int8Array(N+1), om=new Int16Array(N+1), pp=new Uint8Array(N+1); mu[1]=1;
  for(let n=2;n<=N;n++){ const p=spf[n], m=n/p; om[n]=om[m]+1; mu[n]=(spf[m]===p)?0:-mu[m]; pp[n]=(m===1||(pp[m]&&spf[m]===p))?1:0; }
  SV={spf,mu,om,pp,N}; return SV; }
function cumulativeAt(x){ const S=buildSieve(); x=Math.min(Math.floor(x),S.N); let pi=0,theta=0,psi=0,M=1,L=1,race=0;
  for(let n=2;n<=x;n++){ const p=S.spf[n]; if(p===n){pi++;theta+=Math.log(n); if(n%4===1)race-=1; else if(n%4===3)race+=1;}
    if(S.pp[n]) psi+=Math.log(p); M+=S.mu[n]; L+=(S.om[n]&1)?-1:1; }
  return {x,pi,theta,psi,M,L,race}; }
function factor(n){const f=[];let d=2;while(d*d<=n){while(n%d===0){const last=f[f.length-1];if(last&&last[0]===d)last[1]++;else f.push([d,1]);n/=d;}d+=(d===2)?1:2;}if(n>1)f.push([n,1]);return f;}
const isPrime=n=>{ if(n<2) return false; const f=factor(n); return f.length===1&&f[0][1]===1; };
const sigmaOf=f=>f.reduce((s,[p,e])=>s*((Math.pow(p,e+1)-1)/(p-1)),1);
const sopfrOf=f=>f.reduce((s,[p,e])=>s+p*e,0);
const OmegaOf=f=>f.reduce((s,[,e])=>s+e,0);
const muOf=f=>f.some(([,e])=>e>1)?0:(f.length%2?-1:1);
function gaussian(p){ if(p===2) return {a:1,b:1,deg:45,cos4:-1,kind:'ramified: 2 = −i(1+i)²'}; if(p%4===3) return {kind:'inert: stays prime in ℤ[i], norm p²'};
  for(let b=1;b*b<=p;b++){const a2=p-b*b,a=Math.round(Math.sqrt(a2)); if(a*a===a2&&a>=b){const z=C(a,b);return {a,b,deg:Math.atan2(b,a)*180/Math.PI,cos4:cpowint(z,4).re/(p*p),kind:'split: p = a² + b² = (a+bi)(a−bi)'};}} return null; }
function eisenstein(p){ if(p===3) return {a:2,b:1,deg:30,cos6:-1,kind:'ramified: 3 = −ω²(1−ω)²'}; if(p%3===2) return {kind:'inert: stays prime in ℤ[ω], norm p²'};
  const w=C(-0.5,Math.sqrt(3)/2); const lim=Math.floor(Math.sqrt(p))+2;
  for(let b=1;b<=lim;b++) for(let a=b;a<=lim;a++) if(a*a-a*b+b*b===p){const z=cadd(C(a),cscale(w,b)); let th=((carg(z)*180/Math.PI%60)+60)%60; th=Math.min(th,60-th); return {a,b,deg:th,cos6:cpowint(z,6).re/Math.pow(p,3),kind:'split: p = a² − ab + b² = N(a + bω)'};} return null; }
function harmonic(n){let s=0;for(let k=1;k<=n;k++)s+=1/k;return s;}
function lagariasIneq(n,f){const H=harmonic(n),rhs=H+Math.exp(H)*Math.log(H),sg=sigmaOf(f);return {sigma:sg,H,rhs,margin:rhs-sg};}
const robinMargin=(n,f)=>Math.exp(GAMMA_E)*n*Math.log(Math.log(n))-sigmaOf(f);
function modpow(b,e,m){let r=1;b%=m;while(e>0){if(e&1)r=(r*b)%m;b=(b*b)%m;e=Math.floor(e/2);}return r;}
function hasse(p,a,b){ if((4*a*a*a+27*b*b)%p===0) return {singular:true}; let s=0;
  for(let x=0;x<p;x++){const r=((((x*x)%p*x)%p+a*x+b)%p+p)%p; if(r===0)continue; s+= (modpow(r,(p-1)/2,p)===1)?1:-1;}
  const ap=-s, pts=p+1-ap, disc=4*p-ap*ap; const alpha=C(ap/2, Math.sqrt(Math.max(disc,0))/2);
  return {ap,pts,alpha,absAlphaOverSqrtP:cabs(alpha)/Math.sqrt(p),ok:Math.abs(ap)<=2*Math.sqrt(p),disc}; }

// ---------------- symmetric eigenvalues: cyclic two-sided Jacobi (relative accuracy ~ kappa(D^-1/2 A D^-1/2) * eps, Demmel–Veselić) ----------------
function symEigs(Ain, opts={}){
  const n=Ain.length; const A=Ain.map(r=>Array.from(r)); const maxSweeps=opts.maxSweeps||80; const relTol=opts.relTol||1e-40;
  for(let sweep=0;sweep<maxSweeps;sweep++){
    let off=0, dg=0; for(let p=0;p<n;p++){ dg+=A[p][p]*A[p][p]; for(let q=p+1;q<n;q++) off+=A[p][q]*A[p][q]; }
    if(off===0 || off<=relTol*dg) break;
    for(let p=0;p<n-1;p++) for(let q=p+1;q<n;q++){
      const apq=A[p][q]; if(apq===0) continue;
      const app=A[p][p], aqq=A[q][q]; const theta=(aqq-app)/(2*apq);
      const t=(theta>=0?1:-1)/(Math.abs(theta)+Math.sqrt(theta*theta+1));
      const c=1/Math.sqrt(t*t+1), s=t*c;
      for(let r=0;r<n;r++){ if(r===p||r===q) continue; const arp=A[r][p], arq=A[r][q]; const nrp=c*arp-s*arq, nrq=s*arp+c*arq; A[r][p]=nrp; A[p][r]=nrp; A[r][q]=nrq; A[q][r]=nrq; }
      A[p][p]=app-t*apq; A[q][q]=aqq+t*apq; A[p][q]=0; A[q][p]=0;
    }
  }
  return A.map((r,i)=>r[i]).sort((a,b)=>a-b);
}
const lamMin=(A,opts)=>symEigs(A,opts)[0];
function scaledCondition(A){ // kappa of D^-1/2 A D^-1/2 (via its extreme eigenvalues)
  const n=A.length; const d=A.map((r,i)=>Math.sqrt(Math.abs(r[i]))); const Cm=A.map((r,i)=>r.map((v,j)=>v/(d[i]*d[j])));
  const ev=symEigs(Cm); return {lamMinC:ev[0], kappa:ev[n-1]/Math.abs(ev[0])}; }

// ---------------- live Hankel (TSF Theorem 4.1) ----------------
function hankelQ(M, plant){ // q_m = sum gamma^{-2(m+1)} + tail; plant = {g0, delta, replace} : the pair z = g0 -/+ i delta enters; if replace (default) the on-line zero g0 leaves
  const q=[]; for(let m=0;m<M;m++){ const n=m+1; let v=0; for(const g of ZEROS) v+=Math.pow(g,-2*n); v+=momentTail(n,T);
    if(plant){ if(plant.replace!==false) v-=Math.pow(plant.g0,-2*n); v+=2*cpowint(cinv(C(plant.g0,-plant.delta)),2*n).re; } q.push(v); } return q; }
function hankelMatrix(q,k){ return Array.from({length:k},(_,i)=>Array.from({length:k},(_,j)=>q[i+j])); }
function hankelReport(q, K){ const rows=[]; for(let k=1;k<=K;k++){ const H=hankelMatrix(q,k); const lm=lamMin(H); const sc=scaledCondition(H); const floor=sc.kappa*EPS;
    rows.push({k, lamMin:lm, corrLamMin:sc.lamMinC, relFloor:floor, certified: floor<1e-2 && lm>0, sign: lm>0?'+':'−'}); } return rows; }

// ---------------- Paper 14: A(m) = POLE(m) − PRIME(m) + ARCH(m) (ehecke_census.py port) ----------------
let PSIEVE=null;
function byteSieve(n){ if(PSIEVE&&PSIEVE.n>=n) return PSIEVE.s; const s=new Uint8Array(n+1).fill(1); s[0]=0;s[1]=0; for(let i=2;i*i<=n;i++){ if(s[i]){ for(let j=i*i;j<=n;j+=i) s[j]=0; } } PSIEVE={n,s}; return s; }
function tentHits(n, delta, M){ const ln=Math.log(n), m0=Math.trunc(ln/delta); const out=[]; for(const m of [m0-1,m0,m0+1,m0+2]){ if(m>=0&&m<=M){ const w=1-Math.abs(ln-m*delta)/delta; if(w>0) out.push([m,w]); } } return out; }
function primeArrayZeta(delta, M, cutoff){ const s=byteSieve(cutoff); const arr=new Float64Array(M+1);
  for(let p=2;p<=cutoff;p++){ if(!s[p]) continue; const lp=Math.log(p); let n=p; while(n<=cutoff){ const c=lp/Math.sqrt(n); for(const [m,w] of tentHits(n,delta,M)) arr[m]+=c*w; n*=p; } } return arr; }
function S_series(a0, step, t){ if(t===0) return PSI1_QUARTER/(step*step); let s=0, j=0; for(;;){ const a=a0+step*j; const term=Math.exp(-a*t)/(a*a); s+=term; if(term<1e-22 && a*t>1) break; j++; if(j>500000) break; } return s; }
function archArrayZeta(delta, M){ const a0=0.5, step=2.0, cst=PSI_QUARTER-LOGPI; const S=[]; for(let m=0;m<M+2;m++) S.push(S_series(a0,step,m*delta));
  const arr=new Float64Array(M+1); arr[0]=cst+(2.0/delta)*(S[0]-S[1]); for(let m=1;m<=M;m++) arr[m]=-(1.0/delta)*(S[m-1]-2*S[m]+S[m+1]); return arr; }
function poleArray(delta, M){ const Cd=(8.0/delta)*(Math.cosh(delta/2)-1.0); const arr=new Float64Array(M+1); for(let m=0;m<=M;m++) arr[m]=2*Cd*Math.cosh(m*delta/2); return arr; }
function toeplitz(arr, Mi){ return Array.from({length:Mi},(_,i)=>Array.from({length:Mi},(_,j)=>arr[Math.abs(i-j)])); }
function muMin(arr, Mi){ return lamMin(toeplitz(arr,Mi), {relTol:1e-34}); }
function wallScan(arr, delta, Wcap, undec=1e-8){ const Mcap=Math.round(Wcap/delta); for(let Mi=1;Mi<=Mcap;Mi++){ const mu=muMin(arr,Mi); if(mu<-undec) return {wall:Mi*delta, mu}; } return {wall:null, mu:null}; }
function aZetaSide(delta, M){ const arr=new Float64Array(M+1); for(let m=0;m<=M;m++){ let s=0; for(const g of ZEROS){ const B=2*(1-Math.cos(delta*g))/(delta*g*g); s+=2*B*Math.cos(m*delta*g); } arr[m]=s; } return arr; }
const ACACHE={};
function buildA(delta, M, cutoff){ const key=`${delta}|${M}|${cutoff}`; if(ACACHE[key]) return ACACHE[key]; const POLE=poleArray(delta,M), ARCH=archArrayZeta(delta,M), PRIME=primeArrayZeta(delta,M,cutoff); const A=new Float64Array(M+1); for(let m=0;m<=M;m++) A[m]=POLE[m]-PRIME[m]+ARCH[m]; ACACHE[key]={POLE,ARCH,PRIME,A}; return ACACHE[key]; }

// ---------------- self-check against the twin ----------------
function selfCheck(REF){
  const items=[]; const ok=(name,got,want,tol)=>{ let pass; if(typeof want==='number'&&typeof got==='number') pass=Math.abs(got-want)<=tol*Math.max(1,Math.abs(want)); else if(Array.isArray(want)&&Array.isArray(got)&&typeof want[0]==='number'&&typeof tol==='number') pass=got.length===want.length&&got.every((v,i)=>Math.abs(v-want[i])<=tol*Math.max(1,Math.abs(want[i]))); else pass=JSON.stringify(got)===JSON.stringify(want); items.push({name,got,want,pass}); };
  const R=REF, deg=r=>r*180/Math.PI;
  ok('rim angle θ₁ (deg)', deg(rimTheta(ZEROS[0])), R.rim_theta_deg_zero1, 1e-12);
  ok('λ₁…λ₁₀ from 10,000 zeros + tail', [1,2,3,4,5,6,7,8,9,10].map(liLambda), R.li_lambda_n1_10, 1e-10);
  ok('p₂ over 10,000 zeros + tail', ZEROS.reduce((s,g)=>s+1/(g*g),0)+momentTail(1,T), R.p2_plus_tail, 1e-12);
  ok('count of zeros ≤ 100', countZerosLE(100), R.count_zeros_le_100, 0);
  ok('count of zeros ≤ 1000', countZerosLE(1000), R.count_zeros_le_1000, 0);
  ok('N(100) main term', NMain(100), R.N_main_100, 1e-12);
  ok('unfolded gap γ₁→γ₂', (ZEROS[1]-ZEROS[0])*Math.log(ZEROS[0]/TWO_PI)/TWO_PI, R.unfolded_gap_1_2, 1e-12);
  ok('explicit-formula wave of ρ₁ at x = 1000', efTerm(1000,0.5,ZEROS[0]), R.ef_term_zero1_x1000, 1e-10);
  const rho=C(R.dh_beta,R.dh_gamma), r2=R.dh_beta**2+R.dh_gamma**2;
  ok('traitor: reciprocal-circle defect', (1-2*R.dh_beta)/r2, R.dh_recip_circle_defect, 1e-12);
  ok('traitor: pair forced sign', (2*R.dh_beta-1)*(1/cabs(csub(C(1),conj(rho)))**2-1/r2), R.dh_pair_forced_sign, 1e-9);
  ok('traitor: first negative Li share n', firstNegativeLi(rho), R.dh_li_first_negative_n, 0);
  ok('traitor: Li share at n = 80', liQuadContrib(rho,80), R.dh_li_contrib_n80, 1e-9);
  ok('traitor: Disc radius |1 − 1/ρ|', cabs(csub(C(1),cinv(rho))), R.dh_disc_radius, 1e-12);
  ok('traitor: Weil leg angle (a = 0.1)', deg(0.01*R.dh_gamma*(2*R.dh_beta-1)), R.dh_weil_leg_deg_a0p1, 1e-12);
  ok('traitor: visibility a_min', Math.sqrt(Math.PI/(4*(R.dh_beta-0.5)*R.dh_gamma)), R.dh_visibility_a_min, 1e-12);
  ok('traitor: true zeros below its height', countZerosLE(R.dh_gamma), R.dh_N_below, 0);
  const g13=gaussian(13), e7=eisenstein(7), e13=eisenstein(13);
  ok('13 in ℤ[i]: a, b', [g13.a,g13.b], [R.gauss_13[0],R.gauss_13[1]], 0); ok('13 in ℤ[i]: angle', g13.deg, R.gauss_13[2], 1e-9); ok('13 in ℤ[i]: cos4θ = −119/169', g13.cos4, R.gauss_13[3], 1e-12);
  ok('7 in ℤ[ω]: cos6θ = −143/343', e7.cos6, R.eis_7[3], 1e-9); ok('7 in ℤ[ω]: folded angle', e7.deg, R.eis_7[2], 1e-9); ok('13 in ℤ[ω]: cos6θ = 253/2197', e13.cos6, R.eis_13[3], 1e-9);
  const fp=[]; for(let n=2;n<=100;n++){ if(sopfrOf(factor(n))===n) fp.push(n); } ok('sopfr fixed points ≤ 100', fp, R.sopfr_fixed_points_le_100, undefined);
  ok('Lagarias margin at 5040', lagariasIneq(5040,factor(5040)).margin, R.lagarias_5040.margin, 1e-9);
  ok('Robin margin at 10080', robinMargin(10080,factor(10080)), R.robin_margin_10080, 1e-9);
  for(const x of [10000,100000,1000000]){ const c=cumulativeAt(x), t=R.sieve_table[String(x)];
    ok(`π(${x})`, c.pi, t.pi, 0); ok(`M(${x})`, c.M, t.M, 0); ok(`L(${x})`, c.L, t.L, 0); ok(`race(${x})`, c.race, t.race, 0);
    ok(`ψ(${x})`, c.psi, t.psi, 1e-11); ok(`θ(${x})`, c.theta, t.theta, 1e-11); ok(`li(${x})`, li(x), t.li, 1e-11); }
  ok('ψ₁₀ₖ(1000.5)', psiExplicit(1000.5), R.psi_explicit_x1000, 1e-10);
  ok('ψ₁₀ₖ(10000.5)', psiExplicit(10000.5), R.psi_explicit_x10000, 1e-10);
  ok('ψ₁₀ₖ(100000.5)', psiExplicit(100000.5), R.psi_explicit_x100000, 1e-10);
  ok('Hasse a_p, y² = x³ + x + 1 over F₁₀₀₉', hasse(1009,1,1).ap, R.hasse_1009.a_p, 0);
  ok('Hasse a_p over F₁₀₀₀₇', hasse(10007,1,1).ap, R.hasse_10007.a_p, 0);
  for(const [k,v] of Object.entries(R.zeta_checks)){ const z=zetaEM(C(v.s[0],v.s[1])); const pass=Math.hypot(z.re-v.zeta_em[0], z.im-v.zeta_em[1])<=1e-11*Math.max(1,Math.hypot(v.zeta_em[0],v.zeta_em[1])); items.push({name:`ζ at ${k}`, got:[z.re,z.im], want:v.zeta_em, pass}); }
  ok('ζ(2) = π²/6', zetaEM(C(2)).re, Math.PI**2/6, 1e-13);
  ok('2ξ(½)', 2*xi(C(0.5)).re, 2*R.xi_half, 1e-11);
  ok('|ξ(ρ₁)| ≈ 0', cabs(xi(C(0.5,ZEROS[0]))), 0, 1e-10);
  ok('Riemann–Siegel θ(γ₁)', rsTheta(ZEROS[0]), R.rs_theta_at_gamma1, 1e-10);
  ok('Riemann–Siegel θ(100)', rsTheta(100), R.rs_theta_100, 1e-10);
  ok('Z(100.5)', cmul(cexp(C(0,rsTheta(100.5))),zetaEM(C(0.5,100.5))).re, R.Z_100p5, 1e-10);
  const L=lagariasSum(C(0.75,50)); items.push({name:'ξ′/ξ(0.75+50i), 10,000 zeros + tail', got:[L.re,L.im], want:R.lagarias_sum_0p75_50, pass: Math.hypot(L.re-R.lagarias_sum_0p75_50[0], L.im-R.lagarias_sum_0p75_50[1])<=1e-10});
  ok('its angle (deg)', deg(carg(L)), R.lagarias_angle_deg, 1e-9);
  // live Hankel
  const q=hankelQ(19); ok('Hankel moments q₀…q₁₈ (10,000 zeros + tail)', q, R.hankel_q, 1e-10);
  const hr=hankelReport(q,7);
  ok('λ_min(H_k) live, k = 1…5 (Jacobi vs LAPACK)', hr.slice(0,5).map(r=>r.lamMin), R.hankel_lmin_live.slice(0,5), 1e-6);
  ok('λ_min(H₆) live (κε floor ~1e-8)', hr[5].lamMin, R.hankel_lmin_live[5], 1e-4);
  ok('λ_min(H₇) live (κε floor ~1e-5)', hr[6].lamMin, R.hankel_lmin_live[6], 1e-2);
  ok('correlation-normalized λ_min, k = 1…7', hr.map(r=>r.corrLamMin), R.hankel_corr_lmin, 1e-6);
  const qp=hankelQ(19,{g0:ZEROS[0],delta:0.1}); const hp=hankelReport(qp,6); ok('displaced ρ₁ to β = 0.6: first negative k', hp.findIndex(r=>r.lamMin<0)+1, R.hankel_plant_zero1_beta0p6_kstar, 0);
  ok('displaced ρ₁ to β = 0.6: λ_min(H₅)', hp[4].lamMin, R.hankel_plant_zero1_beta0p6_lmin[4], 1e-4);
  // Paper 14 A(m)
  const delta=R.A_delta, MB=300; const {POLE,ARCH,PRIME,A}=buildA(delta,MB,R.A_cutoff);
  ok('A(m), m = 0…11 (δ = 0.05, cutoff 3,445,000)', Array.from(A.slice(0,12)), R.A_first_12, 1e-10);
  ok('POLE(m), m = 0…5', Array.from(POLE.slice(0,6)), R.POLE_first_6, 1e-12); ok('ARCH(m), m = 0…5', Array.from(ARCH.slice(0,6)), R.ARCH_first_6, 1e-10); ok('PRIME(m), m = 0…5', Array.from(PRIME.slice(0,6)), R.PRIME_first_6, 1e-10);
  ok('lag 300: POLE, PRIME, A', [POLE[300],PRIME[300],A[300]], [R.A_lag300.POLE,R.A_lag300.PRIME,R.A_lag300.A], 1e-9);
  const Az=aZetaSide(delta,MB); ok('zero side Aζ(m), m = 0…11', Array.from(Az.slice(0,12)), R.A_zeta_first_12, 1e-9);
  let mx=0, ss=0; for(let m=0;m<=MB;m++){ const d=Math.abs(A[m]-Az[m]); mx=Math.max(mx,d); ss+=d*d; } ok('max |A − Aζ| over 301 lags', mx, R.A_vs_Azeta_max_abs_301, 1e-8); ok('rms |A − Aζ|', Math.sqrt(ss/(MB+1)), R.A_vs_Azeta_rms_301, 1e-8);
  ok('μ_min(T_M) at window 2.5', muMin(A,50), R.A_mu_min_by_window['2.5'], 1e-8); ok('μ_min at window 5', muMin(A,100), R.A_mu_min_by_window['5.0'], 1e-8); ok('μ_min at window 7.5', muMin(A,150), R.A_mu_min_by_window['7.5'], 1e-7);
  const P60=poleArray(delta,60), R60=archArrayZeta(delta,60);
  ok('wall, POLE only (= 2δ)', wallScan(P60,delta,3).wall, R.wall_pole_only, 1e-12); ok('wall, ARCH_ζ only', wallScan(R60,delta,3).wall, R.wall_arch_zeta_only, 1e-12);
  const scaf=new Float64Array(61); for(let m=0;m<=60;m++) scaf[m]=P60[m]+R60[m]; ok('wall, POLE + ARCH_ζ scaffold (E-WALL 0.80)', wallScan(scaf,delta,3).wall, R.wall_pole_plus_arch_zeta, 1e-12);
  const PR21=primeArrayZeta(delta,60,21); const bond=[]; const bondKeys=Object.keys(R.A_prime_bond_table);
  for(const key of bondKeys){ const Mi=Math.round(parseFloat(key)/delta); const np_=new Float64Array(61), wp=new Float64Array(61); for(let m=0;m<=60;m++){ np_[m]=P60[m]+R60[m]; wp[m]=P60[m]-PR21[m]+R60[m]; } bond.push({key, Mi, noPrimes:muMin(np_,Mi), primesLe21:muMin(wp,Mi), all:muMin(A,Mi)}); }
  ok('prime–pole bond table (no primes / ≤ 21 / all) at windows 0.7, 0.8, 1.5, 3', bond.map(b=>[b.noPrimes,b.primesLe21,b.all]).flat(), bondKeys.map(k=>R.A_prime_bond_table[k]).map(v=>[v.no_primes,v.primes_le_21,v.all_primes]).flat(), 1e-8);
  const s21=byteSieve(21); const run=new Float64Array(61); const walls=[]; const pps=[]; for(let p=2;p<=21;p++){ if(!s21[p]) continue; let n=p; while(n<=21){ pps.push([n,p]); n*=p; } } pps.sort((a,b)=>a[0]-b[0]);
  for(const [n,p] of pps){ const c=Math.log(p)/Math.sqrt(n); for(const [m,w] of tentHits(n,delta,60)) run[m]+=c*w; const arr=new Float64Array(61); for(let m=0;m<=60;m++) arr[m]=P60[m]-run[m]+R60[m]; walls.push(wallScan(arr,delta,3).wall); }
  ok('per-prime-power restoration walls', walls, R.A_restore_walls.map(r=>r.wall), undefined);
  ok('ψ(¼), ψ′(¼), ln π', [PSI_QUARTER,PSI1_QUARTER,LOGPI], [R.psi_quarter,R.psi1_quarter,R.log_pi], 1e-14);
  const pass=items.filter(i=>i.pass).length;
  return {items, pass, total:items.length, bond, walls, A, Az, POLE, ARCH, PRIME};
}

root.RH = { setZeros, get ZEROS(){return ZEROS;}, get T(){return T;}, TWO_PI, GAMMA_E, EPS, PSI_QUARTER, PSI1_QUARTER, LOGPI,
  C, cadd, csub, cmul, cdiv, cabs, carg, conj, cexp, clog, cscale, cpowr, cinv, cpowint,
  clgamma, zetaEM, xi, rsTheta, li, tailInvG2, momentTail, rimTheta, liLambda, NMain, efTerm, psiExplicit, lagariasSum, liQuadContrib, firstNegativeLi, countZerosLE,
  buildSieve, cumulativeAt, factor, isPrime, sigmaOf, sopfrOf, OmegaOf, muOf, gaussian, eisenstein, harmonic, lagariasIneq, robinMargin, hasse, SIEVE_CAP,
  symEigs, lamMin, scaledCondition, hankelQ, hankelMatrix, hankelReport,
  byteSieve, tentHits, primeArrayZeta, archArrayZeta, poleArray, toeplitz, muMin, wallScan, aZetaSide, buildA, selfCheck };
})(typeof window !== 'undefined' ? window : globalThis);
