\\                              PORTHOS

\\ A package of Pari/GP routines for working with orthogonal polynomials
\\   Dirk Laurie   dlaurie@na-net.ornl.gov   

\\ Quick help: inside Pari/GP, type 'help'.  Everything useful in here
\\   can be found there too. If all else fails: read the TeX documentation.
\\   TODO: write the TeX documentation.

\\ !! Development notes:
\\   1. If a routine has an "addhelp", that routine has been briefly tested 
\\      during an interactive session.
\\   2. The TeX documentation does not exist yet.

\\ Global variables
{if(type(PORTHOS_VERSION)!="t_INT",
  global(WARN_SHIFT,WARN_NONPOSITIVE,WARN_LOSS,NORMALIZED,FOOLHARDY,
    PORTHOS_VERSION, KNOWN_POLYNOMIALS, SHORT_POLYNOMIALS, DEBUG, LASTCS))}
WARN_SHIFT=NORMALIZED=FOOLHARDY=DEBUG=0; WARN_NONPOSITIVE=WARN_LOSS=1;
{KNOWN_POLYNOMIALS=["Jacobi","Gegenbauer","Legendre","Chebyshev","Laguerre",
  "Hermite","Refinable"]}
{SHORT_POLYNOMIALS=["Jac","Geg","Leg","Che","Lag","Her","Ref"]}
\\ WARN_LOSS has not been implemented yet.  It should issue a warning
\\ whenever something is done that discards information, e.g. calculating
\\ nodes and weights from recursion coefficients defined for even degree.
PORTHOS_VERSION=60;
\\ 60: 1 May 2016, incorporates GP language changes: 
\\    '||' for '|', 'my' for 'local'
\\ 54: improvements to help only
\\ 53: 4 December 2008

gp_version = subst(Pol(version),x,10)
incompatible = 0
if (gp_version<240, if(PORTHOS_VERSION>=60, incompatible="old"))
if (gp_version>=240, if(PORTHOS_VERSION<60, incompatible="new"))
{ if (incompatible, print("Your version of Pari-GP is too "incompatible
" for your Porthos.\n Both [Pari>=2.4,Porthos>=0.60] or neither must hold.")) } 

\\ Routines that create or modify two-term coefficients

{qe_radau(qe,x0,n)=my(x,P); if(n==0, n=length(qe[2])+1);
  if(length(qe[1])<n, qe[1]=concat(qe[1],0)); 
  P=stqd(qe,x0); qe[1][n]=qe[1][n]-P[1][n]; qe}
{addhelp(qe_radau,
   "qe_radau(qe,{x0=0},{n=0}): modify qe to give n-point Radau formula with node at x0.  n=0 means largest possible n")}

\\ !! Present version, using stqd, might fail needlessly.  Rewrite this to use
\\ pu_qe instead.  TODO I no longer understand this comment. DPL 2017-05
{qe_lobatto(qe,x0,n)=my(x,P); if(n==0, n=length(qe[1])+1);
  if(length(qe[1])<n, qe[1]=concat(qe[1],0));
  if(length(qe[2])<n-1, qe[2]=concat(qe[2],0)); 
    if(x0==0, qe[2][n-1]=qe[2][n-1]+qe[1][n],
      P=stqd(qe,x0); qe[2][n-1]=x0/(1-qe[1][n-1]/P[1][n-1]));
      qe[1][n]=0; qe}
{addhelp(qe_lobatto,
 "qe_lobatto(qe,x0,{n=0}): modify qe to give n-point Lobatto formula with nodes at 0 and x0.  n=0 means largest possible n")}

{qe_antigauss(qe)=my(n); n=length(qe[1]);
   qe[1][n]=qe[1][n]-qe[2][n-1]; qe[2][n-1]=2*qe[2][n-1]; qe}  
{addhelp(qe_antigauss,
   "qe_antigauss(qe): modify qe to give anti-Gauss formula")} 

{qe_shift(qe,int1,int2)=my(f,h); if(int2==0, int2=[0,2]);
   f=(int1[2]-int1[1])/(int2[2]-int2[1]); h=int1[1]-f*int2[1];
   stqd(qe*f,-h)}
addhelp(qe_shift, "qe_shift(qe,int1,{int2=[0,2]}): shift two-term coefficients from interval int2 to interval int1")

{qe_jacobi(m,alf,bet)=my(kap,n); n=floor((m+1)/2);   
  [vector(n,k, kap=2*k+alf+bet; 2*(k+bet)/kap*if(k>1,(kap-k)/(kap-1),1)),
   vector(m-n,k, kap=2*k+alf+bet; 2*k*(k+alf)/(kap*(kap+1))),
   if(NORMALIZED,1,2^(alf+bet+1)*fac(alf)*fac(bet)/fac(alf+bet+1))]}   
{addhelp(qe_jacobi,
   "qe_jacobi(m,{alf=0},{bet=0}): two-term Jacobi coefficients for degree m over [0,2]")}

{qe_laguerre(m,alf)=my(n); n=floor((m+1)/2);
  [vector(n,k, k+alf), vector(m-n,k,k), if(NORMALIZED,1,fac(alf))]}
{addhelp(qe_laguerre,
   "qe_jacobi(m,{alf=0}): two-term Laguerre coefficients for degree m")}

{qe_anti(qe,Z,n)=my(mu);
  if(n==0, n=length(qe[1]); if(exactness(qe,Z)<2*n-1, n=n-1)); 
  if(length(Z[1])<2*n+1, error("qe_anti: too few coefficients in Z")); 
  mu=mu_qe(qe,Z); for(k=2,2*n+2,mu[k]=-mu[k]); qe=qe_mu(mu,Z,n+1);
  concat(qe,Z[3])}
{addhelp(qe_anti,"qe_anti(qe,Z,{n=0}): n-point anti-formula of Z "
 "w.r.t. weight qe. n=0 means smallest possible n")}

{vqe_stratify(qe,Z,k)=my(qe1,first,m,n,mu1,mu2);
  if(Z==0, Z=[[qe[1][1]],[],qe[3]]; k=2, 
    m=length(Z[1]); k=if(exactness(qe,Z)==2*m-1, m+1, m-1)); 
  qe1=[Z]; first=1; m=length(qe[1])+1; n=length(Z[1])+k;  
  while(n<m, mu2=mu_qe(Z,qe); 
    if(first,   mu1=mu2; first=0,   mu1=(mu1+mu2)/2; mu2=mu1);
    for(j=k+1,2*k, mu2[j]=-mu2[j]);
    Z=qe_mu(mu2,qe,k); qe1=concat(qe1,[Z]); k=2*k; n=n+k); qe1}
{addhelp(vqe_stratify,
  "vqe_stratify(qe,{Z=Z0},{k=k0}): recursion coefficients for stratified"
  " sequence of formulas for qe, starting at Z, adding k points the first"
  " time round.  Z0 is the one-point Gaussian formula of qe, k0 is the"
  " smallest possible k.  The result is a vector, each component of which"
  " is a qe.")} 

\\ Routines that create or modify three-term coefficients

{ab_radau(ab,xn)=my(n,p); n=length(ab[2]);
  if(length(ab[1])<n, ab[1]=concat(ab[1],0));
  p=p_ab(ab,xn); ab[1][n]=ab[1][n]+p[n+1]/p[n]; ab}
{addhelp(ab_radau,
   "ab_radau(ab,{b=0},{n=0}): modify qe to give n-point Radau formula with node at b.  n=0 means largest possible n")}

{ab_shift(ab,int1,int2)=my(f,h,b0); if(int2==0, int2=[-1,1]);
   f=(int1[2]-int1[1])/(int2[2]-int2[1]); h=int1[1]-f*int2[1];
   ab[1]=vector(length(ab[1]),j,f*ab[1][j]+h); 
   b0=f*ab[2][1]; ab[2]=f^2*ab[2]; ab[2][1]=b0; ab} 
addhelp(ab_shift, "ab_shift(ab,int1,{int2=[-1,1]}): shift three-term coefficients from interval int2 to interval int1")

{ab_kronrod(ab,n)=my(a,b,nmax,s,t,swap,sj,sk,j,k,l,m,k0,k1);  
  nmax=floor((length(ab[1])+length(ab[2]))/3);
  if(n==0, n=nmax, if(n>nmax, n=nmax)); a=vector(2*n+1); b=a;
  for(k=0,floor(3*n/2), a[k+1]=ab[1][k+1]);
  for(k=1,ceil(3*n/2+1), b[k]=ab[2][k]); b[1]=1;
  s=vector(floor(n/2)+2); t=s; t[2]=b[n+2]; 
  for(m=0,n-2, sk=0;
    forstep(k=floor((m+1)/2),0,-1, l=m-k;
      sk=sk+(a[k+n+2]-a[l+1])*t[k+2]+b[k+n+2]*s[k+1]-b[l+1]*s[k+2];
      s[k+2]=sk);
    swap=s; s=t; t=swap); 
  forstep(j=floor(n/2),0,-1, s[j+2]=s[j+1]); 
    for(m=n-1,2*n-3, k0=m+1-n; k1=floor((m-1)/2); sj=0;
    for(k=k0,k1, l=m-k; j=n-1-l;
      sj=sj-(a[k+n+2]-a[l+1])*t[j+2]-b[k+n+2]*s[j+2]+b[l+1]*s[j+3];
      s[j+2]=sj); 
    k=floor((m+1)/2); if(m%2, b[k+n+2]=s[j+2]/s[j+3],
      a[k+n+2]=a[k+1]+(s[j+2]-b[k+n+2]*s[j+3])/t[j+3]);
    swap=s; s=t; t=swap);
  a[2*n+1]=a[n]-b[2*n+1]*s[2]/t[2]; [a,b]}
{addhelp(ab_kronrod,
   "ab_kronrod(ab,{n=0}): recursion coefficients of Kronrod formula."
   " n=0 means largest possible n.")} 

{ab_wmul1(ab,z)=my(n,a,b,r,A,B,trunc,mu); a=ab[1]; b=drop(ab[2],1); 
  n=length(a); r=vector(n+1); trunc=n==length(b); 
  if(trunc,a=concat(a,0); n=n+1); A=drop(a,-1); B=b;
  r[1]=z-a[1]; r[2]=z-a[2]-b[1]/r[1]; A[1]=a[2]+r[2]-r[1]; mu=-r[1]*ab[2][1];
  for(k=2,n, if(k<n, r[k+1]=z-a[k+1]-b[k]/r[k]; A[k]=a[k+1]+r[k+1]-r[k]);
    B[k-1]=B[k-1]*r[k]/r[k-1]); 
  if(trunc,B=drop(B,-1)); [A,concat(mu,B)]}
{addhelp(ab_wmul1,
   "ab_wmul1(ab,{b=0}): modify ab for weight function w(x)*(x-b)")}

{ab_wdiv1(ab0,mu0,z)=my(n,r,nu,ab,aug); n=length(ab0[2]); 
  aug=length(ab0[1])==n;
  ab0[1]=concat(ab0[1],0); if(aug,ab0[2]=concat(ab0[2],0); n=n+1);
  ab=[vector(n-aug),vector(n)]; r=vector(n+1); nu=0; r[1]=-mu0;
  for(k=1,n, r[k+1]=z-ab0[1][k]-ab0[2][k]/r[k]);
  ab[1][1]=ab0[1][1]+r[2]; ab[2][1]=-r[1];
  for(k=2,n, if(k<=length(ab[1]), ab[1][k]=ab0[1][k]+r[k+1]-r[k]); 
  ab[2][k]=ab0[2][k-1]*r[k]/r[k-1]); ab[2][1]=mu0; ab}
{addhelp(ab_wdiv1,
   "ab_wdvi1(ab,mu0,{b=0}): modify ab for weight function w(x)/(x-b)")}

ab_radkron(ab,xn,n)=ab_radau(ab_wdiv1(ab_kronrod(ab_wmul1(ab,xn),n-1),ab[2][1],xn),xn) 

\\ Routines that create or modify quadrature formulas


\\ Direct conversions  ab <=> qe <=> dw <=> xw.     

{qe_ab(ab)=my(a,b,n); n=length(ab[1]);
  for(i=2,n,
   if(WARN_ROUNDOFF&ab[1][i-1]<0, print( "QE_AB: negative coefficient found."));
   if(ab[1][i-1]==0,error("QE_AB: zero divisor found."));
   ab[2][i]=ab[2][i]/ab[1][i-1]; ab[1][i]=ab[1][i]-ab[2][i]);
   if(length(ab[2])==n+1, ab[2][n+1]=ab[2][n+1]/ab[1][n]); 
   [ab[1],drop(ab[2],1),ab[2][1]]}
addhelp(qe_ab, "qe_ab(ab)");

{ab_qe(qe)=
  [vector(length(qe[1]),k,qe[1][k]+if(k>1,qe[2][k-1])), 
   concat(qe[3],vector(length(qe[2]),k,qe[2][k]*qe[1][k]))]}
addhelp(ab_qe, "ab_qe(qe)");

\\ Observation: the last weights are not accurate to full relative precision
{dw_qe(qe)=my(d,w,x,pu,nu,n); d=gaps(qe); n=length(d); dw=matrix(n,2); x=0;
  for(j=1,n, x=x+d[j]; dw[j,1]=d[j]; dw[j,2]=weight(qe,x));
  [0,dw]} 
addhelp(dw_qe, "dw_qe(qe)");

{qe_dw(dw)=my(d,w,u,v,n,mu0); dw[2][1,2]=dw[2][1,2]+dw[1]; dw=dw[2]; 
  n=matsize(dw)[1]; mu0=if(NORMALIZED, 1, sum(k=1,n,dw[k,2]));
  u=[[dw[n,2]],[]]; v=[dw[n,1],[]];
  for(k=1,n-1, u=prqd([concat(dw[n-k,2],v[1]),concat(u[1][1],v[2]),mu0]);
    v=stqd([concat(u[2],0),drop(u[1],1),mu0],-dw[n-k,1])); v } 
\\    concat(v,mu0)}
addhelp(qe_dw, "qe_dw(dw)");

{dw_xw(xw)=my(n); n=matsize(xw)[1];
  forstep(k=n,2,-1, xw[k,1]=xw[k,1]-xw[k-1,1]); [0,xw]}
addhelp(dw_xw, "dw_xw(xw)");

{xw_dw(dw)=my(n); dw[2][1,1]=dw[2][1,1]+dw[1]; dw=dw[2]; n=matsize(dw)[1];
  for(k=2,n, dw[k,1]=dw[k,1]+dw[k-1,1]); dw}
addhelp(xw_dw, "xw_dw(dw)");

{xw_ab(ab)=my(x,xw,n); x=qlpoles(ab); n=length(x); xw=matrix(n,2);
  DEBUG=1;
  for(k=1,n, xw[k,1]=x[k]; xw[k,2]=weight_ab(ab,x[k]));
\\ prql(ab,x[k]); xw[k,2]=LASTCS[1]*ab[2][1];
\\ if (precision(xw[k,2])<precision(0.), print("Precision loss"
\\ if(precision(LASTCS[1])<precision(0.)," in LASTCS"))));
  xw}

{ab_dw(dw)=my(h); h=dw[2][1,1]; dw[2][1,1]=0; 
  Shift(ab_qe(qe_dw(dw)),[0,h],[-h,0])}
    
\\ Shorthand for some indirect conversions
xw_qe(qe)=xw_dw(dw_qe(qe));
qe_xw(xw)=qe_dw(dw_xw(xw));
ab_xw(xw)=ab_dw(dw_xw(xw));  \\ Should be Gragg-Harrod
dw_ab(ab)=dw_xw(xw_ab(ab));  \\ Could be QL with explicit shifts

\\ Operations on continued and partial fractions
\\ stqd: f(z) -> f(z-sig)
\\ prqd: f(z) -> z*f(z-sig)-k (k determined by f^(infinity)=0)
\\ QD algorithms: Parlett, Acta Numerica 1995, p466, 468

{stqd(qe,sig)=my(q,e,n,t,u,f); q=qe[1]; e=qe[2]; n=length(q); t=-sig; 
 for(i=1,n-1, u=q[i]; q[i]=q[i]+t; f=e[i]/q[i]; e[i]=u*f; t=t*f-sig);
 q[n]=q[n]+t; [q,e,qe[3]]}
addhelp(stqd, "stqd(qe,{sig=0}): stationary qd algorithm");

{prqd(qe,sig)=my(q,e,n,c,d,f,t); q=qe[1]; e=qe[2]; n=length(q); t=-sig;
 c=qe[3]*q[1]; d=q[1]-sig;
 for(i=1,length(e), q[i]=d+e[i]; f=q[i+1]/q[i]; e[i]=e[i]*f; d=d*f-sig);
 if(length(e)<n,q[n]=d); [q,e,c]}
addhelp(prqd, "prqd(qe,{sig=0}): progressive qd algorithm");

{prql(T,sig,CS)=my(n,alf,gam,oldgam,oldC,BB,P,R,C,S,CSgiven);
  sig=precision(sig,precision(0.));
  n=length(T[1]); gam=T[1][n]-sig; oldC=1; BB=T[2][n]; CSgiven=CS!=[0,0];
  if(CSgiven, C=CS[1]; S=CS[2], P=gam*gam; R=P+BB; C=P/R; S=BB/R);
  forstep(i=n-1,1,-1,
    oldgam=gam; alf=T[1][i]; gam=precision(C*(alf-sig)-S*oldgam,precision(0.));    
\\    if(DEBUG,print(i" gam "gam));
    T[1][i+1]=oldgam+(alf-gam); P=if(C==0, oldC*BB, (gam*gam)/C);
\\    if(DEBUG,print("T[1][i+1] "T[1][i+1]));
    BB=T[2][i]; R=P+BB; oldC=C; 
    if(i>1, T[2][i+1]=S*R; C=P/R; S=BB/R));
  T[1][1]=sig+gam; T[2][2]=S*P; LASTCS=[C,S]; T}

{prqr(T,sig,CS)=my(n,alf,gam,oldgam,oldC,BB,P,R,C,S,CSgiven);
  n=length(T[1]); gam=T[1][1]-sig; oldC=1; BB=T[2][2]; CSgiven=CS!=[0,0];
  if(CSgiven, C=CS[1]; S=CS[2], P=gam*gam; R=P+BB; C=P/R; S=BB/R);
  for(i=2,n,
    oldgam=gam; alf=T[1][i]; gam=C*(alf-sig)-S*oldgam;
    T[1][i-1]=oldgam+(alf-gam); P=if(C==0, oldC*BB, (gam*gam)/C);
    if(i<n, BB=T[2][i+1]; R=P+BB; oldC=C; T[2][i]=S*R; C=P/R; S=BB/R));
  T[1][n]=sig+gam; T[2][n]=S*P; LASTCS=[C,S]; T}

{r_qe(qe,x)=my(f,n,c); q=qe[1]; e=qe[2]; n=length(q); c=qe[3]; r=0; 
  forstep(j=n,1,-1, r=if(j>1,e[j-1],c)/(x-q[j]/(1-r))); r}

{r_ab(ab,x)=my(r,n); a=ab[1]; b=ab[2]; n=length(a);
  r=0; forstep(j=n,1,-1, r=b[j]/(x-a[j]-r)); r}

{r_xw(xw,t)=sum(j=1,matsize(xw)[1], xw[j,2]/(t-xw[j,1]))}

{stieltjes_function(u,x)=my(typ); typ=Type(u);
  if(typ=="AB", r_ab(u,x),  if(typ=="QE", r_qe(u,x), if(typ=="XW", r_xw(u,x),
    if(typ=="DW", r_xw(xw_dw(u),x),
  error("Can't evaluate stieltjes_function for Type "typ)))))}

{o_pol(u,x)=my(typ); typ=Type(u);
  if(typ=="XW", u=ab_xw(u); typ="AB");
  if(typ=="DW", u=qe_dw(u); typ="QE");
  if(typ=="AB", return(p_ab(u,x)));
  if(typ=="QE", return(pu_qe(u,x)[1]));
  error("Can't evaluate orthogonal polynomials defined by Type "typ)}
 
{pu_qe(qe,x)=my(q,e,p,u,n,c);  
  q=qe[1]; e=qe[2]; n=length(q); u=vector(n+1); p=u; u[1]=1; p[1]=1; 
  for(k=1,n, u[k+1] = x*p[k]-if(k>1,e[k-1]*u[k]);
             p[k+1] = u[k+1]-q[k]*p[k]);  [p,u]}

{p_ab(ab,x)=my(a,b,n,p); a=ab[1]; b=ab[2]; n=length(a); p=vector(n+1); 
  p[1]=1; for(k=1,n, p[k+1]=(x-a[k])*p[k]-if(k>1,b[k]*p[k-1])); p}
addhelp(p_ab, "p_ab(ab,x): evaluate polynomials p_k(x)");

{qe_r(r)=my(n,q,e,c,u,p); u=numerator(r); p=denominator(r); n=length(p)-1; 
  q=vector(n); e=vector(n-1); c=polcoeff(u,n-1)/polcoeff(p,n); p=x*u-c*p;
  for(j=1,n-1, q[j]=polcoeff(p,n-j)/polcoeff(u,n-j); u=p-q[j]*u;
    e[j]=polcoeff(u,n-j-1)/polcoeff(p,n-j); p=x*u-e[j]*p); q[n]=p/u; [q,e,c]}

\\ Operations involving modified moments.  
\\ In all of these, the inner product is the one associated with qe, not QE.  
\\ Relative to DPL(1999) p140, the equations differ because we do not use any
\\ origin 0 vectors or indices here.  Thus, the indices on q, e, Q (=psi) and 
\\ E (=eta) are one less.  All modified moments are NORMALIZED so that mu[1]=1.

{mu_qe(qe,Z,flag)=my(Q,E,N,q,e,n,mu,ro,d,m,s);
  Q=Z[1]; E=Z[2]; N=length(Q);  q=qe[1]; e=qe[2]; n=length(q);
  mu=if(flag, matrix(N+1,N+1), vector(N+1)); ro=[qe[3]]; d=0; s=1;
  for(l=1,N+1, m=s; s=min(m+1,N-l+2); nu=vector(s,k,
      if(k<=m,ro[k]) + if(k>1,q[k-1]*ro[k-1]) - if(k<=d,Q[l-1]*nu[k]));
    if(flag, for(j=1,s, mu[j,l]=nu[j]), mu[l]=nu[1]); if (l==N+1, break);
    d=s; s=vecmin([d,n,N-l+1]); ro=vector(s,k,
      if(k<d,nu[k+1]) + if(k>1,e[k-1]*nu[k]) - if(k<=m&l>1,E[l-1]*ro[k]))); mu} 

{qe_mu(nu,Z,n)=my(q,e,Q,E,m,ro,nu1); Q=Z[1]; E=Z[2];
  if(n==0, n=floor(min(length(nu),length(Q)+1)/2)); 
  q=vector(n); e=vector(n-1); m=2*n; nu1=nu[1]; nu=concat(0,nu); 
  for(k=1,n, 
    ro=vector(m,l, nu[l+1] + if(k+l>2,nu[l]*Q[k+l-2]) - if(k>1,q[k-1]*ro[l+1]));
    if(k>1, e[k-1]=ro[1]/nu[1]); m=m-1; 
    nu=vector(m,l, ro[l+1] + if(k+l>2,ro[l]*E[k+l-2]) - if(k>1,e[k-1]*nu[l+1]));
    q[k]=nu[1]/ro[1]; m=m-1); [q,e,nu1]}

{mu_ab="Not implemented: mu_ab"}

{mu_xw(xw,R)=my(p,pk); 
  for(k=1,matsize(xw)[1], p=if(k>1,p)+xw[k,2]*o_pol(R,xw[k,1])); p}

{fs_xw(xw)=my(n); n=matsize(xw)[1];
    for(k=1,n-1, xw[k,2]=sum(j=k,n,xw[j,2]); 
        for(j=k+1,n, xw[j,2]*=(xw[j,1]-xw[k,1]))); xw}
{r_fs(xw,t)=my(f); f=0;
   forstep(k=matsize(xw)[1],1,-1, f=(f+xw[k,2])/(t-xw[k,1])); f }


\\ Calculation of nodes and weights 

poles(qe)=my(x); x=gaps(qe); for(k=2,length(x), x[k]=x[k]+x[k-1]); x

\\ Shift strategy: Wilkinson shift slightly decreased to reduce the risk
\\ of overshoot
{qdshift(qe,force)=my(q1,q2,e1); q1=qe[1][1]; q2=qe[1][2]; e1=qe[2][1];
  if(force||FOOLHARDY, poles(qe)[1],  
  real(poles(qe)[1])/(1.+abs(2*e1)/(abs(q1)+abs(q2))))}

{gaps(qe)=my(q,e,n,x,h,qm,s1,s,u,v,qditer,overshoot); 
  q=qe[1]; n=length(q); e=qe[2]=take(qe[2],n-1); 
  if(n==1, return([q[1]]*1.));
  if(n==2, x=(q[1]+e[1]+q[2])*0.5; qm=(q[1]+q[2])*0.5; 
     h=sqrt(((q[1]-q[2])^2+e[1]^2)/4+e[1]*qm); 
     if(abs(x+h)<abs(x-h),h=-h); return([q[1]*q[2]/(x+h),2*h]));
  qditer=4; s=0.; x=[];
  for(k=1,n-1, if(q[k]<=0||e[k]<=0, 
    if(!FOOLHARDY, error(
    "Some coefficients are not positive."  
    "\nTo continue regardless, set 'FOOLHARDY=1' and try again."
  )))); 
  while(n>2, 
    s1=if(qditer>0, qditer=qditer-1; 0., qdshift(take(qe,-2))); 
    qe=prqd(qe,s1); s=s+s1;  overshoot=real(s1*qe[1][n])<0;
    if(abs(qe[1][n])+qe[1][1]==qe[1][1]&abs(qe[2][n-1])+qe[2][1]==qe[2][1],
      x=concat(x,s); n=n-1; qe=drop(qe,-1); s=0, \\ else
      if(overshoot&WARN_SHIFT, print(qe); print("Overshoot in 'gaps' for n=",n,
      "; shift=", qe[1][n], ", factor=",precision(qe[1][n]/s1,9)))));
  concat(x,gaps(qe))}

{weight(qe,x)=my(n,pu,nsq,w); n=length(qe[1]); 
  nsq=1; w=1; pu=o_pol(qe,x); 
  for(j=1,n-1, nsq=nsq*qe[1][j]*qe[2][j]; w=w+pu[j+1]^2/nsq);
  qe[3]/w}

{weight_ab(ab,x)=my(n,pu,nsq,w); n=length(ab[1]); 
  nsq=1; w=1; pu=o_pol(ab,x); 
  for(j=1,n-1, nsq=nsq*ab[2][j+1]; w=w+pu[j+1]^2/nsq);
  ab[2][1]/w}

{xw_combine(F1,F2,eps,g)=my(j,m,n,x,w,X,W); 
  x=F1[1]; w=F1[2]; X=F2[1]; W=F2[2]; m=length(x)+length(X);
  if(!eps, eps=eps=2*m^2*
     10.0^(-vecmin([precision(F1),precision(F2),precision(0.)]))); 
  if(!g, g=[1/2,1/2]); F1=[concat(x,X),concat(g[1]*w,g[2]*W)];
  n=vecsort(F1[1],,1); x=vecextract(F1[1],n); w=vecextract(F1[2],n); n=m; 
  for(k=1,n-1, if(abs(x[k]-x[k+1])<eps, w[k]=w[k]+w[k+1]; 
    if(abs(w[k])<eps, w[k]=0; m=m-1); w[k+1]=0; m=m-1));
  X=vector(m); W=vector(m); j=0;
  for(k=1,n, if(w[k], j=j+1; X[j]=x[k]; W[j]=w[k])); [X,W]
  }

{xw_shift(xw,int1,int2)=my(f,h); if(int2==0, int2=[-1,1]);
   f=(int1[2]-int1[1])/(int2[2]-int2[1]); h=int1[1]-f*int2[1];
   xw=xw*f; for(j=1,matsize(xw)[1],xw[j,1]=xw[j,1]+h); xw}
addhelp(xw_shift, "xw_shift(xw,int1,{int2=[-1,1]}): shift formula from interval int2 to interval int1")

{dw_shift(dw,int1,int2)=my(f,h); if(int2==0, int2=[0,2]);
   f=(int1[2]-int1[1])/(int2[2]-int2[1]); h=int1[1]-f*int2[1];
   dw=dw*f; dw[1]=dw[1]+h; dw}
addhelp(xw_shift, "xw_shift(xw,int1,{int2=[0,2]}): shift formula from interval int2 to interval int1")

{qlpoles(T)=my(n,d,a1,a2,b,iter,x);
  n=length(T[1]); T[2]=take(T[2],n); a1=T[1][1]; x=[]; if(n==1, return([a1]));
  b=T[2][2]; a2=T[1][2];
  if(n==2, d=(a2-a1)/2; h=b/(d+if(d<0,-1,1)*sqrt(d^2+b)); return([a1-h,a2+h]));
  while(n>2, d=abs(T[1][1])+abs(T[1][2]);
    while(abs(T[2][2])+d>d, sig=qlpoles(AB(T,3))[1]; T=prql(T,sig);
      d=abs(T[1][1])+abs(T[1][2]) );
    x=concat(x,T[1][1]); n=n-1; T=[take(T[1],-n),take(T[2],-n)]);
  concat(x,qlpoles(T)) }

\\ Operations on orthogonal expansions
{check_length(ab,n,name,dif)=if(length(ab[1])<n, 
  error(name ": length of ab (" length(ab[1]) ") should be at least " n))}
{o_mulx(c,ab)=my(n); n=length(c); check_length(ab,n,"o_mulx"); 
  vector(n+1,l, 
    if(l>1,c[l-1]) + if(l<n,ab[2][l]*c[l+1]) + if(l<=n,ab[1][l]*c[l]))}
{o_ser(c,ab,x)=my(n,f,f0,f1);  n=length(c); check_length(ab,n-1,"o_ser");
  f1=0; f=c[n]; forstep(k=n-1,1,-1, 
    f0=f1; f1=f; f=(x-ab[1][k])*f1-if(k<n-1,ab[2][k]*f0)+c[k]); 
  f}
{o_dot(c1,c2,nu)=my(n); n=min(length(c1),length(c2)); 
  check_length(ab,n,"o_dot"); sum(i=1,n, c1[i]*c2[i]*nu[i])}
 
\\ Recursion coefficients for refinable functions
{ab_refinable(n,g)=my(ab,ckk,p0,p1,p2,nu,M,N,symm); M=length(g); 
 nu=vector(n); ab=[vector(n),vector(n-1)]; g=2*g/sum(j=1,M, g[j]);
 symm=1; for(i=1,floor(M/2), if(g[i]!=g[M+1-i], symm=0; break));
 N=if(symm, N=floor((M+1)/2), M); 
 p1=vector(N,j,[]); p2=vector(N,j,[1]); ckk=1; 
 for(k=1,n, p0=p1; p1=p2;
   if(k==1,  nu[k]=1, ;
     nu[k]=sum(j=1,N, if(symm&(2*j<=M), 2, 1)* 
       g[j]*o_dot(p1[j],p1[j],nu))/(2*(1-ckk^2));
     ab[2][k-1]=nu[k]/nu[k-1] ); 
   for(j=1,N, p2[j]=o_mulx(p1[j],ab)+(j-1)*concat(p1[j],0) );  
   ckk=ckk/2; ab[1][k]=if(symm, (M-1)/2,
     sum(j=1,N, g[j]*o_dot(p1[j],p2[j],nu))/(4*nu[k]*(1-2*ckk^2)));
   if(k==n, break);
   for(j=1,N, p2[j][k]=p2[j][k]+2*ckk*ab[1][k]; p2[j] = p2[j]/2 - 
       ab[1][k]*concat(p1[j],0) - if(k>1,ab[2][k-1]*concat(p0[j],[0,0])) ) ); 
 ab[2]=concat(1,ab[2]); ab} 


\\ Utilities

{drop(x,k)=my(n); n=length(x); 
  if(type(x)!="t_VEC", return(x));
  if(type(x[1])=="t_VEC", for(j=1,n,x[j]=drop(x[j],k)); return(x));
  if(k>=0, vector(n-k,j,x[j+k]), vector(n+k,j,x[j]))}

{take(x,k)=my(n); n=length(x); 
  if(type(x)!="t_VEC", return(x));
  if(type(x[1])=="t_VEC", 
     for(j=1,n,x[j]=take(x[j],k-(j-1)*sign(k))); return(x));
  if(k>=0, vector(k,j,x[j]), vector(-k,j,x[n+k+j]))}

{reverse(x)=my(n); n=length(x);
  if(type(x)!="t_VEC", return(x));
  if(type(x[1])=="t_VEC", for(j=1,n,x[j]=reverse(x[j])); return(x));
  vector(n,j,x[n+1-j])}

fac(x)=if(type(x)=="t_INT", x!, gamma(1+x));
addhelp(fac, "fac(x): x!, taking into account whether x is integer") 

inside(x,interval)= x>=interval[1] & x<=interval[2];

\\ Various
{exactness(Z1,Z2,eps)=my(n,N,d,j);
  n=length(Z1[1]); N=length(Z2[1]); d=0; j=1;
  if(eps==0, eps=2*(min(n,N)^2)*
     10.0^(-vecmin([precision(Z1),precision(Z2),precision(0.)]))); 
  if(abs(Z1[3]-Z2[3])>eps, return(-1));
   while(1, if(j>n||j>N||abs(Z1[1][j]-Z2[1][j])>eps, break); d=d+1;
    if(j>=n||j>=N||abs(Z1[2][j]-Z2[2][j])>eps, break); d=d+1; j=j+1); d}

vxw_vqe(vqe)=vector(length(vqe),n,xw_qe(vqe[n]))

integrate(fct,xw)=my(s); for(j=1,matsize(xw)[1], x=xw[j,1]; s=s+xw[j,2]*eval(fct)); s

{stratified_integrate(vxw)=my(m,y); m=length(vxw); y=vector(m);
  for(k=1,m, y[k]=integrate(vxw[k]); if(k>1, y[k]=(y[k]+y[k-1])/2)); y}

\\ {qe_wmul1(qe,b)= qe=prqd(qe,b); qe[1]=drop(qe[1],-1); 
\\  print("qe_wmul1 is not debugged"); qe}

{findpolynomial(u)=my(pk); for(k=1,length(SHORT_POLYNOMIALS), 
  pk=SHORT_POLYNOMIALS[k]; 
  if(u>=pk & u<=concat(pk,"}"), return(pk)))}
  
\\ AB, QE, XW, DW: There is a lot of duplication, which might be
\\ avoidable by clever programming.

{AB(u,deg,p1,p2)=my(typ,m,n); typ=Type(u); 
  if (typ=="t_INT", deg=u; u="Leg"; typ="t_STR");
  if (typ=="t_STR", typ=findpolynomial(u);
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", p1=0; p2=0; typ="Jac");
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", return(ab_shift(AB(qe_jacobi(deg,p1,p2)),[-1,1],[0,2])));
    if (typ=="Lag", return(AB(qe_laguerre(deg,p1))));
    if (typ=="Her", return(ab_hermite(deg,p1)));
    if (typ=="Ref", return(ab_refinable(deg,p1))));
  if (deg==0 & typ=="AB", return(u));
  if (typ=="AB", return(Take(u,deg+1,typ)));
  if (typ=="QE", return(ab_qe(u)));
  if (typ=="XW", return(ab_xw(u)));
  if (typ=="DW", return(ab_dw(u)))
}

{QE(u,deg,p1,p2)=my(typ,m,n); typ=Type(u);
  if (typ=="t_INT", deg=u; u="Leg"; typ="t_STR");
  if (typ=="t_STR", typ=findpolynomial(u); 
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", p1=0; p2=0; typ="Jac");
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", return(qe_jacobi(deg,p1,p2)));
    if (typ=="Lag", return(qe_laguerre(deg,p1)));
    if (typ=="Her", error("Not implemented: QE Hermite"));
    if (typ=="Ref", return(QE(ab_refinable(dim,p1)))));
  if (deg==0 & typ=="QE", return(u)); 
  if (typ=="QE", return(Take(u,deg+1,typ)));
  if (typ=="AB", return(qe_ab(u)));
  if (typ=="XW", return(qe_xw(u)));
  if (typ=="DW", return(qe_dw(u)))
}

{XW(u,npts,p1,p2)=my(typ,m,n,deg); typ=Type(u); 
  if (typ=="t_INT", npts=u; u="Leg"; typ="t_STR");
  deg=2*npts-1;
  if (typ=="t_STR", typ=findpolynomial(u); 
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", p1=0; p2=0; typ="Jac");
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", 
      return(XW(dw_shift(dw_qe(qe_jacobi(deg,p1,p2)),[-1,1],[0,2]))));
    if (typ=="Lag", return(XW(qe_laguerre(deg,p1))));
    if (typ=="Her", return(XW(ab_hermite(deg,p1))));
    if (typ=="Ref", return(XW(ab_refinable(npts,p1)))));
  if (deg==0 & typ=="XW", return(u)); 
  if (typ=="XW", 
    return(vecextract(u,vector(npts,k,k),vector(matsize(u)[2],k,k))));
  if (typ=="AB", return(xw_ab(u)));
  if (typ=="QE", return(xw_dw(dw_qe(u))));
  if (typ=="DW", return(xw_dw(u)))
}

{DW(u,npts,p1,p2)=my(typ,m,n,deg); typ=Type(u);
  if (typ=="t_INT", npts=u; u="Leg"; typ="t_STR");
  deg=2*npts-1;
  if (typ=="t_STR", typ=findpolynomial(u); 
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", p1=0; p2=0; typ="Jac");
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", return(dw_qe(qe_jacobi(deg,p1,p2))));
    if (typ=="Lag", return(dw_qe(qe_laguerre(deg,p1))));
    if (typ=="Her", return(dw_xw(xw_ab(ab_hermite(deg,p1)))));
    if (typ=="Ref", return(dw_xw(xw_ab(ab_refinable(npts,p1))))));
  if (deg==0 & typ=="DW", return(u)); 
  if (typ=="DW", return([u[1],XW(u[2],npts)]));
  if (typ=="AB", return(dw_xw(xw_ab(u))));
  if (typ=="QE", return(dw_qe(u)));
  if (typ=="XW", return(dw_xw(u)))
}

{MU(U,R)=my(typu,typr,mu); typu=Type(U); typr=Type(R);
  if(typu=="QE" & typr=="QE", return(mu_qe(U,R)));
  if(typu=="DW", U=XW(U); typu="XW");
  if(typu=="XW", return(mu_xw(U,R)))} 

{Shift(u,I1,I2)=my(typ); typ=Type(u);
  if (typ=="AB", return(ab_shift(u,I1,I2)));
  if (typ=="XW", return(xw_shift(u,I1,I2)));
  if (typ=="DW", return(dw_shift(u,I1,I2)));
  if (typ=="MU", error("Not implemented: Shift MU"));
  if (typ=="QE", return(qe_shift(u,I1,I2)))}

{Dimension(u,typ)= if (typ==0, typ=Type(u));
  if(typ=="MU", return(length(u)));
  if(typ=="XW", return(2*matsize(u)[1]));
  if(typ=="DW", return(2*matsize(u[2])[1]));
  if(typ=="AB", return(length(u[1])+length(u[2])));
  if(typ=="QE", return(length(u[1])+length(u[2])+1))}
addhelp(Dimension,"Dimension(u,{typ=Type(u)}): number of significant elements in u.")

{Take(R,dim,typ)=my(n); if(dim==0, return(R));
  if (typ==0, typ=Type(R));
  if (typ=="MU", return(take(R,dim)), n=sign(dim)*floor(abs(dim)/2));
  if (typ=="QE", R[1]=take(R[1],n); R[2]=take(R[2],dim-n-sign(dim)),
    if (typ=="AB", R[1]=take(R[1],n); R[2]=take(R[2],dim-n),
      error("Can't Take from Type "typ)));
     R}

{Type(u)=my(typ,typ2,len,rtyp); rtyp=0; typ=type(u); len=length(u);
  if (typ=="t_MAT", rtyp="XW",
    if (typ!="t_VEC", rtyp="?",
      if(len<2, rtyp="MU")));
  if (rtyp==0,
    typ2=type(u[2]);
    if (type(u[1])!="t_VEC", rtyp=if(len==2 & typ2=="t_MAT", "DW", "MU"),
    if (length(u)==2 & typ2=="t_VEC", rtyp="AB",
    if (length(u)==3 & typ2=="t_VEC" & type(u[3])!="t_VEC", rtyp="QE",
    rtyp="MU"))));
  Check(u,rtyp); if(rtyp=="?", typ, rtyp)}

{Check(u,typ,force)=my(OK); if (force+DEBUG==0, return);
  OK=1;
  if(OK&typ=="XW", OK=inside(matsize(u)[2],[2,3]));
  if(OK&typ=="DW", Check(u[2],"XW"));
  if(OK&typ=="AB", OK=inside(length(u[2])-length(u[1]),[0,1]));
  if(OK&typ=="QE", OK=inside(length(u[1])-length(u[2]),[0,1]));
  if(OK&typ=="MU", for(j=1,length(u), if(type(u[j])=="t_VEC", OK=0; break)));
  if(!OK, error(u" does not have the right shape for Type "typ))}

{ print1("   <PORTHOS "floor(PORTHOS_VERSION/100)"."
   floor((PORTHOS_VERSION%100)/10)""PORTHOS_VERSION%10">"); 
  read("~/pari/porthos-help.gp"); }

ab_hermite(deg,p)=my(m); m=ceil(deg/2); [vector(m,n,0),concat(sqrt(Pi),vector(deg-m,n,n/2))];


