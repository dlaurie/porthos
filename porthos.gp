\\                              PORTHOS

\\ A package of Pari/GP routines for working with orthogonal polynomials
\\   Dirk Laurie   dlaurie@na-net.ornl.gov   

\\ Quick help: inside Pari/GP, type 'help'.  Everything useful in here
\\   can be found there too. If all else fails: read the TeX documentation.
\\   TODO: write the TeX documentation.

\\ !! Development notes:
\\   1. If a routine has an "addhelp" in "porthos-help.gp", that routine 
\\      is considered to be stable. 
\\   2. If a routine has an "addhelp" here only, that routine has been 
\\      briefly tested during an interactive session.
\\   3. If a routine has no "addhelp", it is either an oversight, or the
\\      routine is untested.
\\   3. The TeX documentation does not exist yet.

PORTHOS_VERSION=60;
\\ 60: 1 May 2016, incorporates GP language changes after 2009: 
\\    '||' for '|', 'my' for 'local'
\\ 54: improvements to help only
\\ 53: 4 December 2008

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

incompatible=0;
gp_version = subst(Pol(version),x,10)
if (gp_version<240, if(PORTHOS_VERSION>=60, incompatible="old"))
if (gp_version>=240, if(PORTHOS_VERSION<60, incompatible="new"))
{ if (incompatible, print("Your version of Pari-GP is too "incompatible
" for your Porthos.\n Both [Pari>=2.4,Porthos>=0.60] or neither must hold.")) } 

\\ QE: Routines that create or modify two-term coefficients

is_scalar(gen)= type(gen)!="t_VEC" && type(gen)!="t_MAT"
is_vector(gen,len)= type(gen)=="t_VEC" && #gen==len

{is_QE(gen)= is_vector(gen,3) && type(gen[1])=="t_VEC" && type(gen[2])=="t_VEC"
  && #gen[1]>=#gen[2] && #gen[1]<=#gen[2]+1 && is_scalar(gen[3]) }
  
{qe_sf(r)=my(n,q,e,c,u,p); u=numerator(r); p=denominator(r); n=length(p)-1; 
  q=vector(n); e=vector(n-1); c=polcoeff(u,n-1)/polcoeff(p,n); p='z*u-c*p;
  for(j=1,n-1, q[j]=polcoeff(p,n-j)/polcoeff(u,n-j); u=p-q[j]*u;
    e[j]=polcoeff(u,n-j-1)/polcoeff(p,n-j); p='z*u-e[j]*p); 
  q[n]=polcoeff(p/u,0); [q,e,c]}

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

{qe_modify(qe,modify)=
  if(!modify, return(qe));
  if(modify=="anti", return(qe_antigauss(qe)));
  if(modify=="radau", return(qe_radau(qe)));
  if(modify=="lobatto", return(qe_lobatto(qe,2)));
  error("Unknown modification '"modify"' in qe_modify")}

\\{qe_shift(qe,int1,int2)=my(f,h); if(int2==0, int2=[0,2]);
\\   f=(int1[2]-int1[1])/(int2[2]-int2[1]); h=int1[1]-f*int2[1];
\\   stqd(qe*f,-h)}
\\addhelp(qe_shift, "qe_shift(qe,int1,{int2=[0,2]}): shift two-term coefficients from interval int2 to interval int1")

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

\\ AB: Routines that create or modify three-term coefficients

{is_AB(gen)= is_vector(gen,2) && type(gen[1])=="t_VEC" && type(gen[2])=="t_VEC"
  && #gen[1]<=#gen[2] && #gen[1]>=#gen[2]-1 } 
  
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

ab_hermite(deg,p)=my(m); m=ceil(deg/2); [vector(m,n,0),concat(sqrt(Pi),vector(deg-m,n,n/2))];  

{ab_modify(ab,modify)=
  if(!modify, return(ab));
  if(modify=="anti", ab[#ab,2]*=2; return(ab));
  if(modify=="radau", return(ab_radau(ab)));
  error("Unknown modification '"modify"' in qe_modify")}
\\ Routines that create or modify quadrature formulas

is_XW(gen)= type(gen)=="t_MAT" && (#gen==2 || #gen==3)
is_DW(gen)= is_vector(gen,2) && is_scalar(gen[1]) && is_XW(gen[2])

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
{dw_qe(qe)=my(d,w,x,pu,nu,n);  
  d=gaps(qe); n=length(d); dw=matrix(n,2); x=0;
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

{xw_ab(ab)=my(x,xw,n); x=vecsort(qlpoles(ab)); n=length(x); xw=matrix(n,2);
  DEBUG=1;
  for(k=1,n, xw[k,1]=x[k]; xw[k,2]=weight_ab(ab,x[k]));
\\ prql(ab,x[k]); xw[k,2]=LASTCS[1]*ab[2][1];
\\ if (precision(xw[k,2])<precision(0.), print("Precision loss"
\\ if(precision(LASTCS[1])<precision(0.)," in LASTCS"))));
  xw}

{ab_dw(dw)=my(h); h=dw[2][1,1]; dw[2][1,1]=0; 
  Shift(ab_qe(qe_dw(dw)),-h)}
    
\\ Shorthand for some indirect conversions
\\xw_qe(qe)=xw_dw(dw_qe(qe));
\\qe_xw(xw)=qe_dw(dw_xw(xw));
\\ab_xw(xw)=ab_dw(dw_xw(xw));  \\ TODO: Should be Gragg-Harrod
\\dw_ab(ab)=dw_xw(xw_ab(ab));  \\ Could be QL with explicit shifts

\\ Operations on continued fractions
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

\\ Operations on Stieltjes functions

{is_SF(gen)= type(gen)=="t_RFRAC" && 
  poldegree(numerator(gen)) == poldegree(denominator(gen))-1 }

{sf_qe(qe)=my(f,n,c); q=qe[1]; e=qe[2]; n=length(q); c=qe[3]; r=0; 
  forstep(j=n,1,-1, r=if(j>1,e[j-1],c)/('z-q[j]/(1-r))); r}

{sf_ab(ab)=my(r,n); a=ab[1]; b=ab[2]; n=length(a);
  r=0; forstep(j=n,1,-1, r=b[j]/('z-a[j]-r)); r}

{sf_xw(xw)=sum(j=1,#xw~, xw[j,2]/('z-xw[j,1]))}

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

\\ Operations involving modified moments, type 'NU'. These are moments 
\\ 'nu' with respect to the classical Chebyshev polynomials (not monic,
\\ over [-1,1]), and moments 'mu' with respect to an arbitrary monic basis 
\\ generated by an object Z of Type QE.

{ is_NU(gen)= my(k,ok); k=1; ok=type(gen)=="t_VEC" && #gen[1]>0;
  while(ok&&k<=#gen, ok = is_scalar(gen[k]); k=k+1); ok }

\\ moments relative to classical Chebyshev polynomials
{nu_xw(xw,n)= my(nu,T0,T1,T2,x,w);
  if(n==0,n=Dimension(xw));
  x=xw[,1]; w=xw[,2]~;
  nu=vector(n); T0=vector(#x,j, 1)~; T1=x;
  nu[1]=w*T0; nu[2]=w*T1;
  for(k=3,n, T2=vector(#x,j, 2*x[j]*T1[j]-T0[j])~; nu[k]=w*T2;
    T0=T1; T1=T2);
  nu }

\\ Recursion formulas from Chebyshev moments as in JCAM(5)1979-p235-243.pdf
{ ab_nu(nu,n) = my(a,b,nrm,old,p0,p1,p2,q,r,
    \\ multiply T-series by Tk
    tkpk(a) = my(b); b=vector(2*#a-1); b[#a]=a[1];
      for(i=1,#a-1,b[#a+i]=a[i+1]/2;
      b[#a-i]=a[i+1]/2); b,
    \\ multiply T-series by T1
    t1pk(a) = my(b); b=vector(#a+1); 
      b[2]=a[1]; if (#a==1, b[1]=0; return(b));
      b[1]=a[2]/2; for (i=2,#a, b[i+1]+=a[i]/2);
      for (i=3,#a, b[i-1]+=a[i]/2); b,
    \\ scalar product of x with a possibly offset part of y
    scalp(x,y,d) = sum(j=1,#x, x[j]*y[j+d]));
  if (n==0,n=floor(#nu/2)); a=vector(n); b=vector(n); 
  p1=[1]; b[1]=nu[1];
  for (k=1,n,
    q=tkpk(p1); p2=t1pk(p1);  
    nrm = scalp(q,nu); 
    r=0; if (k>1, b[k]=nrm/old; r=p1[k-1]/p1[k]);
    if(k==2, r=2*r);
    a[k] = (scalp(q,nu,1)/nrm+r)/2;
    for(i=1,#p1, p2[i]-=a[k]*p1[i]);
    if (k>2, b[k]=b[k]/2); 
    for(i=1,#p0, p2[i]-=b[k]*p0[i]);
    p0=p1; p1=p2; old=nrm;
  );
  [a,b]
}

nu_qe(qe)=my(k=1,mu=mu_qe(qe,"Che")); for(j=3,#mu, k=k*2; mu[j]=k*mu[j]); mu
  
\\ Mixed moments.
\\ In mu_qe, qe_mu and mu_xw, the inner product is the one associated 
\\ with qe, not Z, which is used purely to generate the basis polynomials. 
\\ Relative to DPL(1999) p140, the equations differ because we do not use any
\\ origin 0 vectors or indices here.  Thus, the indices on q, e, Q (=psi) and 
\\ E (=eta) are one less.  All modified moments are NORMALIZED so that mu[1]=1.

{mixed(qe,Z,flag)=my(Q,E,N,q,e,n,mu,ro,d,m,s);
  if(!Z, return(nu_qe(qe)));
  Q=Z[1]; E=Z[2]; N=length(Q);  q=qe[1]; e=qe[2]; n=length(q);
  mu=if(flag, matrix(N+1,N+1), vector(N+1)); ro=[qe[3]]; d=0; s=1;
  for(l=1,N+1, m=s; s=min(m+1,N-l+2); nu=vector(s,k,
      if(k<=m,ro[k]) + if(k>1,q[k-1]*ro[k-1]) - if(k<=d,Q[l-1]*nu[k]));
    if(flag, for(j=1,s, mu[j,l]=nu[j]), mu[l]=nu[1]); if (l==N+1, break);
    d=s; s=vecmin([d,n,N-l+1]); ro=vector(s,k,
      if(k<d,nu[k+1]) + if(k>1,e[k-1]*nu[k]) - if(k<=m&l>1,E[l-1]*ro[k]))); mu}
{addhelp(mixed, 
"mixed(qe,Z,flag): Mixed moments of monic orthogonal polynomials generated"
"\nby qe and Z, in the inner product of qe. flag=0 means only the first"
"\nrow, anything else means the full matrix");}

mu_qe=mixed

{qe_mu(nu,Z,n)=my(q,e,Q,E,m,ro,nu1);
  Q=Z[1]; E=Z[2];
  if(n==0, n=floor(min(length(nu),length(Q)+1)/2)); 
  q=vector(n); e=vector(n-1); m=2*n; nu1=nu[1]; nu=concat(0,nu); 
  for(k=1,n, 
    ro=vector(m,l, nu[l+1] + if(k+l>2,nu[l]*Q[k+l-2]) - if(k>1,q[k-1]*ro[l+1]));
    if(k>1, e[k-1]=ro[1]/nu[1]); m=m-1; 
    nu=vector(m,l, ro[l+1] + if(k+l>2,ro[l]*E[k+l-2]) - if(k>1,e[k-1]*nu[l+1]));
    q[k]=nu[1]/ro[1]; m=m-1); [q,e,nu1]}
addhelp(qe_mu, "qe_mu(nu,Z,n)");

{mu_xw(xw,R)=my(p,pk);
  if(type(R)=="t_INT", return(nu_xw(xw)));
  for(k=1,matsize(xw)[1], p=if(k>1,p)+xw[k,2]*o_pol(R,xw[k,1])); p}
addhelp(mu_xw, "mu_xw(xw,R)");

\\ Operations involving factorial series TODO: WHY IS THIS HERE?

is_FS(gen)= type(gen)=="t_MAT" && matsize(gen)[1]==2

{fs_xw(xw)=my(n); n=matsize(xw)[1];
    for(k=1,n-1, xw[k,2]=sum(j=k,n,xw[j,2]); 
        for(j=k+1,n, xw[j,2]*=(xw[j,1]-xw[k,1]))); xw}
addhelp(fs_xw, "fs_xw(xw)");

{xw_fs(fs)=my(n); n=matsize(fs)[1];} \\ TODO

{sf_fs(fs,t)=my(f); f=0;
   forstep(k=matsize(fs)[1],1,-1, f=(f+fs[k,2])/(t-fs[k,1])); f }
addhelp(sf_fs, "sf_fs(fs)");

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

{qlpoles(T)=my(n,d,a1,a2,b,iter=0,x);
  n=length(T[1]); T[2]=take(T[2],n); a1=T[1][1]; x=[]; if(n==1, return([a1]));
  b=T[2][2]; a2=T[1][2];
  if(n==2, d=(a2-a1)/2; h=b/(d+if(d<0,-1,1)*sqrt(d^2+b)); return([a1-h,a2+h]));
  while(n>2, d=abs(T[1][1])+abs(T[1][2]); 
    while(abs(T[2][2])+d>d || iter<10, sig=qlpoles(Take(T,4))[1]; T=prql(T,sig);
      d=abs(T[1][1])+abs(T[1][2]); iter=iter+1 );
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
addhelp(drop,"drop(x,k): drop first k or last -k elements of vector x, or of each element if x is a vector of vectors")

{take(x,k)=my(n); n=length(x); 
  if(type(x)!="t_VEC", return(x));
  if(type(x[1])=="t_VEC", 
     for(j=1,n,x[j]=take(x[j],k-(j-1)*sign(k))); return(x));
  if(k>=0, vector(k,j,x[j]), vector(-k,j,x[n+k+j]))}
addhelp(take,"take(x,k): keep only first k or last -k elements of vector x, or of each element if x is a vector of vectors")

{reverse(x)=my(n); n=length(x);
  if(type(x)!="t_VEC", return(x));
  if(type(x[1])=="t_VEC", for(j=1,n,x[j]=reverse(x[j])); return(x));
  vector(n,j,x[n+1-j])}
addhelp(reverse,"reverse(x): reverse the vector x, or each element if x is a vector of vectors")

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
addhelp(exactness,"exactness(Z1,Z2,eps): degree of polynomial to which Z1 and Z2 agree to optional tolerance");

vxw_vqe(vqe)=vector(length(vqe),n,xw_qe(vqe[n]))

{integrate(fct,xw)= if(type(fct)!="t_CLOSURE", 
  error("first parameter to 'integrate' must be a t_CLOSURE")); 
  sum(j=1,matsize(xw)[1], xw[j,2]*fct(xw[j,1])) }

{stratified_integrate(vxw)=my(m,y); m=length(vxw); y=vector(m);
  for(k=1,m, y[k]=integrate(vxw[k]); if(k>1, y[k]=(y[k]+y[k-1])/2)); y}

\\ {qe_wmul1(qe,b)= qe=prqd(qe,b); qe[1]=drop(qe[1],-1); 
\\  print("qe_wmul1 is not debugged"); qe}

{findpolynomial(u)=my(pk); for(k=1,length(SHORT_POLYNOMIALS), 
  pk=SHORT_POLYNOMIALS[k]; 
  if(u>=pk & u<=concat(pk,"}"), return(pk)))}
  
\\ --- Moment expansion from weight function ---
\\ Only actually works if w = (x)->1 because Pari-GP can't do a definite
\\ integral of a power series
\\ Will need something like Lewanowicz's algorithm if we are to get
\\ anything useful out of knowing w.
{mx_fct(w,n) = my(serp=default(seriesprecision),s);
  if(n,default(seriesprecision,n),n=serp);
\\ we want the expansion
\\ intformal(w(x)/(z-x),'x)) = mu0/z + mu1/z^2 + ... + mun/z^n + O(z^-(n+1))
\\ but GP (a) can't do negative power asymptotics (b) will not expand the
\\ quotient.
\\ rewrite as
\\ intformal(t*w(x)/(1-x*t),'x)) = mu0*t + mu1*t^2 + ... + mun*t^n + O(t^(n+1))
  s=sum(k=0,n-1,('x*'t)^k)+O('t^(n+1));
  s='t*intformal(w('x)*s,'x);
  default(seriesprecision,serp); 
  subst(s,x,1) }

\\ AB, QE, XW, DW: There is a lot of duplication, which might be
\\ avoidable by clever programming.


\\ -------------- High-level frontend ---------------

{AB(u,deg,p1,p2)=my(typ,m,n); typ=Type(u); 
  if (is_AB(u), return(u));
  if (is_XW(u), u=dw_xw(u));    
  if (is_SF(u), u=qe_sf(u)); \\ TODO: implement ab_sf
  if (is_DW(u), u=qe_dw(u));
  if (is_QE(u), return(ab_qe(u)));
  if (typ=="t_INT", deg=u; u="Leg"; typ="t_STR");
  if (typ=="t_STR", typ=findpolynomial(u);
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", p1=0; p2=0; typ="Jac");
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", return(ab_shift(ab_qe(qe_jacobi(deg,p1,p2)),[-1,1],[0,2])));
    if (typ=="Lag", return(ab_qe(qe_laguerre(deg,p1))));
    if (typ=="Her", return(ab_hermite(deg,p1)));
    if (typ=="Ref", return(ab_refinable(deg,p1))));
}

{QE(u,deg,p1,p2)=my(typ,m,n); typ=Type(u);
  if (is_QE(u), return(u));
  if (is_AB(u), return(qe_ab(u)));
  if (is_XW(u), u=dw_xw(u));
  if (is_DW(u), return(qe_dw(u)));
  if (is_SF(u), return(qe_sf(u)));
  if (is_NU(u), return(qe_mu(u)));
  if (typ=="t_INT", deg=u; u="Leg"; typ="t_STR");
  if (typ=="t_STR", typ=findpolynomial(u); 
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", p1=0; p2=0; typ="Jac");
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", return(qe_jacobi(deg,p1,p2)));
    if (typ=="Lag", return(qe_laguerre(deg,p1)));
    if (typ=="Her", error("Not implemented: QE Hermite"));
    if (typ=="Ref", return(QE(ab_refinable(dim,p1)))));
}

{XW(u,npts,p1,p2,modify)=my(deg,shifted=0,family);
print("Type = ",Type(u));
  if (is_XW(u), return(u));
  if (is_DW(u), return(xw_dw(u)));
  if (is_NU(u), u=ab_nu(u));
  if (is_SF(u), u=qe_sf(u));
\\ --- convert string or numeric arguments to AB or QE
  if (type(modify)!="t_STR",
    if (type(p2)=="t_STR", modify=p2; p2=0,
      if (type(p1)=="t_STR", modify=p1; p1=0,
        if (type(npts)=="t_STR", modify=npts))));
  if (type(u)=="t_INT", npts=u; u="Leg");
  deg=2*npts-1;
  if (type(u)=="t_STR", family=findpolynomial(u); 
    if (family=="Che", p1=-1/2; p2=p1; family="Jac");
    if (family=="Leg", p1=0; p2=0; family="Jac");
    if (family=="Geg", p1=p1-1/2; p2=p1; family="Jac");
    if (family=="Jac", u=qe_jacobi(deg,p1,p2); shifted=1);
    if (family=="Lag", u=qe_laguerre(deg,p1));
    if (family=="Her", u=ab_hermite(deg,p1));
    if (family=="Ref", u=ab_refinable(npts,p1)) );
\\ ---
  if (is_AB(u), return(xw_ab(ab_modify(u))));
  if (is_QE(u), u=dw_qe(qe_modify(u)); if (shifted, u[1]=u[1]-1);
    return(xw_dw(u)),
    error("Programming error: by now data should be QE, but is ",Type(u)));
}

{DW(u,npts,p1,p2)=my(m,n,deg); 
  if (is_DW(u), return(u)); 
  if (is_AB(u), return(dw_xw(xw_ab(u))));
  if (is_SF(u), u=qe_sf(u));
  if (is_QE(u), return(dw_qe(u)));
  if (is_XW(u), return(dw_xw(u)));
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
}

MX(U,n)= { my(typ,serp=default(seriesprecision));
  if(n, default(seriesprecision,n));
  if(type(U)=="t_STR", typ=findpolynomial(U);
    if (typ=="Che", p1=-1/2; p2=p1; typ="Jac");
    if (typ=="Leg", return(log((1+t)/(1-t)))/t);
    if (typ=="Geg", p1=p1-1/2; p2=p1; typ="Jac");
    if (typ=="Jac", return(dw_qe(qe_jacobi(deg,p1,p2))));
    if (typ=="Lag", return(dw_qe(qe_laguerre(deg,p1))));
    if (typ=="Her", return(dw_xw(xw_ab(ab_hermite(deg,p1))))));
  if(is_SF(U), subst(U,z,1/t)+O(t^(Dimension(U)+1)));
}

{MU(U,R)=my(typ,mu); typ=Type(R);
  if(is_QE(U) & typ=="QE", return(mu_qe(U,R)));
  if(is_DW(U), U=XW(U));
  if(is_XW(U), return(mu_xw(U,R)))} 

{NU(U,n,family)=my(Z);
  if(is_NU(U),return(u));
  if (type(family)!="t_STR",
    if (type(n)=="t_STR", family=n; n=0));
  if(!n, n=Dimension(U));
print("n="n"  "U);
  if (is_DW(U), U=XW(U));
  if (is_XW(U), return(nu_xw(U));
  if (is_QE(U), return(nu_qe(U)));
  error("Not implemented: NU("Type(U)",n)"))}

{SF(u)=
  if(is_SF(u), return(u));
  if(is_DW(u), u=xw_dw(u),
  if(is_NU(u), u=qe_mu(u)));
  if(is_XW(u), return(sf_xw(u)));
  if(is_AB(u), return(sf_ab(u)));  
  if(!is_QE(u), u=QE(u));
  if(is_QE(u), return(sf_qe(u)))
}

{ Convert(U,V) =
  if (is_QE(V), QE(U),
  if (is_AB(V), AB(U),
  if (is_XW(V), XW(U),
  if (is_DW(V), DW(U),
  if (is_SF(V), SF(U),
  if (is_NU(V), NU(U)))))))
}
addhelp(Convert,"Convert(U,V): convert U to Type(V)")

\\{Shift(u,I1,I2)=my(typ); typ=Type(u);
\\  if (typ=="AB", return(ab_shift(u,I1,I2)));
\\  if (typ=="XW", return(xw_shift(u,I1,I2)));
\\  if (typ=="DW", return(dw_shift(u,I1,I2)));
\\  if (typ=="MU", error("Not implemented: Shift MU"));
\\  if (typ=="QE", return(qe_shift(u,I1,I2)))}
{Shift(u,d)=
  if (is_AB(u), for(k=1,#u[1], u[1][k]+=d));
  if (is_XW(u), for(k=1,matsize(u)[1], u[k,1]+=d));
  if (is_DW(u), u[1]+=d);
  if (is_QE(u), u=stqd(u,-d));
  if (is_SF(u), u=subst(u,z,z-d));
  u
}
addhelp(Shift,"Shift(u,d): shift origin of support interval by d")

{Scale(u,h)=my(c);
  if (is_AB(u), u[1]*=h; u[2]*=h^2; u[2][1]/=h);
  if (is_XW(u), u*=h);
  if (is_DW(u), u*=h);
  if (is_QE(u), u*=h);
  if (is_NU(u), c=1; for(k=2,#u, c*=h; u[k]*=c));
  if (is_SF(u), u=subst(u,x,x/h));
  u
}
addhelp(Scale,"Scale(u,h): scale length of support interval by h")

{Dimension(u)=
  if(is_NU(u), return(#u));
  if(is_XW(u), return(2*matsize(u)[1]));
  if(is_DW(u), return(2*matsize(u[2])[1]));
  if(is_AB(u), return(#u[1]+#u[2]));
  if(is_QE(u), return(#u[1]+#u[2]+1));
  if(is_SF(u), return(#numerator(u)+#denominator(u)-1));
}
addhelp(Dimension,"Dimension(u,{typ=Type(u)}): number of significant elements in u.")

{Degree(u,v)=
  if(!v,v=QE(Dimension(u)));
  if(!is_QE(u), u=QE(u));
  if(!is_QE(v), v=QE(v));
  exactness(u,v)}
addhelp(Degree,"Degree(u[,v]): degree of u, taking v as exact. Weight w(x)=1 over [0,2] by default.")

{Take(R,dim)=my(n=sign(dim)*floor(abs(dim)/2));
  if(dim==0, return(R));
  if (is_NU(R), return(take(R,dim)));
  if (is_QE(R), R[1]=take(R[1],n); R[2]=take(R[2],dim-n-sign(dim)),
    if (is_AB(R), R[1]=take(R[1],n); R[2]=take(R[2],dim-n),
      error("Can't Take from Type "Type(R))));
     R}

porthos_types=["QE","AB","XW","DW","NU","SF"]
{Type(gen)= my(name); 
  for(i=1,#porthos_types, name=porthos_types[i];
    if(eval(concat(["is_",name,"(gen)"])),return(name)) );
  type(gen) }

{ print1("   <PORTHOS "floor(PORTHOS_VERSION/100)"."
   floor((PORTHOS_VERSION%100)/10)""PORTHOS_VERSION%10">"); 
  read("porthos-help.gp"); }

default(breakloop,0);
