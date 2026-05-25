function roe_solver
global gamma ga gb gc gd ge gf gg gh gi gj % era USE PIgamma_PAR
global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx                     % era USE GEOM_PAR
global   time dt dtodx dtitest4 timemax timek     % era USE TIME_PAR
global  xsh10 xsh20 p1 rho1 u1 p2 rho2 u2 p3 rho3 u3 p4 rho4 u4 % era USE INIT_PAR
global  vsh1 vsh2 vsh3 vsh4 vsh5 vsh3l vsh3r vsh4l vsh4r % era USE INIT_PAR
global  xsh1 xsh2 xsh3 xsh4 xsh5 xsh3l xsh3r xsh4l xsh4r % era USE INIT_PAR
global p t u s rho a e amach ptot ttot flow flht h htot % era USE VARS
global w1 w2 w3 f1 f2 f3 phi1 phi2 phi3                    % era USE VARS
global pxeno uxeno hxeno ppxeno hhxeno ppt ut hht % era USE ENO

enuo2 = 0.5*dt/dx;

ncmm = ncm-1;
for n=2:ncmm
    nm=n;
    np=nm+1;
    
    pa  =  p(nm);
    pb  =  p(np);
    ua  =  u(nm);
    ub  =  u(np);
    ha  =  h(nm);
    hb  =  h(np);
    
    if n == 2
        p002   =  pa;
        u002   =  ua;
        h002   =  ha;
        rho002 =  p002/h002*ga;
        a002   =  sqrt(gamma*p002/rho002);
    end
    
    if n == ncmm
        pncm   =  pb;
        uncm   =  ub;
        hncm   =  hb;
        rhoncm =  pncm/hncm*ga;
        ancm   =  sqrt(gamma*pncm/rhoncm);
    end
    
    
    if iord ~= 1
        ppa = log(pa);
        ppb = log(pb);
        hha = log(ha);
        hhb = log(hb);
        ppa  =  ppa + .5* ppxeno(n);
        ppb  =  ppb - .5* ppxeno(n+1);
        ua   =  ua  + .5* uxeno(n);
        ub   =  ub  - .5* uxeno(n+1);
        hha  =  hha + .5* hhxeno(n);
        hhb  =  hhb - .5* hhxeno(n+1);
        
        ppa = ppa + enuo2*ppt(nm);
        ua  = ua  + enuo2*ut(nm);
        hha = hha + enuo2*hht(nm);
        
        if n == 2
            pp002  =  log(p(nm));
            u002   =  u(nm);
            hh002  =  log(h(nm));
            pp002  =  pp002 - .5* ppxeno(n);
            u002   =  u002  - .5* uxeno(n);
            hh002  =  hh002 - .5* hhxeno(n);
            pp002  =  pp002 + enuo2*ppt(n);
            u002   =  u002  + enuo2*ut(n);
            hh002  =  hh002 + enuo2*hht(n);
            p002   =  exp(pp002);
            h002   =  exp(hh002);
            rho002 =  p002/h002*ga;
            a002   =  sqrt(gamma*p002/rho002);
        end
        
        ppb = ppb + enuo2*ppt(np);
        ub  = ub  + enuo2*ut(np);
        hhb = hhb + enuo2*hht(np);
        
        if n == ncmm
            ppncm  =  log(p(np));
            uncm   =  u(np);
            hhncm  =  log(h(np));
            ppncm  =  ppncm    + .5* ppxeno(n+1);
            uncm   =  uncm     + .5* uxeno(n+1);
            hhncm  =  hhncm    + .5* hhxeno(n+1);
            ppncm  =  ppncm    + enuo2*ppt(n+1);
            uncm   =  uncm     + enuo2*ut(n+1);
            hhncm  =  hhncm    + enuo2*hht(n+1);
            pncm   =  exp(ppncm);
            hncm   =  exp(hhncm);
            rhoncm =  pncm/hncm*ga;
            ancm   =  sqrt(gamma*pncm/rhoncm);
        end
        
        pa = exp(ppa);
        pb = exp(ppb);
        ha = exp(hha);
        hb = exp(hhb);
    end
    

    rhoa = pa/ha*ga;
    rhob = pb/hb*ga;
    aa   = sqrt(gamma*pa/rhoa);
    ab   = sqrt(gamma*pb/rhob);
    
    %icalc=0;
    
  
    
    %{
    if pc <= 0.0 || pd <= 0.0 || hc<= 0.0 || hd <= 0.0
        icalc=1;
    end
    
    if icalc == 1
        %ppa = log(pa);
        %ppb = log(pb);
        %hha = log(ha);
        %hhb = log(hb);
        %r3a = ppa + gamma/aa*ua;
        %r2a = hha-ppa/ga;
        %r2b = hhb-ppb/ga;
        %r1b = ppb-gamma/ab*ub;
        uc = (r3a-r1b)/(gamma/aa+gamma/ab);
        ud = uc;
        ppc= r3a-gamma/aa*uc;
        ppd= ppc;
        hhc= ppc/ga+r2a;
        hhd= ppd/ga+r2b;
        pc = exp(ppc);
        pd = exp(ppd);
        hc = exp(hhc);
        hd = exp(hhd);
        fprintf('warning icalc=1 at k=%i and n=%i',k,n)
    end
    
    %}
    
    [f1a,f2a,f3a] = decod(pa,ua,ha);
    [f1b,f2b,f3b] = decod(pb,ub,hb);

  

    %inizio del metodo di ROE 
    ubar = (sqrt(rhoa)*ua + sqrt(rhob)*ub)/(sqrt(rhoa)+sqrt(rhob));
    hbar= (sqrt(rhoa)*ha + sqrt(rhob)*hb)/(sqrt(rhoa)+sqrt(rhob));
    abarquadro = (gamma-1)*(hbar - 0.5*ubar^2);
    lambdabar1 = ubar-sqrt(abarquadro);
    lambdabar2= ubar;
    lambdabar3= ubar+sqrt(abarquadro);

   

   
    
    %edivisorhobar = hbar-abarquadro/gamma;

%{
    A = [0 1 0; 
        -(3-gamma)*ubar.^2/2 (3-gamma)*ubar gamma-1; 
        (gamma-1)*ubar.^3-gamma*ubar*edivisorhobar gamma*edivisorhobar-1.5*(gamma-1)*ubar.^2 gamma*ubar];
    lambdabar1= ubar-sqrt(abarquadro);
    lambdabar2= ubar;
    lambdabar3=ubar+sqrt(abarquadro);
%}
    Ua = [rhoa;
      rhoa*ua;
      pa*gb+.5*rhoa*ua^2];

Ub = [rhob;
      rhob*ub;
      pb*gb+.5*rhob*ub^2];
    
deltaU = Ub-Ua;
    

    Rbar= [1 1 1; 
        lambdabar1 lambdabar2 lambdabar3; 
        hbar-sqrt(abarquadro)*ubar ubar.^2/2 hbar+sqrt(abarquadro)*ubar];
    
    ws = Rbar\deltaU;
%{
    if any(~isfinite(rho)) || any(~isfinite(u)) || any(~isfinite(p))
    error('NaN/Inf trovato: instabilità esplosa al passo %d', k);
end

if any(rho <= 0)
    warning('densità negativa al passo %d', k);
end

if any(p <= 0)
    warning('pressione negativa al passo %d', k);
end
%}
   %{
 rbar1 = Rbar(1,:)';
    rbar2 = Rbar(2,:)';
    rbar3 = Rbar(3,:)';
   %}
   

   %{ 
DNF1= f1b-f1a;
    DNF2= f2b-f2a;
    DNF3= f3b-f3a;
    DNF = [DNF1,DNF2,DNF3]';
    %}
    %product_matrix = Rbar*lambda_matrix;
    


   % if any(isnan(product_matrix(:)))
    %error('NaN in product_matrix');
%end


    %ws = product_matrix\DNF;
   % if rcond(product_matrix) < 1e-12
    %fprintf('ROE FALLBACK k=%d n=%d rcond=%e\n',k,n,rcond(product_matrix));

    % fallback semplice e stabile
    %ws = pinv(product_matrix) * DNF;
%else
    %ws = product_matrix \ DNF;
%end
%{
    DNFL = 0.5*(lambdabar1-abs(lambdabar1))*rbar1*ws(1)+0.5*(lambdabar2-abs(lambdabar2))*rbar2*ws(2)+0.5*(lambdabar3-abs(lambdabar3))*rbar3*ws(3);
    DNF1L=DNFL(1);
    DNF2L=DNFL(2);
    DNF3L =DNFL(3);
    DNFR = 0.5*(lambdabar1+abs(lambdabar1))*rbar1*ws(1)+0.5*(lambdabar2+abs(lambdabar2))*rbar2*ws(2)+0.5*(lambdabar3+abs(lambdabar3))*rbar3*ws(3);
    DNF1R=DNFR(1);
    DNF2R=DNFR(2);
    DNF3R=DNFR(3);

    phi1(n)=f1a+DNF1L;    %Fn+1/2
    phi2(n)=f2a+DNF2L;
    phi3(n)=f3a+DNF3L;
     %}
Lambda_abs = diag(abs([lambdabar1 lambdabar2 lambdabar3]));

Rabs = Rbar * Lambda_abs;
Fa = [f1a f2a f3a]';
Fb = [f1b f2b f3b]';
flux_diff = Rabs * ws; 
F = 0.5*( Fa + Fb )- 0.5*flux_diff;

phi1(n) = F(1);  %Fn+1/2
phi2(n) = F(2);
phi3(n) = F(3);




end
if(itest == 1)  %REFLECTING WALL B.C.
    r1dum  = p002-rho002*a002*u002;
    r2dum  = h002-p002/rho002;
    pin    = r1dum;
    uin    = 0.0;
    hin    = r2dum+pin/rho002;
    [phi1(1),phi2(1),phi3(1)] = decod(pin,uin,hin);
    r3dum  = pncm+rhoncm*ancm*uncm;
    r2dum  = hncm-pncm/rhoncm;
    pex    = r3dum;
    uex    = 0.0;
    hex    = r2dum+pex/rhoncm;
    [phi1(ncm),phi2(ncm),phi3(ncm)] = decod(pex,uex,hex);
end
if(itest >= 2)  % REFLECTING WALL B.C.
    r1dum  = p002-rho002*a002*u002;
    r2dum  = h002-p002/rho002;
    pin    = p002;
    uin    = u002;
    hin    = h002;
    [phi1(1),phi2(1),phi3(1)] = decod(pin,uin,hin);
    r3dum  = pncm+rhoncm*ancm*uncm;
    r2dum  = hncm-pncm/rhoncm;
    pex    = pncm;
    uex    = uncm;
    hex    = hncm;
    [phi1(ncm),phi2(ncm),phi3(ncm)] = decod(pex,uex,hex);
end
end