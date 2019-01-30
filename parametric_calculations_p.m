function [level, hata, Ids_app] = parametric_calculations_p ( temp, vg, vd, vs, vb, vth, vth_app, Wdrawn, Ldrawn, model_address, command, Ids_aim, Ids, I, I_aim, vdsat, vdsat_difference, ocean_address_vgs, ocean_address_low_vds, ocean_address_vds, result_address_vgs, result_address_vgs_low_vds, result_address_vds, ub1_n_max, step, level, baslangic, son, short_ch, type )
% detayli yorumlar icin parametric_calculations_n fonksiyonuna bakabilirsiniz, algoritma ayni mantikla calisir sadece denklemler PMOS icin degistirilmistir
    hata = 0;
%%% Dosya Okuma    
    file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);


for i=baslangic:son
    
     k=strfind(C{1}{i},'kt1=');
     if(k==1)
        kt1_n = textscan(C{1}{i},'kt1=%f');
        kt1_n = kt1_n{1};
        kt1_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'kt1l=');
     if(k==1)
        kt1l_n = textscan(C{1}{i},'kt1l=%f');
        kt1l_n = kt1l_n{1};
        kt1l_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'kt2=');
     if(k==1)
        kt2_n = textscan(C{1}{i},'kt2=%f');
        kt2_n = kt2_n{1};
        kt2_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'ute=');
     if(k==1)
        ute_n = textscan(C{1}{i},'ute=%f');
        ute_n = ute_n{1};
        ute_n_char = C{1}{i};        
     end 
     
     k=strfind(C{1}{i},'u0=');
     if(k==1)
        u0_n = textscan(C{1}{i},'u0=%f');
        u0_n = u0_n{1};
     end
     
     k=strfind(C{1}{i},'ua=');
     if(k==1)
        ua_n = textscan(C{1}{i},'ua=%f');
        ua_n = ua_n{1};
     end
     
     k=strfind(C{1}{i},'ua1=');
     if(k==1)
        ua1_n = textscan(C{1}{i},'ua1=%f');
        ua1_n = ua1_n{1};
        ua1_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'ub=');
     if(k==1)
        ub_n = textscan(C{1}{i},'ub=%f');
        ub_n = ub_n{1};     
     end
     
     k=strfind(C{1}{i},'ub1=');
     if(k==1)
        ub1_n = textscan(C{1}{i},'ub1=%f');
        ub1_n = ub1_n{1};
        ub1_n_char = C{1}{i};        
     end
     
     k=strfind(C{1}{i},'uc=');
     if(k==1)
        uc_n = textscan(C{1}{i},'uc=%f');
        uc_n = uc_n{1};
     end
     
     k=strfind(C{1}{i},'uc1=');
     if(k==1)
        uc1_n = textscan(C{1}{i},'uc1=%f');
        uc1_n = uc1_n{1};
        uc1_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'at=');
     if(k==1)
        at_n = textscan(C{1}{i},'at=%f');
        at_n = at_n{1};
        at_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'prt=');
     if(k==1)
        prt_n = textscan(C{1}{i},'prt=%f');
        prt_n = prt_n{1};
        prt_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'delta=');
     if(k==1)
        delta_n = textscan(C{1}{i},'delta=%f');
        delta_n = delta_n{1};
        delta_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'vth0=');
     if(k==1)
        vth0_n = textscan(C{1}{i},'vth0=%f');
        vth0_n = vth0_n{1};
     end
         
     k=strfind(C{1}{i},'rdsw=');
     if(k==1)
        rdsw_n = textscan(C{1}{i},'rdsw=%f');
        rdsw_n = rdsw_n{1};
     end
     
     k=strfind(C{1}{i},'vsat=');
     if(k==1)
        vsat_n = textscan(C{1}{i},'vsat=%f');
        vsat_n = vsat_n{1};
     end
     
     k=strfind(C{1}{i},'k1=');
     if(k==1)
        k1_n = textscan(C{1}{i},'k1=%f');
        k1_n = k1_n{1};
     end
     
     k=strfind(C{1}{i},'k2=');
     if(k==1)
        k2_n = textscan(C{1}{i},'k2=%f');
        k2_n = k2_n{1};
     end
     
     k=strfind(C{1}{i},'a0=');
     if(k==1)
        a0 = textscan(C{1}{i},'a0=%f');
        a0 = a0{1};
     end
     
     k=strfind(C{1}{i},'a1=');
     if(k==1)
        a1 = textscan(C{1}{i},'a1=%f');
        a1 = a1{1};
     end
     
     k=strfind(C{1}{i},'a2=');
     if(k==1)
        a2 = textscan(C{1}{i},'a2=%f');
        a2 = a2{1};
     end
     
     k=strfind(C{1}{i},'ags=');
     if(k==1)
        ags = textscan(C{1}{i},'ags=%f');
        ags = ags{1};
     end
     
     k=strfind(C{1}{i},'b0=');
     if(k==1)
        b0 = textscan(C{1}{i},'b0=%f');
        b0 = b0{1};
     end
     
     k=strfind(C{1}{i},'b1=');
     if(k==1)
        b1 = textscan(C{1}{i},'b1=%f');
        b1 = b1{1};
     end
     
     k=strfind(C{1}{i},'cdsc=');
     if(k==1)
        cdsc = textscan(C{1}{i},'cdsc=%f');
        cdsc = cdsc{1};
     end
     
     k=strfind(C{1}{i},'cdscb=');
     if(k==1)
        cdscb = textscan(C{1}{i},'cdscb=%f');
        cdscb = cdscb{1};
     end
     
     k=strfind(C{1}{i},'cdscd=');
     if(k==1)
        cdscd = textscan(C{1}{i},'cdscd=%f');
        cdscd = cdscd{1};
     end
     
     k=strfind(C{1}{i},'cit=');
     if(k==1)
        cit = textscan(C{1}{i},'cit=%f');
        cit = cit{1};
     end
     
     k=strfind(C{1}{i},'drout=');
     if(k==1)
        drout = textscan(C{1}{i},'drout=%f');
        drout = drout{1};
     end
     
     k=strfind(C{1}{i},'dsub=');
     if(k==1)
        dsub = textscan(C{1}{i},'dsub=%f');
        dsub = dsub{1};
     end
     
     k=strfind(C{1}{i},'dvt0=');
     if(k==1)
        dvt0 = textscan(C{1}{i},'dvt0=%f');
        dvt0 = dvt0{1};
     end
     
     k=strfind(C{1}{i},'dvt0w=');
     if(k==1)
        dvt0w = textscan(C{1}{i},'dvt0w=%f');
        dvt0w = dvt0w{1};
     end
     
     k=strfind(C{1}{i},'dvt1=');
     if(k==1)
        dvt1 = textscan(C{1}{i},'dvt1=%f');
        dvt1 = dvt1{1};
     end
     
     k=strfind(C{1}{i},'dvt1w=');
     if(k==1)
        dvt1w = textscan(C{1}{i},'dvt1w=%f');
        dvt1w = dvt1w{1};
     end
     
     k=strfind(C{1}{i},'dvt2=');
     if(k==1)
        dvt2 = textscan(C{1}{i},'dvt2=%f');
        dvt2 = dvt2{1};
     end
     
     k=strfind(C{1}{i},'dvt2w=');
     if(k==1)
        dvt2w = textscan(C{1}{i},'dvt2w=%f');
        dvt2w = dvt2w{1};
     end
     
     k=strfind(C{1}{i},'dwb=');
     if(k==1)
        dwb = textscan(C{1}{i},'dwb=%f');
        dwb = dwb{1};
     end
     
     k=strfind(C{1}{i},'dwg=');
     if(k==1)
        dwg = textscan(C{1}{i},'dwg=%f');
        dwg = dwg{1};
     end
     
     k=strfind(C{1}{i},'eta0=');
     if(k==1)
        etao = textscan(C{1}{i},'eta0=%f');
        etao = etao{1};
     end
     
     k=strfind(C{1}{i},'etab=');
     if(k==1)
        etab = textscan(C{1}{i},'etab=%f');
        etab = etab{1};
     end
     
     k=strfind(C{1}{i},'k3=');
     if(k==1)
        k3 = textscan(C{1}{i},'k3=%f');
        k3 = k3{1};
     end
     
     k=strfind(C{1}{i},'k3b=');
     if(k==1)
        k3b = textscan(C{1}{i},'k3b=%f');
        k3b = k3b{1};
     end
     
     k=strfind(C{1}{i},'keta=');
     if(k==1)
        keta = textscan(C{1}{i},'keta=%f');
        keta = keta{1};
     end
          
     k=strfind(C{1}{i},'lint=');
     if(k==1)
        lint = textscan(C{1}{i},'lint=%f');
        lint = lint{1};
     end
     
     k=strfind(C{1}{i},'ll=');
     if(k==1)
        ll = textscan(C{1}{i},'ll=%f');
        ll = ll{1};
     end
     
     k=strfind(C{1}{i},'lln=');
     if(k==1)
        lln = textscan(C{1}{i},'lln=%f');
        lln = lln{1};
     end
     
     k=strfind(C{1}{i},'lw=');
     if(k==1)
        lw = textscan(C{1}{i},'lw=%f');
        lw = lw{1};
     end
     
     k=strfind(C{1}{i},'lwn=');
     if(k==1)
        lwn = textscan(C{1}{i},'lwn=%f');
        lwn = lwn{1};
     end
     
     k=strfind(C{1}{i},'lwl=');
     if(k==1)
        lwl = textscan(C{1}{i},'lwl=%f');
        lwl = lwl{1};
     end
     
    
     k=strfind(C{1}{i},'nlx=');
     if(k==1)
        nlx = textscan(C{1}{i},'nlx=%f');
        nlx = nlx{1};
     end
     
     k=strfind(C{1}{i},'w0=');
     if(k==1)
        w0 = textscan(C{1}{i},'w0=%f');
        w0 = w0{1};
     end
          
     k=strfind(C{1}{i},'wint=');
     if(k==1)
        wint = textscan(C{1}{i},'wint=%f');
        wint = wint{1};
     end
     
     k=strfind(C{1}{i},'wl=');
     if(k==1)
        wl = textscan(C{1}{i},'wl=%f');
        wl = wl{1};
     end
     
     k=strfind(C{1}{i},'wln=');
     if(k==1)
        wln = textscan(C{1}{i},'wln=%f');
        wln = wln{1};
     end
     
     k=strfind(C{1}{i},'ww=');
     if(k==1)
        ww = textscan(C{1}{i},'ww=%f');
        ww = ww{1};
     end
     
     k=strfind(C{1}{i},'wwn=');
     if(k==1)
        wwn = textscan(C{1}{i},'wwn=%f');
        wwn = wwn{1};
     end
     
     k=strfind(C{1}{i},'wwl=');
     if(k==1)
        wwl = textscan(C{1}{i},'wwl=%f');
        wwl = wwl{1};
     end
          
     k=strfind(C{1}{i},'mobmod=');
     if(k==1)
        mobmod = textscan(C{1}{i},'mobmod=%f');
        mobmod = mobmod{1};
     end
          
     k=strfind(C{1}{i},'nch=');
     if(k==1)
        nch = textscan(C{1}{i},'nch=%f');
        nch = nch{1};
     end
              
     k=strfind(C{1}{i},'nfactor=');
     if(k==1)
        nfactor = textscan(C{1}{i},'nfactor=%f');
        nfactor = nfactor{1};
     end
          
     k=strfind(C{1}{i},'prwg=');
     if(k==1)
        prwg = textscan(C{1}{i},'prwg=%f');
        prwg = prwg{1};
     end
          
     k=strfind(C{1}{i},'prwb=');
     if(k==1)
        prwb = textscan(C{1}{i},'prwb=%f');
        prwb = prwb{1};
     end
          
     k=strfind(C{1}{i},'pscbe1=');
     if(k==1)
        pscbe1 = textscan(C{1}{i},'pscbe1=%f');
        pscbe1 = pscbe1{1};
     end
               
     k=strfind(C{1}{i},'pscbe2=');
     if(k==1)
        pscbe2 = textscan(C{1}{i},'pscbe2=%f');
        pscbe2 = pscbe2{1};
     end
          
     k=strfind(C{1}{i},'pdiblc1=');
     if(k==1)
        pdiblc1 = textscan(C{1}{i},'pdiblc1=%f');
        pdiblc1 = pdiblc1{1};
     end
          
     k=strfind(C{1}{i},'pdiblc2=');
     if(k==1)
        pdiblc2 = textscan(C{1}{i},'pdiblc2=%f');
        pdiblc2 = pdiblc2{1};
     end
          
     k=strfind(C{1}{i},'pdiblcb=');
     if(k==1)
        pdiblcb = textscan(C{1}{i},'pdiblcb=%f');
        pdiblcb = pdiblcb{1};
     end
          
     k=strfind(C{1}{i},'pvag=');
     if(k==1)
        pvag = textscan(C{1}{i},'pvag=%f');
        pvag = pvag{1};
     end
          
     k=strfind(C{1}{i},'pclm=');
     if(k==1)
        pclm = textscan(C{1}{i},'pclm=%f');
        pclm = pclm{1};
     end
          
     k=strfind(C{1}{i},'voff=');
     if(k==1)
        voff = textscan(C{1}{i},'voff=%f');
        voff = voff{1};
     end
          
     k=strfind(C{1}{i},'wr=');
     if(k==1)
        wr = textscan(C{1}{i},'wr=%f');
        wr = wr{1};
     end
          
     k=strfind(C{1}{i},'xj=');
     if(k==1)
        xj = textscan(C{1}{i},'xj=%f');
        xj = xj{1};
     end
     
     k=strfind(C{1}{i},'tnom=');
     if(k==1)
        tnom = textscan(C{1}{i},'tnom=%f');
        tnom = tnom{1};
     end
     
     k=strfind(C{1}{i},'tox=');
     if(k==1)
        tox = textscan(C{1}{i},'tox=%f');
        tox = tox{1};
     end
     
     k=strfind(C{1}{i},'toxm=');
     if(k==1)
        toxm = textscan(C{1}{i},'toxm=%f');
        toxm = toxm{1};
     end
end


    if exist('wr','var') == 0
        wr = 1;
    elseif isempty(wr) == 1
        wr = 1;
    end
    
    if exist('delta','var') == 0
        delta_n = 0;
    end

fclose(file);
%%%

l = 5*10^-6;
w = 5*10^-6;

% parametrelerin dosyadan okunan ilk değerleri
kt1_n_def   = kt1_n;
ua1_n_def   = ua1_n;
ute_n_def   = ute_n;
prt_n_def   = prt_n;
at_n_def    = at_n ;
delta_n_def = delta_n;
ub1_n_def = ub1_n;   

%if step==0
% -------------------------- Denklemler -----------------------------------
vth = abs(vth);
   nch = nch*10^6;
% i=1;   
% for vd=1.8:-0.01:0
vgs = -(vg - vs);
vds = -(vd - vs);
vbs = (vb - vs);
c_temp = (temp+273.15)/(tnom+273.15);

k1=1.38064852*10^-23;
q=1.602176565*10^-19;
Esi=8.85*10^-12*11.7055;
Eox=8.85*10^-12*3.9018;
%Eg=1.160-((7.02*10^-4*(tnom+273.15)^2)/(tnom+273.15+1108));
%ni=5.29*10^19*(c_temp)^2.54*exp(-6726/(273+temp))*10^6;
%ni = 1.45*(10^10)*(tnom/300.15)*sqrt(tnom/300.15)*exp(21.5565981-(Eg)/(2*k1*(tnom+273.15)));
ni = 3.1*10^16*(temp+273.15)^1.5*exp(-7000/(temp+273.15));
Cox=Eox/tox;

Vt = k1*(temp+273.15)/q;
phi_s = 2*k1*(tnom+273.15)/q*log(nch/ni);
Vsat = vsat_n - at_n*(c_temp-1);

dL = -lint + ll/(l^lln) + lw/(w^lwn) + lwl/(l^lln*w^lwn);
Leff = Ldrawn - 2*dL;

Vx = 0.9*(phi_s-(k1_n^2/(4*k2_n^2)));
if Vx>-3
    Vbc = -3;
elseif Vx<-30
    Vbc = -30;
else
    Vbc = Vx;
end
Vbseff = Vbc + 0.5*(vbs - Vbc- 0.001 + sqrt((vbs-Vbc-0.001)^2 - 0.004*Vbc));

Xdep = sqrt((2*Esi*(phi_s-Vbseff))/(q*nch));
Xdepo = sqrt(2*Esi*phi_s/(q*nch));

Lt = sqrt(Esi*Xdep/Cox)*(1+dvt2*Vbseff);
lto = sqrt(Esi*Xdepo/Cox);

Cd = Esi/Xdep;
Cdepo = Esi/Xdepo;

n = 1+nfactor*Cd/Cox+(((cdsc+cdscd*vds+cdscb*Vbseff)*(exp(-dvt1*Leff/(2*Lt))+2*exp(-dvt1*Leff/Lt)))/Cox)+cit/Cox;
Vgsteff = (2*n*Vt*log(1+exp((vgs-vth)/(2*n*Vt))))/(1+2*n*Cox*sqrt((2*phi_s)/(q*Esi*nch))*exp(-(vgs-vth-2*voff)/(2*n*Vt)));
%Vgsteff = (2*n*Vt*log(1+exp((vgs-vth)/(2*n*Vt))))/(1+2*n*(Cox/Cdepo)*exp(-(vgs-vth-2*voff)/(2*n*Vt)));

dW1 = -wint + wl/(l^wln) + ww/(w^wwn) + wwl/(l^wln*w^wwn);
Weff1 = Wdrawn -2*dW1;
dW = dW1 + dwg*Vgsteff + dwb*(sqrt(phi_s-Vbseff)-sqrt(phi_s));
if Vbseff>0
    dW = dW1 + dwg*Vgsteff + dwb*((sqrt(phi_s)^3/(phi_s+Vbseff/2))-sqrt(phi_s));
end
Weff = Wdrawn - 2*dW;

if Weff<2*10^-8
   Weff = 2*10^-8*((4e-08-Weff)/(6e-08-2*Weff)); 
end
if Weff1<2*10^-8
   Weff1 = 2*10^-8*((4e-08-Weff1)/(6e-08-2*Weff1)); 
end

Rdsw = rdsw_n+prt_n*(c_temp-1);
Rds = (Rdsw*(1+prwg*Vgsteff+prwb*(sqrt(phi_s-Vbseff)-sqrt(phi_s))))/((10^6*Weff)^wr);

if mobmod==1
Tmp = (((ua_n+ua1_n*(c_temp-1))+(uc_n+uc1_n*(c_temp-1))*Vbseff)*((Vgsteff+2*vth)/tox)+(ub_n+ub1_n*(c_temp-1))*(((Vgsteff+2*vth)/tox)^2));
elseif mobmod==2
Tmp = (((ua_n+ua1_n*(c_temp-1))+(uc_n+uc1_n*(c_temp-1))*Vbseff)*((Vgsteff)/tox)+(ub_n+ub1_n*(c_temp-1))*(((Vgsteff)/tox)^2));    
elseif mobmod==3
Tmp = (((ua_n+ua1_n*(c_temp-1))*((Vgsteff+2*vth)/tox)+(ub_n+ub1_n*(c_temp-1))*(((Vgsteff+2*vth)/tox)^2))*(1+(uc_n+uc1_n*(c_temp-1))*Vbseff));    
end

if Tmp >= -0.8
    denominator = 1+Tmp;
elseif Tmp <= -0.8
    denominator = (0.6+Tmp)/(7+10*Tmp);
end

ute_app = ute_n_def;

mu = (u0_n*((c_temp)^ute_app))*10^-4/denominator;

Esat = 2*Vsat/mu;

Abulk_temp = (1+(k1_n/(2*sqrt(phi_s-Vbseff)))*(((a0*Leff)/(Leff+2*sqrt(xj*Xdep)))*(1-ags*Vgsteff*(Leff/(Leff+2*sqrt(xj*Xdep)))^2)+b0/(Weff+b1)));
if Abulk_temp >= 0.1
Abulk = Abulk_temp*(1/(1+keta*Vbseff));
elseif Abulk_temp <= 0.1
Abulk =((0.2-Abulk_temp)/(3-20*Abulk_temp))*(1/(1+keta*Vbseff));
end


lambda = a1*Vgsteff + a2;
a = Abulk^2*Weff*Vsat*Cox*Rds + (1/lambda-1)*Abulk;
b = -((Vgsteff+2*Vt)*(2/lambda-1) + Abulk*Esat*Leff + 3*Abulk*(Vgsteff+2*Vt)*Weff*Vsat*Cox*Rds);
c = (Vgsteff+2*Vt)*Esat*Leff + 2*((Vgsteff+2*Vt)^2)*Weff*Vsat*Cox*Rds;

Vdsat = (-b-sqrt(b^2-4*a*c))/(2*a);
Vdsat_app = Vdsat;
Vdseff = Vdsat - 0.5*(Vdsat-vds-delta_n+sqrt((Vdsat-vds-delta_n)^2+4*delta_n*Vdsat));
    
  Litl = sqrt((Esi*tox*xj)/Eox);

Vasat = (Esat*Leff + Vdsat + 2*Rds*Vsat*Cox*Weff*Vgsteff*(1-(Abulk*Vdsat)/(2*(Vgsteff+2*Vt))))/(2/lambda - 1 + Rds*Vsat*Cox*Weff*Abulk);

Vascbe = ((Leff/pscbe2)*exp((pscbe1*Litl)/(vds-Vdseff)));

theta_rout = pdiblc1*(exp(-drout*Leff/(2*lto))+2*exp(-drout*Leff/lto)) + pdiblc2;
Vadiblc = ((Vgsteff+2*Vt)/(theta_rout*(1+pdiblcb*Vbseff)))*(1-(Abulk*Vdsat)/(Abulk*Vdsat + Vgsteff + 2*Vt));
Vaclm = (Abulk*Esat*Leff+Vgsteff)*(vds-Vdseff)/(pclm*Abulk*Esat*Litl);
Va = Vasat + (1+(pvag*Vgsteff)/(Esat*Leff))*((Vaclm*Vadiblc)/(Vaclm+Vadiblc));
Idso = (Weff*mu*Cox*Vgsteff*(1-(Abulk*Vdseff)/(2*(Vgsteff+2*Vt)))*Vdseff)/(Leff*(1+Vdseff/(Esat*Leff)));

% FITTING FUNCTIONS
%[ratio, distance ] = pmos_width_length_shifting( vg, vd, -vbs, Wdrawn, Ldrawn, temp );

Ids_app = ((Idso/(1+Rds*Idso/Vdseff))*(1+(vds-Vdseff)/Va)*(1+(vds-Vdseff)/Vascbe));

%Ids_app(i,1) = ((Idso/(1+Rds*Idso/Vdseff))*(1+(vds-Vdseff)/Va)*(1+(vds-Vdseff)/Vascbe));
% i = i+1;
% end



% -------------------------------------------------------------------------

% ----------------------- Parametre Çıkarma--------------------------------

%%%%%%% KT1 %%%%%%%
if vbs == 0 && level == 1 
    
%27 derece default analiz
% const2 = abs(vth_def)-abs(vth0_n);
% const = ((vth-const2)-abs(vth0_n))/(c_temp-1)-kt1_n_def;
% kt1_n_app1 = ((vth_app-const2)-abs(vth0_n))/(c_temp-1)-const

%default analize gerek olmadan

const2 = vth-vth0_n-(kt1_n + kt1l_n/Leff + kt2_n*vbs)*(c_temp-1);
kt1_n_app2 = ((vth_app-const2)-vth0_n)/(c_temp-1)-kt1l_n/Leff-kt2_n*vbs

    kt1_n = kt1_n_app2;

%denklem ile
% Vth0ox = vth0_n + (kt1_n+kt1l_n/Leff+kt2_n*Vbseff)*(c_temp-1) - k1*sqrt(phi_s);
% K1ox = k1_n*tox/toxm;
% K2ox = k2_n*tox/toxm;
% lt = sqrt(Esi*Xdep/Cox)*(1+dvt2*Vbseff);
% ltw = sqrt(Esi*Xdep/Cox)*(1+dvt2w*Vbseff);
% nds = 10^20;
% Vbi = Vt * log(nch*nds/ni^2);
% 
% vth_equ = Vth0ox + K1ox*sqrt(phi_s-Vbseff)-K2ox*Vbseff + K1ox*(sqrt(1+nlx/Leff)-1)*sqrt(phi_s)...
%  + (k3+k3b*Vbseff)*(tox/(Weff+w0))*phi_s - dvt0w*(exp(-dvt1w*Weff*Leff/(2*ltw))...
%  +2*exp(-dvt1w*Weff*Leff/ltw))*(Vbi-phi_s) - dvt0*(exp(-dvt1*Leff/(2*lt))+2*exp(-dvt1*Leff/lt))*(Vbi-phi_s)...
%  -(exp(-dsub*Leff/(2*lto))+2*exp(-dsub*Leff/lto))*(etao + etab*Vbseff)*0;
  
%---- KT1 Dosyaya yazdırma
    kt1_n_char_new = sprintf('kt1=%d', kt1_n);
    kt1_n_char_old = kt1_n_char;
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması   
    execute = sprintf(command,kt1_n_char,kt1_n_char_new,model_address);
    system(execute);    
    % Parametrelerin dosyadan tekrar okunması (hata almamak için)
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file,formatSpec);
        for i=baslangic:son
    
     k=strfind(C{1}{i},'kt1=');
     if(k==1)
        kt1_n = textscan(C{1}{i},'kt1=%f');
        kt1_n = kt1_n{1};
        kt1_n_char = C{1}{i};
     end
        end  
        fclose(file);
%----
end

%%%%%%% KT1L %%%%%%%
if vbs == 0 && level == 1 && short_ch == 1 
const2 = vth-vth0_n-(kt1_n + kt1l_n/Leff + kt2_n*vbs)*(c_temp-1);
kt1l_app = (((vth_app-const2)-vth0_n)/(c_temp-1)-kt1_n-kt2_n*vbs)*Leff;

if isempty(kt1l_app)
        hata = 1;
       return
end
kt1l_n = kt1l_app;


%---- KT1L Dosyaya yazdırma
    kt1l_n_char_new = sprintf('kt1l=%d', kt1l_n);
    kt1l_n_char_old = kt1l_n_char;
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması   
    execute = sprintf(command,kt1l_n_char,kt1l_n_char_new,model_address);
    system(execute);    
    % Parametrelerin dosyadan tekrar okunması (hata almamak için)
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file,formatSpec);
        for i=baslangic:son
     k=strfind(C{1}{i},'kt1l=');
     if(k==1)
        kt1l_n = textscan(C{1}{i},'kt1l=%f');
        kt1l_n = kt1l_n{1};
        kt1l_n_char = C{1}{i};
     end
        end  
        fclose(file);
        
                    
        %%%%%% UB1 in KT1L %%%%%%          
    while max(abs(diff(I)))>mean(abs(diff(I((size(I,1)-20):size(I,1)))))*2
        if abs(ub1_n)<abs(ub1_n_max) 
        min_diff = min(diff(I,2));
        ub1_n = ub1_n+(ub1_n_max-ub1_n_def)/step; 

        % Değişkenlerin tekrardan karaktere dönüştürülmesi
        ub1_n_char_new = sprintf('ub1=%d', ub1_n);
        ub1_n_char_old = ub1_n_char;

        % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması  
        execute = sprintf(command,ub1_n_char,ub1_n_char_new,model_address);
        system(execute);

    %% Parametrelerin dosyadan tekrar okunması (hata almamak için)

    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file,formatSpec);

        for i=baslangic:son

             k=strfind(C{1}{i},'ub1=');
             if(k==1)
                ub1_n = textscan(C{1}{i},'ub1=%f');
                ub1_n = ub1_n{1};
                ub1_n_char = C{1}{i};
             end
        end  
    fclose(file);
    ub1_n

    %%

        %Yeni parametreler ile analiz
        ocean_command = sprintf('ocean -restore %s',ocean_address_vgs); 
        unix(ocean_command)

        result = importdata(result_address_vgs);
        new_data = result.data;
        V = new_data(:,1);
        I = new_data(:,2);

        %%

            if min(diff(I,2))<min_diff
               % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
               execute = sprintf(command,ub1_n_char_new,ub1_n_char_old,model_address);
               system(execute);
               break            
            end
           if abs(ub1_n)>=abs(ub1_n_max)
                break
           end     
        end

    end
end
%----

%%%%%% AT %%%%%%%%
% if vbs == 0 && level == 1
% vdsat_diff = vdsat-Vdsat_app;
% 
% syms at_n
% Vsat = vsat_n - at_n*(c_temp-1);
% 
% mu = (u0_n*((c_temp)^ute_app))*10^-4/denominator;
% Esat = 2*Vsat/mu;
% Abulk_temp = (1+(k1_n/(2*sqrt(phi_s-Vbseff)))*(((a0*Leff)/(Leff+2*sqrt(xj*Xdep)))*(1-ags*Vgsteff*(Leff/(Leff+2*sqrt(xj*Xdep)))^2)+b0/(Weff1+b1)));
% if Abulk_temp >= 0.1
% Abulk = Abulk_temp*(1/(1+keta*Vbseff));
% elseif Abulk_temp <= 0.1
% Abulk =((0.2-Abulk_temp)/(3-20*Abulk_temp))*(1/(1+keta*Vbseff));
% end
% 
% lambda = a1*Vgsteff + a2;
% a = Abulk^2*Weff*Vsat*Cox*Rds + (1/lambda-1)*Abulk;
% b = -((Vgsteff+2*Vt)*(2/lambda-1) + Abulk*Esat*Leff + 3*Abulk*(Vgsteff+2*Vt)*Weff*Vsat*Cox*Rds);
% c = (Vgsteff+2*Vt)*Esat*Leff + 2*((Vgsteff+2*Vt)^2)*Weff*Vsat*Cox*Rds;
% 
% Vdsat = vdsat_diff+(-b-sqrt(b^2-4*a*c))/(2*a);
% 
% eqn = abs(vdsat) + vdsat_difference - Vdsat;
% at_app = double(solve(eqn,at_n));
%     
%     at_n = at_app
% end
% %%%%%%%%%%%%%%%%%%%
% %%---- AT Dosyaya yazdırma
%     at_n_char_new = sprintf('at=%d', at_n);
%     at_n_char_old = at_n_char;
%     % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması 
%     execute = sprintf(command,at_n_char,at_n_char_new,model_address);
%     system(execute);   
%     %% Parametrelerin dosyadan tekrar okunması (hata almamak için)
%     file = fopen(model_address);
%     formatSpec = '%s';
%     C = textscan(file,formatSpec);
%         for i=baslangic:son
%     
%      k=strfind(C{1}{i},'at=');
%      if(k==1)
%         at_n = textscan(C{1}{i},'at=%f');
%         at_n = at_n{1};
%         at_n_char = C{1}{i};
%      end
%         end  
%         fclose(file);
% %%----
% [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
% [ V, I ] = simulation( ocean_address_vgs, result_address_vgs );
% [ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds );
% vth = importdata('./sonuclar/vth.dat');
% vdsat = importdata('./sonuclar/vdsat.dat');


%%%%%%% UTE %%%%%%%
if vbs == 0 && level == 1
    
[Ids_boyut,~] = size(Ids);
%ratio=Ids(Ids_boyut,1)/Ids_app;
%ratio_vdsat = vdsat/Vdsat_app;

fark=Ids(Ids_boyut,1)-Ids_app;
%fark_vdsat = vdsat-Vdsat_app;
Vsat = vsat_n - at_n*(c_temp-1);
    syms ute_app
mu = (u0_n*((c_temp)^ute_app))*10^-4/denominator;

Esat = 2*Vsat/mu;

Abulk_temp = (1+(k1_n/(2*sqrt(phi_s-Vbseff)))*(((a0*Leff)/(Leff+2*sqrt(xj*Xdep)))*(1-ags*Vgsteff*(Leff/(Leff+2*sqrt(xj*Xdep)))^2)+b0/(Weff+b1)));
if Abulk_temp >= 0.1
Abulk = Abulk_temp*(1/(1+keta*Vbseff));
elseif Abulk_temp <= 0.1
Abulk =((0.2-Abulk_temp)/(3-20*Abulk_temp))*(1/(1+keta*Vbseff));
end


lambda = a1*Vgsteff + a2;
a = Abulk^2*Weff*Vsat*Cox*Rds + (1/lambda-1)*Abulk;
b = -((Vgsteff+2*Vt)*(2/lambda-1) + Abulk*Esat*Leff + 3*Abulk*(Vgsteff+2*Vt)*Weff*Vsat*Cox*Rds);
c = (Vgsteff+2*Vt)*Esat*Leff + 2*((Vgsteff+2*Vt)^2)*Weff*Vsat*Cox*Rds;

Vdsat = (-b-sqrt(b^2-4*a*c))/(2*a);
%Vdsat = 0.863;
Vdseff = Vdsat - 0.5*(Vdsat-vds-delta_n+sqrt((Vdsat-vds-delta_n)^2+4*delta_n*Vdsat));
    
Litl = sqrt((Esi*tox*xj)/Eox);

Vasat = (Esat*Leff + Vdsat + 2*Rds*Vsat*Cox*Weff*Vgsteff*(1-(Abulk*Vdsat)/(2*(Vgsteff+2*Vt))))/(2/lambda - 1 + Rds*Vsat*Cox*Weff*Abulk);

Vascbe = ((Leff/pscbe2)*exp((pscbe1*Litl)/(vds-Vdseff)));

theta_rout = pdiblc1*(exp(-drout*Leff/(2*lto))+2*exp(-drout*Leff/lto)) + pdiblc2;
Vadiblc = ((Vgsteff+2*Vt)/(theta_rout*(1+pdiblcb*Vbseff)))*(1-(Abulk*Vdsat)/(Abulk*Vdsat + Vgsteff + 2*Vt));
Vaclm = (Abulk*Esat*Leff+Vgsteff)*(vds-Vdseff)/(pclm*Abulk*Esat*Litl);
Va = Vasat + (1+(pvag*Vgsteff)/(Esat*Leff))*((Vaclm*Vadiblc)/(Vaclm+Vadiblc));
Idso = (Weff*mu*Cox*Vgsteff*(1-(Abulk*Vdseff)/(2*(Vgsteff+2*Vt)))*Vdseff)/(Leff*(1+Vdseff/(Esat*Leff)));

Ids_app2 = fark+((Idso/(1+Rds*Idso/Vdseff))*(1+(vds-Vdseff)/Va)*(1+(vds-Vdseff)/Vascbe));

eqn = Ids_app2 - Ids_aim(Ids_boyut,1);
ute_app = double(solve(eqn,ute_app));

    if isempty(ute_app)
        syms ute_app
        ute_app = double(solve(simplify(eqn),ute_app));
        if isempty(ute_app)
            [ ute_n, Ids, I ] =ute_iteration_algorithm(command, model_address,step, I, I_aim, Ids, Ids_aim, ocean_address_vds, ocean_address_vgs, result_address_vds, result_address_vgs, baslangic, son, type );
            hata = 3;
            return            
        end
    end
    ute_n = ute_app

    
        %%%%%%%%% UB1 %%%%%%%%%%
    while max(abs(diff(I)))>mean(abs(diff(I((size(I,1)-20):size(I,1)))))*2
        if abs(ub1_n)<abs(ub1_n_max) 
        min_diff = min(diff(I,2));
        ub1_n = ub1_n+(ub1_n_max-ub1_n_def)/step; 
   
    % Değişkenlerin tekrardan karaktere dönüştürülmesi
    ub1_n_char_new = sprintf('ub1=%d', ub1_n);
    ub1_n_char_old = ub1_n_char;
    
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması  
    execute = sprintf(command,ub1_n_char,ub1_n_char_new,model_address);
    system(execute);
    
%% Parametrelerin dosyadan tekrar okunması (hata almamak için)

file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);

for i=baslangic:son
    
     k=strfind(C{1}{i},'ub1=');
     if(k==1)
        ub1_n = textscan(C{1}{i},'ub1=%f');
        ub1_n = ub1_n{1};
        ub1_n_char = C{1}{i};
     end
end  
fclose(file);
ub1_n

%%

    %Yeni parametreler ile analiz
    ocean_command = sprintf('ocean -restore %s',ocean_address_vgs); 
    unix(ocean_command)
       
    result = importdata(result_address_vgs);
    new_data = result.data;
    V = new_data(:,1);
    I = new_data(:,2);

    %%
   
        if min(diff(I,2))<min_diff
           % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
           execute = sprintf(command,ub1_n_char_new,ub1_n_char_old,model_address);
           system(execute);
           break            
        end
       if abs(ub1_n)>=abs(ub1_n_max)
            break
       end     
    end

    
end    
        %%%%%%%%%%%%%%%%%%%%%%%%
    
end


%%%%%%% KT2 %%%%%%%
if vbs < 0 && level == 2
const2 = vth-vth0_n-(kt1_n + kt1l_n/Leff + kt2_n*vbs)*(c_temp-1);
kt2_app = (((vth_app-const2)-vth0_n)/(c_temp-1)-kt1l_n/Leff-kt1_n)/Vbseff

kt2_n = kt2_app

%---- KT2 Dosyaya yazdırma
    kt2_n_char_new = sprintf('kt2=%d', kt2_n);
    kt2_n_char_old = kt2_n_char;
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması   
    execute = sprintf(command,kt2_n_char,kt2_n_char_new,model_address);
    system(execute);    
    % Parametrelerin dosyadan tekrar okunması (hata almamak için)
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file,formatSpec);
        for i=baslangic:son
    
     k=strfind(C{1}{i},'kt2=');
     if(k==1)
        kt2_n = textscan(C{1}{i},'kt2=%f');
        kt2_n = kt2_n{1};
        kt2_n_char = C{1}{i};
     end
        end  
        fclose(file);
%----

% [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
% [ V, I ] = simulation( ocean_address_vgs, result_address_vgs );
% [ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds );
% vth = importdata('./sonuclar/vth.dat');
% vdsat = importdata('./sonuclar/vdsat.dat');

 end    
%%%%%%% UC1 %%%%%%%
 if vbs < 0 && level == 2
[Ids_boyut,~] = size(Ids);
%ratio=Ids(Ids_boyut,1)/Ids_app;
%ratio_vdsat = vdsat/Vdsat_app;

fark=Ids(Ids_boyut,1)-Ids_app;
%fark_vdsat = vdsat-Vdsat_app;

syms uc1_n

if mobmod==1
Tmp = (((ua_n+ua1_n*(c_temp-1))+(uc_n+uc1_n*(c_temp-1))*Vbseff)*((Vgsteff+2*vth_app)/tox)+(ub_n+ub1_n*(c_temp-1))*(((Vgsteff+2*vth_app)/tox)^2));
elseif mobmod==2
Tmp = (((ua_n+ua1_n*(c_temp-1))+(uc_n+uc1_n*(c_temp-1))*Vbseff)*((Vgsteff)/tox)+(ub_n+ub1_n*(c_temp-1))*(((Vgsteff)/tox)^2));    
elseif mobmod==3
Tmp = (((ua_n+ua1_n*(c_temp-1))*((Vgsteff+2*vth_app)/tox)+(ub_n+ub1_n*(c_temp-1))*(((Vgsteff+2*vth_app)/tox)^2))*(1+(uc_n+uc1_n*(c_temp-1))*Vbseff));    
end

    denominator = 1+Tmp;

mu = (u0_n*((c_temp)^ute_app))*10^-4/denominator;

Esat = 2*Vsat/mu;

Abulk_temp = (1+(k1_n/(2*sqrt(phi_s-Vbseff)))*(((a0*Leff)/(Leff+2*sqrt(xj*Xdep)))*(1-ags*Vgsteff*(Leff/(Leff+2*sqrt(xj*Xdep)))^2)+b0/(Weff+b1)));
if Abulk_temp >= 0.1
Abulk = Abulk_temp*(1/(1+keta*Vbseff));
elseif Abulk_temp <= 0.1
Abulk =((0.2-Abulk_temp)/(3-20*Abulk_temp))*(1/(1+keta*Vbseff));
end


lambda = a1*Vgsteff + a2;
a = Abulk^2*Weff*Vsat*Cox*Rds + (1/lambda-1)*Abulk;
b = -((Vgsteff+2*Vt)*(2/lambda-1) + Abulk*Esat*Leff + 3*Abulk*(Vgsteff+2*Vt)*Weff*Vsat*Cox*Rds);
c = (Vgsteff+2*Vt)*Esat*Leff + 2*((Vgsteff+2*Vt)^2)*Weff*Vsat*Cox*Rds;

Vdsat = (-b-sqrt(b^2-4*a*c))/(2*a);
%Vdsat = 0.863;
Vdseff = Vdsat - 0.5*(Vdsat-vds-delta_n+sqrt((Vdsat-vds-delta_n)^2+4*delta_n*Vdsat));
    
Litl = sqrt((Esi*tox*xj)/Eox);

Vasat = (Esat*Leff + Vdsat + 2*Rds*Vsat*Cox*Weff*Vgsteff*(1-(Abulk*Vdsat)/(2*(Vgsteff+2*Vt))))/(2/lambda - 1 + Rds*Vsat*Cox*Weff*Abulk);

Vascbe = ((Leff/pscbe2)*exp((pscbe1*Litl)/(vds-Vdseff)));

theta_rout = pdiblc1*(exp(-drout*Leff/(2*lto))+2*exp(-drout*Leff/lto)) + pdiblc2;
Vadiblc = ((Vgsteff+2*Vt)/(theta_rout*(1+pdiblcb*Vbseff)))*(1-(Abulk*Vdsat)/(Abulk*Vdsat + Vgsteff + 2*Vt));
Vaclm = (Abulk*Esat*Leff+Vgsteff)*(vds-Vdseff)/(pclm*Abulk*Esat*Litl);
Va = Vasat + (1+(pvag*Vgsteff)/(Esat*Leff))*((Vaclm*Vadiblc)/(Vaclm+Vadiblc));
Idso = (Weff*mu*Cox*Vgsteff*(1-(Abulk*Vdseff)/(2*(Vgsteff+2*Vt)))*Vdseff)/(Leff*(1+Vdseff/(Esat*Leff)));

Ids_app2 = fark+((Idso/(1+Rds*Idso/Vdseff))*(1+(vds-Vdseff)/Va)*(1+(vds-Vdseff)/Vascbe));

eqn = Ids_app2 - Ids_aim(Ids_boyut,1);
uc1_app = double(solve(eqn,uc1_n));

    if isempty(uc1_app)
        hata = 1;
       return
    end
    uc1_n = uc1_app

    
    %%%%%%%%% UB1 in UC1 %%%%%%%%%%%    
   while max(abs(diff(I)))>mean(abs(diff(I((size(I,1)-20):size(I,1)))))*2
        if abs(ub1_n)<abs(ub1_n_max) 
        min_diff = min(diff(I,2));
        ub1_n = ub1_n+(ub1_n_max-ub1_n_def)/step; 
        
    % Değişkenlerin tekrardan karaktere dönüştürülmesi
    ub1_n_char_new = sprintf('ub1=%d', ub1_n);
    ub1_n_char_old = ub1_n_char;
    
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması  
    execute = sprintf(command,ub1_n_char,ub1_n_char_new,model_address);
    system(execute);
    
%% Parametrelerin dosyadan tekrar okunması (hata almamak için)

file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);

for i=baslangic:son
    
     k=strfind(C{1}{i},'ub1=');
     if(k==1)
        ub1_n = textscan(C{1}{i},'ub1=%f');
        ub1_n = ub1_n{1};
        ub1_n_char = C{1}{i};
     end
end  
fclose(file);
ub1_n

%%

    %Yeni parametreler ile analiz
    ocean_command = sprintf('ocean -restore %s',ocean_address_vgs); 
    unix(ocean_command)
       
    result = importdata(result_address_vgs);
    new_data = result.data;
    V = new_data(:,1);
    I = new_data(:,2);

    %%
   
        if min(diff(I,2))<min_diff
           % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
           execute = sprintf(command,ub1_n_char_new,ub1_n_char_old,model_address);
           system(execute);
           break            
        end
       if abs(ub1_n)>=abs(ub1_n_max)
            break
       end     
    end

    
    end   
%%%%%%%%%%%%%%%%%%%%
end


%---- UTE Dosyaya yazdırma
    ute_n_char_new = sprintf('ute=%d', ute_n);
    ute_n_char_old = ute_n_char;
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması 
    execute = sprintf(command,ute_n_char,ute_n_char_new,model_address);
    system(execute);   
    %% Parametrelerin dosyadan tekrar okunması (hata almamak için)
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file,formatSpec);
        for i=baslangic:son
    
     k=strfind(C{1}{i},'ute=');
     if(k==1)
        ute_n = textscan(C{1}{i},'ute=%f');
        ute_n = ute_n{1};
        ute_n_char = C{1}{i};
     end
        end  
        fclose(file);
%----


%---- UC1 Dosyaya yazdırma
    uc1_n_char_new = sprintf('uc1=%d', uc1_n);
    uc1_n_char_old = uc1_n_char;
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması 
    execute = sprintf(command,uc1_n_char,uc1_n_char_new,model_address);
    system(execute);   
    %% Parametrelerin dosyadan tekrar okunması (hata almamak için)
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file,formatSpec);
        for i=baslangic:son
    
     k=strfind(C{1}{i},'uc1=');
     if(k==1)
        uc1_n = textscan(C{1}{i},'uc1=%f');
        uc1_n = uc1_n{1};
        uc1_n_char = C{1}{i};
     end
        end  
        fclose(file);
%----

%--------------------------------------------------------------------------
%end

end


