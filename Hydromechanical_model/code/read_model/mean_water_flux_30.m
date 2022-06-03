function [vz,Mask]=mean_water_flux_30(delt3,stime,ttime,s,filename)

%% ICE
%load ice sheet thickness
% ----------------------------------------------------
xls(1) = 0;
delx = 16000;
icenode = 101;
ntime = 300;% number of time steps
delt = 100;% time steop is 100 yr;
delt0 = 0; % unit in years (0)
delt1 = 19000;% unit in years 
delt2 = 20000;% unit in years
Ho = 3000.0; % max ice sheet thickness, unit in m
Lo = 1300000; % max ice sheet length, unit in m
time = 0;
for n=2:icenode
    xls(n) = xls(n-1)+delx;
end

for it = 1:ntime
    
    time = time + delt;
    
    % Dynamic Ice Sheet Model; Length and Max-thickness Calculation
    
    if  (time<=delt0)
        Lt = 1.0;
        Ht = 0.0;
        
    end
    
    if(time>delt0) && (time<=delt1);
        
        
        Lsl1 = Lo/(delt1-delt0);
        Hsl1 = Ho/(delt1-delt0);
        
        Lt = Lsl1*(time-delt0);
        Ht = Hsl1*(time-delt0);
    end
    
    if(time>delt1) && (time<=delt2)
        Lt = Lo;
        Ht = Ho;
    end
    
    if(time>delt2) && (time<=delt3)
        Lt = Lo;
        Ht = Ho;
        Lsl2 = Lo/(delt3-delt2);
        Hsl2 = Ho/(delt3-delt2);
        Lt = Lo - Lsl2*(time-delt2);
        Ht = Ho - Hsl2*(time-delt2);
    end
    
    if (time>=delt3)
        Lt = 1;
        Ht = 0;
    end

  
    % These commands are for the thickness evolution calculation
    
    for n=1:icenode
        arg = abs( 1.0 - ( xls(n)/Lt )^2 );
        nu_i(it,n) = Ht*( arg^0.5 );
 
        nun(n)= Ht*( arg^0.5 );
        if(xls(n)>Lt)
            nu_i(it,n) = 0.0;
            nun(n) = 0.0;
        end
        xx_ice(it,n) = xls(n);
    end
    
    
end


ice_thick=nu_i';

Mask=ice_thick-ice_thick;
Mask(ice_thick>0)=1;
Mask(Mask==0)=nan;

% Mask2=ice_thick-ice_thick;
% 
% Mask2(15:35,:)=1;
% Mask2(Mask2==0)=nan;
% 
% Mask=Mask.*Mask2;

%%
fid = fopen(filename,'rt');
out_9_f=[];
fgetl(fid);
fgetl(fid);
for i=1:300
    fgetl(fid);
    file = textscan(fid,[repmat('%f',1,s)]);
    kk=cell2mat(file);
    out_9_f=[out_9_f ;kk];
end

fclose(fid);

vz=out_9_f(:,end);
tem=out_9_f(:,2);
vz=reshape(vz,101,299);
tem=reshape(tem,101,299);

Mask(:,300)=[];
vz_m=vz.*Mask;

kk1=mean(vz_m,'omitnan')*1000;

TT=stime:1:ttime;

plot(TT,kk1(stime:ttime));
xlim([stime,ttime]);
hold on