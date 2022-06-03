ni_1t = 1;
nj_1t = 161+2; % 71;
nk_1t = 161+1;
ni_2t = 1;
nj_2t = 2; % 70;
nk_2t = 161+2; % 71;

ni_1tr = 2249;
nj_1tr = cols+2250; % 71;
nk_1tr = cols+2249;
ni_2tr = 2249;
nj_2tr = 2250; % 70;
nk_2tr = cols+2250; % 71;
% 
m = 0;
% 
nelemt = (3-1)*(161-1)*2;
nnodet = 3*161;

for nn_m=1:3-1
    for n_m=1:161-1
        
        m=m+1;
        nit(m,1) = ni_1t;
        njt(m,1) = nj_1t;
        nkt(m,1) = nk_1t;
        m=m+1;
        nit(m,1) = ni_2t;
        njt(m,1) = nj_2t;
        nkt(m,1) = nk_2t;
        
        % update nodes
        
        ni_1t = ni_1t+1;
        nj_1t= nj_1t+1;
        nk_1t = nk_1t+1;
        
        ni_2t = ni_2t+1;
        nj_2t = nj_2t+1;
        nk_2t = nk_2t+1;
        
    end
    ni_1t = ni_1t+1;
    nj_1t = nj_1t+1;
    nk_1t = nk_1t+1;
    
    ni_2t = ni_2t+1;
    nj_2t = nj_2t+1;
    nk_2t = nk_2t+1;
    
    
end


for m=1:640;
    elarray_elhom(1,m)=nit(m);
    elarray_elhom(2,m)=njt(m);
    elarray_elhom(3,m)=nkt(m);
end

% i= 0;
% ndarray = zeros(161*3,1);
% n_first = 2249;
% for r = 1:2
%     
%     for c = 1:161
%         i = i+1;
%         ndarray(i,1) = n_first + c-1;
%     end
% 
%     n_first = 2249 + r*281;
%         
% end

ndarray = zeros(3,161);

for j = 1:size(ndarray,2);
    ndarray(1,j) = 2249+j-1;
end
for i = 1:size(ndarray,1);
    ndarray(i,:)  = ndarray(1,:)+281*(i-1);
end

ndarray = ndarray';
ndarray = ndarray(:);

for it = 1:ntime
    for n=1:nnodet
        dxx(ndarray(n),it) = disx_t(it,ndarray(n));
        dyy(ndarray(n),it) = disy_t(it,ndarray(n));
    end
end

output_tec_filename='CVFEM_Rift2D_tec_LongIS_100k_Temp_0.dat';
fid=fopen(output_tec_filename,'w');
fprintf(fid,...
'VARIABLES = "X","Y","dX","dY","sig_mean_t","dsig_dt","ux","uy","sigxx","sigyy","sig_1","sig_3","temp","fc","sign","tau","dnu_dt"\n'); 

for it = 1:ntime
    
    for n = 1:nnodet
        xarrayt(1,n) = x(ndarray(n));
        xarrayt(2,n) = y(ndarray(n));
        xarrayt(3,n) = dxx(ndarray(n),it);
        xarrayt(4,n) = dyy(ndarray(n),it);
        xarrayt(5,n) = sig_kk_t(ndarray(n),it);
        xarrayt(6,n) = dsig_dt(ndarray(n),it);
        xarrayt(7,n) = -ux_t(ndarray(n),it);
        xarrayt(8,n) = -uy_t(ndarray(n),it);
        xarrayt(9,n) = sigx_t(ndarray(n),it);
        xarrayt(10,n) = sigy_t(ndarray(n),it);
        xarrayt(11,n) = sig_1(ndarray(n),it);
        xarrayt(12,n) = sig_3(ndarray(n),it);
%         xarrayt(13,n) = temp(ndarray(n),it);
        xarrayt(13,n) = 0;
        xarrayt(14,n) = fc(ndarray(n),it);
        xarrayt(15,n) = sign_p(ndarray(n),it);
        xarrayt(16,n) = tau_p(ndarray(n),it);
%         xarrayt(17,n) = dnu_dt(ndarray(n),it);
    end
    fprintf(fid,'ZONE N= %d, E= %d, F=FEPOINT, ET=TRIANGLE\n',nnodet,nelemt);
    fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n',xarrayt);
    fprintf(fid,'%d %d %d\n',elarray_elhom);
end

