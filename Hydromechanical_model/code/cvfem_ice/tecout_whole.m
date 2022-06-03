ni_1 = 1;
nj_1 = cols+2; % 71;
nk_1 = cols+1;
ni_2 = 1;
nj_2 = 2; % 70;
nk_2 = cols+2; % 71;
 
% 
m = 0;
% 
nelem = (rows-1)*(cols-1);
nnode = rows*cols;
% 
for nn_m=1:rows-1
    for n_m=1:cols-1
        
        m=m+1;
        ni(m,1) = ni_1;
        nj(m,1) = nj_1;
        nk(m,1) = nk_1;
        m=m+1;
        ni(m,1) = ni_2;
        nj(m,1) = nj_2;
        nk(m,1) = nk_2;
        
        % update nodes
        
        ni_1 = ni_1+1;
        nj_1 = nj_1+1;
        nk_1 = nk_1+1;
        
        ni_2 = ni_2+1;
        nj_2 = nj_2+1;
        nk_2 = nk_2+1;
        
    end
    ni_1 = ni_1+1;
    nj_1 = nj_1+1;
    nk_1 = nk_1+1;
    
    ni_2 = ni_2+1;
    nj_2 = nj_2+1;
    nk_2 = nk_2+1;
    
    
end
% 
% figure(5)
% tri=[ni,nj,nk];
% triplot(tri,x,y);
% 
% figure(6)
% trisurf(tri,x,y,phi);
% shading interp   % smooth it out
% view(0,90);
% colorbar;
% alpha(0.5);
% 
for it = 1:ntime
    for n=1:nnode
        dxx(n,it) = disx_t(it,n);
        dyy(n,it) = disy_t(it,n);
    end
end
 
output_tec_filename='CVFEM_stress_tec_whole.dat';
fid=fopen(output_tec_filename,'w');
%fprintf(fid,'VARIABLES = "X", "Y","sigx_t","sigy_t","sigxy_t","sigz_t","sig_kk_t","ux","uy","FC","phi_t","dsig_dt","pp","FC_p","matp","dFSM"\n');
fprintf(fid,'VARIABLES = "X", "Y","sigx_t","sigy_t","sigxy_t","sigz_t","sig_kk_t","ux","uy","FC","dsig_dt","matp","dFSM"\n'); 
nelem_tri = (cols-1)*(rows-1)*2;
xarray = 0;
elarray = 0;
 
for it = 1:ntime
    
    for n=1:nnode
        xarray(1,n)=dxx(n,it);
        xarray(2,n)=dyy(n,it);
        xarray(3,n) = sigx_t(n,it);
        xarray(4,n) = sigy_t(n,it);
        xarray(5,n) = sigxy_t(n,it);
        xarray(6,n) = sigz_t(n,it);
        xarray(7,n) = sig_kk_t(n,it);
        xarray(8,n) = -ux_t(n,it);
        xarray(9,n) = -uy_t(n,it);
        xarray(10,n) = fc(n,it);
        %xarray(10,n) = phi_t(n,it);
        xarray(11,n) = dsig_dt(n,it);
        %xarray(13,n) = qx_t(n,it);
        %xarray(14,n) = qy_t(n,it);
        %xarray(14,n) = fc_pp(n,it);
        xarray(12,n) = matp(n);
        xarray(13,n) = dFSM(n);
    end
    
    for m=1:nelem_tri;
        elarray(1,m)=ni(m);
        elarray(2,m)=nj(m);
        elarray(3,m)=nk(m);
    end
    
    fprintf(fid,'ZONE N= %d, E= %d, F=FEPOINT, ET=TRIANGLE\n',nnode,nelem_tri);
    fprintf(fid,'%6.3f %6.3f %12.6f %12.6f %12.6f %12.6f %12.6f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n',xarray);
    fprintf(fid,'%d %d %d\n',elarray);
end


