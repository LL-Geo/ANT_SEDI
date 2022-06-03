function [wflx_v c] = asthen_viscous(it,wflx,wflx_v,icenode,c,cold,bo)


wflx_v(it,:) = -wflx(it,:);

% INITIALIZING INPUT PARAMETERS
nnode = icenode; %number of spatial steps
%ntime = 120; %number of time steps
delx = 1.0e4; %spatial step
delt = 100; %time step (unit in yr)
Da = 5e7; % (unit in m^2/yr)
%xinit = 0; %initial x position
%co = 10000.0;
%t = 0.0;

beta = Da*delt/(delx^2); %Diffusion term multiplier

%x = (xinit:delx:((delx*nnode)-delx))'; %distance vector

% Graphing Variables to take slices of output
iprint = 5;
ipr = 0;
ic = 0;

%Initial thickness of asthenosphere (c)
%for i = 1:nnode
    %c(i,1) =  10000.0;
    %c(i,1) =  0.0;
    %cold(i,1) = c(i,1);
%end

%bo = 0.0*ones(nnode,1);
% form a & b vector matrix
%for it=1:ntime
    %t = t+delt;
    
    for i=1:nnode-1
        for j=1:3
            a(i,j)=0.0;
        end
        b(i)=0.0;
    end
    
    
    for i=1:nnode-1
        % central difference
        a(i,1) = -beta;
        a(i,2)= 2*beta + 1;
        a(i,3)= - beta;
        
        if(i>1) && (i<nnode-1)
            b(i) =  c(i,1) + beta*(-wflx_v(it,i-1)+2*wflx_v(it,i)-wflx_v(it,i+1))+beta*(bo(i-1,1)-2*bo(i,1)+bo(i+1,1));
        end
        
        if(i==1)
            b(i) = c(i,1);
        end
        
        if(i==nnode-1)
            b(i) = c(i,1);
        end
        
    end
    
    % enforce upstream no flux bc
    
    a(1,1) = 0.0;
    a(1,2) = 1+beta;
    
    %   % enforce  down stream no flux bc
    %   a(nnode,2) = 1+beta;
    %   a(nnode,3) = 0.0;
    
    
    % enforce  down stream fixed elevation
    
    b(nnode-1) = b(nnode-1) - a(nnode-1,3)*c(nnode,1);
    a(nnode-1,3) = 0;
    
    aa = a;
    
    
    % solve matrix using tridiagonal solver
    
    
    b(1)=b(1)/a(1,2);
    a(1,2)=a(1,3)/a(1,2);
    
    for  n=2:nnode-1;
        ww=a(n,2)-a(n,1)*a(n-1,2);
        a(n,2)=a(n,3)/ww;
        b(n)=(b(n)-a(n,1)*b(n-1))/ww;
    end
    n=nnode-2;
    
    for i=1:nnode-1
        if(n>0)
            b(n)=b(n)-a(n,2)*b(n+1);
        end
        n = n -1;
    end
    
    
    for i=1:nnode-1
        c(i,1) = b(i);
    end
    
    
    %c(nnode,1) = 10000;
    
    for n=1:nnode
        uplift(n,1) = (c(n,1)-cold(n,1))/delt;
        cold(n,1) = c(n,1);
    end
    
    
    
    %fb(it,1) = max(c(:,1) - 10000); % forebulge thickness calculation
    
    
    % store plot array
    
    
    % store plot array
    
    ipr = ipr + 1;
    if(ipr>=iprint)
        ipr = 0;
        ic = ic+1;
        concplot(:,ic) = c(:,1);
        up_plot(:,ic) = uplift(:,1);
        
    end
    
    % c(nnode,1) = c(nnode-1,1);
    % c(1,1) = c(2,1);
    
    % store displacement for each time step
    %dis_m(it,:) = c(:,1);

    
%end

end