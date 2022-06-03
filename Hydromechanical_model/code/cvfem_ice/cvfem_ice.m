% CVFEM_ice
%--------------------------------------------------------------------------
% Simulate the disturbed stress within the lithosphere due to ice sheet
% loading on top, with viscous mantle flow along the bottom

% Simple nine node quads
% Asuming retangular domain and square elements
%--------------------------------------------------------------------------
% Authors: Yipeng Zhang-Mark Person-Vaughan Voller, 
% Affiliations: Earth and Enviromental Science Dept., New Mexico Tech
% Contact: zhangyipengjd@gmail.com, mark.person@nmt.edu
% 08/20/2018
%--------------------------------------------------------------------------
% Model Validation against Lund (2005) SKB report
% titled with "Effects of deglacition on the crustal stress field and
% implication for endglacial faulting: A parametric0 study of simple Earth 
% and ice models"

clear all
clc
close all
%GRID SET UP-----

%We assumme the global grid is set up and numbered as follows


%Figure

% 11*     12*     13*     14*      15*

%  6*      7*      8*      9*      10*

%  1*      2*      3*      4*       5*

%note we must have ODD numbers of rows and columns (cols)
%In the above we have 5 cols, 3 rows and 2 elements
%We are also currenly assuming that all elemnts are the same size
%BUT the x (DX) and y (DY) dimensions can differ.

% zero out coefficients


%domain size
xdim=2800000; % unit in m
ydim=100000;  % unit in m

cols=281;%901;%51;
rows=11;%11%7;


%size of elements
DX=xdim/(cols-1);
DY=ydim/(rows-1);

%NODES--following number system shown in Figure
nodes=rows*cols;
for ii=1:rows
    for jj=1:cols
        nval=(ii-1)*cols+jj;
        x(nval)=(jj-1)*DX;
        y(nval)=(ii-1)*DY;
        %         if ii==rows-3 %XXXXX EXTRA
        %             y(nval)=y(nval)-DY/2;
        %         end
        %         if ii==rows-2 %XXXXX EXTRA
        %             y(nval)=y(nval)-DY;
        %         end
        %         if ii==rows-1 %XXXXX EXTRA
        %             y(nval)=y(nval)-3*DY/2;
        %         end
        %         if ii==rows
        %             y(nval)=y(nval)-2*DY;
        %         end %XXXX end extra
    end
end

nret=0; % number of elements

%This code creates the array ret(elements,9) used to
%identify the global node with the local element node
%Our exampel figure has two elements
%elemnt 1 has global nodes 1,2,3,8,13,12,11,6,7
%note counterclockwise sence
%first node is always bottom left and last is center node

for ii=1:2:rows-1
    for jj=1:2:cols-1
        nret=nret+1;
        nval=(ii-1)*cols+jj;
        ret(nret,1)=nval;
        ret(nret,2)=nval+1;
        ret(nret,3)=nval+2;
        ret(nret,4)=nval+2+cols;
        ret(nret,5)=nval+2+2*cols;
        ret(nret,6)=nval+1+2*cols;
        ret(nret,7)=nval+2*cols;
        ret(nret,8)=nval+cols;
        ret(nret,9)=nval+1+cols;
        
    end
end

figure(1)
plot(x,y,'o')


nelem = (cols-1)*(rows-1)/4;

for m=1:nelem
    xc(m) = x(ret(m,9));
    yc(m) = y(ret(m,9));
end


figure(2)
plot(xc,yc,'o')

% find the centroid of the element and determine the material properties of
% the element
mat = zeros(1,nret);
matp = zeros(1,nodes);
% mat_h = zeros(1,nret_h);
% mat_hp = zeros(1,nodes_h);

for i = 1:nret
    xc(i) = x(ret(i,9));
    yc(i) = y(ret(i,9));
    if (yc(i)<80000)
        mat(i) = 2;% mat = 2 stands for deeper crust
    else
        mat(i) = 1;% mat = 1 stands for shallower crust
    end
    
    matp(ret(i,1)) = mat(i);
    matp(ret(i,2)) = mat(i);
    matp(ret(i,3)) = mat(i);
    matp(ret(i,4)) = mat(i);
    matp(ret(i,5)) = mat(i);
    matp(ret(i,6)) = mat(i);
    matp(ret(i,7)) = mat(i);
    matp(ret(i,8)) = mat(i);
    matp(ret(i,9)) = mat(i);
end


%EDGES/BOUNDaries
% the code below stores the domain boundaries
% as a set of two noded edges--numbered counter clockwise
% the last entry in the array e identifies a  specfic boundary
nedge=0;
%bound 1 (bottom)
for jj=1:cols-1;
    nedge=nedge+1;
    nval=jj;
    e(1,nedge)=nval;
    e(2,nedge)=nval+1;
    e(3,nedge)=1;
end
%bound 2 right
for ii=1:rows-1
    nedge=nedge+1;
    nval=ii*cols;
    e(1,nedge)=nval;
    e(2,nedge)=nval+cols;
    e(3,nedge)=2;
end
%bound 3 (top)
for jj=1:cols-1;
    nedge=nedge+1;
    nval=rows*cols-(jj-1);
    e(1,nedge)=nval;
    e(2,nedge)=nval-1;
    e(3,nedge)=3;
end
%bound 4 left
for ii=1:rows-1
    nedge=nedge+1;
    nval=1+(rows-1)*cols-(ii-1)*cols;
    e(1,nedge)=nval;
    e(2,nedge)=nval-cols;
    e(3,nedge)=4;
end


%plot nodes
figure (3)
hold on
for k=1:nret
    for ii=1:8
        
        plot([x(ret(k,ii)),x(ret(k,ii+1))],[y(ret(k,ii)),y(ret(k,ii+1))],'O')
    end
end
hold off
%--------------------------------------------------------------------------
%isoparametric element

%   7*        6*        5*

%   8*        9*        4*
%        |---------|
%   1*   |    2*   |    3*

%we assume a square isoparmetric element in psi, eta space
%the coordinates of local node 1 (bottom left) are psi=-1, eta=-1
%the coordinates of local node 5 (top right) are psi=1, eta =1

%Each node will have a number of sections of conrol volume faces
%For example local node 1,3, 5 and 7 have two face sections
%local nodes 2,4,6,8 have four face sections
%and the center node (psi=0, eta=0) node 9 has 8 sections.

%We create the folowing arrays psiver (9,10) and etaver(9,10) to store
%the vertices (in eta, psi space) of these face edges for each node.
%Considering the entries for node 2
%the firts entry psiver(2,1) gives us the number of face (4)
%the entries psiver(2,2) psiver(2,3), psiver(2,4), psiver(2,5), and
%psiver(2,6)then give the psi values at the vertices MOVING in a
%counter-clockwise direction with out interruption of an element edge
%we creat etaver in the same way
%Note extra entry in psiver(9,10), etaver(9,10) to close CV

psiver=[2,-0.5,-0.5,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4, 0.5, 0.5, 0.0,-0.5,-0.5, 0.0, 0.0, 0.0, 0.0;...
    2, 1.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4, 1.0, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0;...
    2, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4,-0.5,-0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0;...
    2,-1.0,-0.5,-0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4,-1.0,-0.5,-0.5,-0.5,-1.0, 0.0, 0.0, 0.0, 0.0;...
    8, 0.0, 0.5, 0.5, 0.5, 0.0,-0.5,-0.5,-0.5, 0.0];

etaver=[2,-1.0,-0.5,-0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4,-1.0,-0.5,-0.5,-0.5,-1.0, 0.0, 0.0, 0.0, 0.0;...
    2,-0.5,-0.5,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4, 0.5, 0.5, 0.0,-0.5,-0.5, 0.0, 0.0, 0.0, 0.0;...
    2, 1.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4, 1.0, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0;...
    2, 0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;...
    4,-0.5,-0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0;...
    8,-0.5,-0.5, 0.0, 0.5, 0.5, 0.5, 0.0,-0.5,-0.5];

%These arrays allow for a more comapct derivation of the
%discrete equations

%ALSO Note:
%Assumming a shifting that sets the x,y coordinates of the
% local node 1 to x=0, y=0 we have the following realtionship
%between the isoparemetric (psi, eta) and gloabl (x,y) positions

%psi=x/lenx-1; eta=y/leny-1                  EQ(1)

%where
%lenx=DX/2 and leny =DY/2                    EQ(2)
% YPZ: not correct lenx = DX, leny = DY;
%
%Then by the product rule we can write
% dT/dx = [dT/d_psi] X [d_psi.dx] = [dT/d_psi] X (1/lenx)     EQ(3)

%END GRID SETUP -------

%SET up Coefficents--like my Book
apx=zeros(nodes,1);
apy=zeros(nodes,1);
ax=zeros(nodes,nodes);
ay=zeros(nodes,nodes);
axs=zeros(nodes,nodes);
ays=zeros(nodes,nodes);

%source terms for boundary
BBx=zeros(nodes,1);
BCx=zeros(nodes,1);
BBy=zeros(nodes,1);
BCy=zeros(nodes,1);

%source term for domain
BBpx=zeros(nodes,1);
BBpy=zeros(nodes,1);
P=zeros(nodes,1);

ux=zeros(nodes,1);
uy=zeros(nodes,1);

%constants
E = 56e9; %Pa;  % Pa 1e9 Pa = 1 GPa
nu = 0.25;  %unitless;

%--------------------------------------------------------------------------
% switch on if it is plane strain problem
E=E/(1-nu^2); 
nu=nu/(1-nu);

%pesudo diffusivities
kxx= E/(1-nu.^2);%??? not sure what is this doing-- Yipeng Zhang on Jan23,2016
kxy=E/(2*(1+nu));
kyy=kxx;

kxxs=kxx*nu;
kxys=kxy;
kyys=kxxs;
%--------------------------------------------------------------------------
E_m(1) = E;% unit in Pa
E_m(2) = E;% unit in Pa

nu_m(1) = nu;
nu_m(2) = nu;


for ele=1:nret %loop on elemnts
    
    E = E_m(mat(ele));
    nu = nu_m(mat(ele));
    
    kxx= E/(1-nu.^2);
    kxy=E/(2*(1+nu));
    kyy=kxx;
    
    kxxs=kxx*nu;
    kxys=kxy;
    kyys=kxxs;
    
    %isoparmetric conversion see Eq(1) to EQ(3)
    lenx=abs(x(ret(ele,1))-x(ret(ele,5)))/2;
    leny=abs(y(ret(ele,1))-y(ret(ele,5)))/2;
    
    for node =1:9
        faces=psiver(node,1);
        for ii=1:faces
            
            %signed (x,y) lengths of faces using EQ(2)
            delx=(psiver(node,ii+2)-psiver(node,ii+1))*lenx;
            dely=(etaver(node,ii+2)-etaver(node,ii+1))*leny;
            
            psi=(psiver(node,ii+1)+psiver(node,ii+2))/2;
            %psi=psinode(node,1)
            eta=(etaver(node,ii+1)+etaver(node,ii+2))/2;
            
            
            %Derivatives of shape functions (in x,y)
            %see bi-quadractic shape functions Eq(3) see page
            %https://engineering.purdue.edu/~abe601/lecture/Ch_26.pdf
            
            Nx(1)=.25*(psi*eta+eta*(psi-1))*(eta-1)/lenx;
            Nx(3)=.25*(psi*eta+eta*(psi+1))*(eta-1)/lenx;
            Nx(5)=.25*(psi*eta+eta*(psi+1))*(eta+1)/lenx;
            Nx(7)=.25*(psi*eta+eta*(psi-1))*(eta+1)/lenx;
            Nx(2)=-psi*eta*(eta-1)/lenx;
            Nx(6)=-psi*eta*(eta+1)/lenx;
            Nx(4)=-.5*(psi+(psi+1))*(eta^2-1)/lenx;
            Nx(8)=-.5*(psi+(psi-1))*(eta^2-1)/lenx;
            Nx(9)=2*psi*(eta^2-1)/lenx;
            
            Ny(1)=.25*(psi*eta+psi*(eta-1))*(psi-1)/leny;
            Ny(3)=.25*(psi*eta+psi*(eta-1))*(psi+1)/leny;
            Ny(5)=.25*(psi*eta+psi*(eta+1))*(psi+1)/leny;
            Ny(7)=.25*(psi*eta+psi*(eta+1))*(psi-1)/leny;
            Ny(4)=-psi*eta*(psi+1)/leny;
            Ny(8)=-psi*eta*(psi-1)/leny;
            Ny(2)=-.5*(eta+(eta-1))*(psi^2-1)/leny;
            Ny(6)=-.5*(eta+(eta+1))*(psi^2-1)/leny;
            Ny(9)=2*eta*(psi^2-1)/leny;
            
            %XXX Set up diffusion Coeficients
            %XXTHIS is where you need to add modifications for linear
            %elastisity--follow the approach in previous traingle ele code
            
            nodeg=ret(ele,node);  %global node_number of current node
            
            for jj=1:9
                nodenb=ret(ele,jj); %global numbers of neigbouring nodes
                BBpx(nodeg,1)=BBpx(nodeg,1)-P(nodeg,1)*Nx(jj)*dely;
                BBpy(nodeg,1)=BBpy(nodeg,1)+P(nodeg,1)*Ny(jj)*delx;
                
                ax(nodeg,nodenb)=ax(nodeg,nodenb)+(kxx*Nx(jj)*dely-kxy*Ny(jj)*delx);
                axs(nodeg,nodenb)=axs(nodeg,nodenb)+(kxxs*Ny(jj)*dely-kxys*Nx(jj)*delx);
                ay(nodeg,nodenb)=ay(nodeg,nodenb)+(kxy*Nx(jj)*dely-kyy*Ny(jj)*delx);
                ays(nodeg,nodenb)=ays(nodeg,nodenb)+(kxys*Ny(jj)*dely-kyys*Nx(jj)*delx);
            end
            
        end
    end
end

for node=1:nodes
    apx(node,1)=-ax(node,node);
    ax(node,node)=0;
    apy(node,1)=-ay(node,node);
    ay(node,node)=0;
end

Ax=sparse(ax);  %stores coefficents as a sparse system
Axs=sparse(axs);
Ay=sparse(ay);
Ays=sparse(ays);


%Set Speficied Stress Boundary Conditions On Top
%--------------------------------------------------------------------------
%Ice sheet loading term evolving with time is added from 0 to 1000 m by
%Yipeng Zhang on 06/18/2015
%--------------------------------------------------------------------------
xls(1) = 0;
delx = 10000;
icenode = 281;
ntime = 300;% number of time steps
delt = 100;% time steop is 100 yr;
delt0 = 0; % unit in years (0)
delt1 = 19000;% unit in years 
delt2 = 20000;% unit in years
delt3 = 29000;% unit in years
Ho = 3000.0; % max ice sheet thickness, unit in m
Lo = 1300000; % max ice sheet length, unit in m
time = 0;


% specify lithosphere properties
Delx_w = 10000;
rhom = 3380.0;% unit in kg/m^3
grav = 9.812*(3600*24*365)*(3600*24*365);
h = ydim; % elastic thickness
% The definition of Rigidity (from Turcotte and Schubert 2002. Page 115)
D = E*h^3/12/(1-nu^2);
D = D*(3600.*24.*365.)*(3600.*24.*365.);
alpha_1 = (4*D/rhom/grav)^0.25;

for n=2:icenode
    xls(n) = xls(n-1)+delx;
end

% for n=2:icenode_h
%     xls_h(n) = xls_h(n-1)+delx_h;
% end

uy_t = zeros(nodes,ntime);
wflx = zeros(ntime,icenode);
dis_m = zeros(ntime,icenode);
wflx_v = zeros(ntime,icenode);
%c = zeros(ntime,icenode);


%--------------------------------------------------------------------------
% The following codes are for initilization of parameters of subroutine
% asthen_viscous
%--------------------------------------------------------------------------

%Initial thickness of asthenosphere (c)
nnode = icenode;
for i = 1:nnode
    %c(i,1) =  10000.0;
    c(i,1) =  0.0;
    cold(i,1) = c(i,1);
end

bo = 0.0*ones(nnode,1);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

for it = 1:ntime
    
    it
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
    
    



    if (it == 1)
        for n = 1:cols
            for m = 1:rows
                dnu_dt(n+281*(m-1),it)= (nu_i(it,n)-0)/delt;
            end
        end
    end

    
    if (it > 1)
        for n = 1:cols
            for m = 1:rows
                dnu_dt(n+281*(m-1),it)= (nu_i(it,n)-nu_i(it-1,n))/delt;
            end
        end
    end
    
       
    
    %----------------------------------------------------------------------
    % The following commands will compute displacement due to ice sheet
    % loading
    %----------------------------------------------------------------------
    
    % move the point loading line (xo = 0) from node 1 to node kk
    for kk = 1:icenode
        xo(kk) = xls(kk);
        % calcualte the solution in the x >= 0 domain using the load at
        % node kk at time it
        for nn=1:icenode
            xbar = abs(xls(nn)-xo(kk))/alpha_1;
            wi(kk,nn) = 910*grav*nu_i(it,kk)*Delx_w*alpha_1^3.0/(8.0*D) ...
                *exp(-xbar)*(cos(xbar)+sin(xbar));
        end
        % calcualte the solution in the x<= 0 domain using the load at
        % node kk at time it
        for jj = 1:icenode
            xbarm = abs(-xo(kk)-xls(jj))/alpha_1;
            wi_m(kk,jj) = 910*grav*nu_i(it,kk)*Delx_w*alpha_1^3.0/(8.0*D) ...
                *exp(-xbarm)*(cos(xbarm)+sin(xbarm));
        end
    end
    % summarize wi at node 1 to node mm from point loading at kk=1 to kk=nnode
    %(x>=0)and from the mirro image kk = 2 to i = nnode (x<0)
    
    for mm =1:icenode
        ws(it,mm) = sum(wi(1:icenode,mm))+sum(wi_m(2:icenode,mm));
        wflx(it,mm) = wflx(it,mm) + ws(it,mm);
    end
    %----------------------------------------------------------------------
    %Call function asthen_viscous to calculate the relative thickness "b"
    %and apply the "b" values to the lower boundary as specified
    %displacement boundary
    
    [wflx_v c] = asthen_viscous(it,wflx,wflx_v,icenode,c,cold,bo);
    
    % store displacement for each time step
    dis_m(it,:) = c(:,1);

%     % store displacement for each time step
%     dis_m(it,:) = c(:,1);
    %--------------------------------------------------------------------------
    % BBx=zeros(nodes,1);
    % BCx=zeros(nodes,1);
    % BBy=zeros(nodes,1);
    % BCy=zeros(nodes,1);
    BBx(1:nodes,1)=0;
    BCx(1:nodes,1)=0;
    BBy(1:nodes,1)=0;
    BCy(1:nodes,1)=0;
    
    
    for k=1:nedge
        nseg=e(3,k);%boundary segment number
        k1=e(1,k); %nodes on boundary seg
        k2=e(2,k);
        
        len=sqrt((x(k1)-x(k2)).^2+(y(k1)-y(k2)).^2);
        
        if nseg==1  % u and v fixed boundary
            %BCy(k1,1)=1e18;
            %BCy(k2,1)=1e18;
            
            %BBy(k1,1) = BBy(k1,1)+1e18*(wflx(it,k1));%*len/2;
            %BBy(k2,1) = BBy(k2,1)+1e18*(wflx(it,k2));%*len/2;
            BBy(k1,1) = BBy(k1,1)+1e18*(-dis_m(it,k1));%*len/2;
            BBy(k2,1) = BBy(k2,1)+1e18*(-dis_m(it,k2));%*len/2;

            
            BCy(k1,1) = BCy(k1,1)+1e18;%*len/2;
            BCy(k2,1) = BCy(k2,1)+1e18;%*len/2;
            
        end
        
        
        if nseg==2  % u and v fixed boundary
            BCx(k1,1)=1e18;
            BCx(k2,1)=1e18;
            %BCy(k1,1)=1e18;
            %BCy(k2,1)=1e18;
        end
        
        if nseg==4
            BCx(k1,1)=1e18;
            BCx(k2,1)=1e18;
            %BCy(k1,1)=1e18;  %only switch on if not exploiting symmetry
            %BCy(k2,1)=1e18;
            % impose a lateral stress term sig4 to the left boundary
            %BBx(k1,1)=BBx(k1,1)+sig4*len./2;
            %BBx(k2,1)=BBx(k2,1)+sig4*len./2;
        end
        
        if nseg==3
            
            %if (k1>=2812)&&(k1<=2801) % for small domain problem
            %             % impose a vertical stress term sig3 to the top boundary
            %
            BBy(k1,1)=BBy(k1,1)+910*9.81*nu_i(it,k1-2811)*len./2;
            BBy(k2,1)=BBy(k2,1)+910*9.81*nu_i(it,k2-2810)*len./2;
            
            %
            %end
            
            %       BBy(k1,1) = BBy(k1,1)+1e18*(wflx(it,k1-2811))*len/2;
            %       BBy(k2,1) = BBy(k2,1)+1e18*(wflx(it,k2-2810))*len/2;
            %BBy(k1,1) = 1e18*(wflx(it,k1-2811));
            %BBy(k2,1) = 1e18*(wflx(it,k2-2810));
            
            %       BCy(k1,1)=BCy(k1,1)+1e18*len/2;
            %       BCy(k2,1)=BCy(k2,1)+1e18*len/2;
            %BCy(k1,1)=1e18;
            %BCy(k2,1)=1e18;
        end
              
    end
    
    ux=zeros(nodes,1);
    uy=zeros(nodes,1);
    
    %Solver
    for iter=1:7000 %50000 %simple Jacobi-solver--may need more
        %
        
        
        for init=1:10
            RHS=Ax*ux+Axs*uy;
            RHS=RHS+BBx+BBpx;
            adiv=apx+BCx;
            ux=ux+.76*(RHS-adiv.*ux)./adiv;% the relaxation coefficient was changed from 0.95 to 0.76 due to refined grid
        end
        
        for init=1:10
            RHS=Ay*uy+Ays*ux;
            RHS=RHS+BBy+BBpy;
            adiv=apy+BCy;
            uy=uy+.76*(RHS-adiv.*uy)./adiv; % the relaxation coefficient was changed from 0.95 to 0.76 due to refined grid
        end
    end
    
    ux_t(:,it) = ux;
    uy_t(:,it) = uy;
    
    % select severl nodes underneath the ice sheet max, right in front of the
    % ice sheet and far away from the ice sheet and store the vertical
    % displacement rate for each time step
    uy_t2811(it) = uy(2811)/delt;
    uy_t2961(it) = uy(2961)/delt;
    uy_t3031(it) = uy(3031)/delt;
    
    %plot nodes
    % disx=x-10*ux';
    % disy=y-10*uy';
    
    disx=x-ux';
    disy=y-uy';
    
    disx_t(it,:) = disx;
    disy_t(it,:) = disy;
    
    %     figure (2)
    %     hold on
    %     for k=1:nret
    %         for ii=1:8
    %             plot([disx(ret(k,ii)),disx(ret(k,ii+1))],[disy(ret(k,ii)),disy(ret(k,ii+1))],'-O')
    %         end
    %     end
    
    %--------------------------------------------------------------------------
    % Tbe following codes compute the stresses given calculated displacement
    %--------------------------------------------------------------------------
    
    %isoparametric element
    
    %   7*        6*        5*
    
    %   8*        9*        4*
    %        |---------|
    %   1*   |    2*   |    3*
    
    %we assume a square isoparmetric element in psi, eta space
    %the coordinates of local nde 1 (bottom left) are psi=-1, eta=-1
    %the coordinates of local nde 5 (top right) are psi=1, eta =1
    
    %Each nde will have a number of sections of conrol volume faces
    %For example local nde 1,3, 5 and 7 have two face sections
    %local nodes 2,4,6,8 have four face sections
    %and the center nde (psi=0, eta=0) nde 9 has 8 sections.
    
    
    psi_el = [-1,0,1,1,1,0,-1,-1,0];
    eta_el = [-1,-1,-1,0,1,1,1,0,0];
    
    
    %for ele = 1:nret
    %E = E_m(mat(ele));
    %nu = nu_m(mat(ele));
    
    % arg1 = E/(1-nu^2);
    % arg2 = nu;
    % arg3 = E/2/(1+nu);
    %end
    
    nelem = (cols-1)*(rows-1)/4;
    
    sigx(1:nodes) = 0;
    sigxy(1:nodes) = 0;
    sigy(1:nodes) = 0;
    sigz(1:nodes) = 0;
    IC(1:nodes)=0;
    
    
    area_nd(1:nodes) = 0.0;
    
    for m=1:nelem
        
        %E = E_m(mat(m));
        %nu = nu_m(mat(m));
        
        % Convert the plane strain coefficient to plane stress conficient
        
        %E = E/(1-nu^2);
        %nu = nu/(1-nu);
        
        arg1 = E/(1-nu^2);
        arg2 = nu;
        arg3 = E/2/(1+nu);
        
        
        nde(1)=ret(m,1);
        nde(2)=ret(m,2);
        nde(3)=ret(m,3);
        nde(4)=ret(m,4);
        nde(5)=ret(m,5);
        nde(6)=ret(m,6);
        nde(7)=ret(m,7);
        nde(8)=ret(m,8);
        nde(9)=ret(m,9);
        
        lenx=abs(x(ret(m,1))-x(ret(m,5)))/2;
        leny=abs(y(ret(m,1))-y(ret(m,5)))/2;
        
        area=lenx*leny;
        
        
        for l=1:9
            %Calculate deritvete of shape functions at each node in an element
            
            psi=psi_el(l);
            eta=eta_el(l);
            
            
            %Derivatives of shape functions (in x,y)
            %see bi-quadractic shape functions Eq(3) see page
            %https://engineering.purdue.edu/~abe601/lecture/Ch_26.pdf
            
            Nx_s(1)=.25*(psi*eta+eta*(psi-1))*(eta-1)/lenx;
            Nx_s(3)=.25*(psi*eta+eta*(psi+1))*(eta-1)/lenx;
            Nx_s(5)=.25*(psi*eta+eta*(psi+1))*(eta+1)/lenx;
            Nx_s(7)=.25*(psi*eta+eta*(psi-1))*(eta+1)/lenx;
            Nx_s(2)=-psi*eta*(eta-1)/lenx;
            Nx_s(6)=-psi*eta*(eta+1)/lenx;
            Nx_s(4)=-.5*(psi+(psi+1))*(eta^2-1)/lenx;
            Nx_s(8)=-.5*(psi+(psi-1))*(eta^2-1)/lenx;
            Nx_s(9)=2*psi*(eta^2-1)/lenx;
            
            Ny_s(1)=.25*(psi*eta+psi*(eta-1))*(psi-1)/leny;
            Ny_s(3)=.25*(psi*eta+psi*(eta-1))*(psi+1)/leny;
            Ny_s(5)=.25*(psi*eta+psi*(eta+1))*(psi+1)/leny;
            Ny_s(7)=.25*(psi*eta+psi*(eta+1))*(psi-1)/leny;
            Ny_s(4)=-psi*eta*(psi+1)/leny;
            Ny_s(8)=-psi*eta*(psi-1)/leny;
            Ny_s(2)=-.5*(eta+(eta-1))*(psi^2-1)/leny;
            Ny_s(6)=-.5*(eta+(eta+1))*(psi^2-1)/leny;
            Ny_s(9)=2*eta*(psi^2-1)/leny;
            
            for i = 1:9
                sigx(nde(l)) =  sigx(nde(l)) + arg1*Nx_s(i)*ux(nde(i))*area;
                sigx(nde(l)) =  sigx(nde(l))  + arg1*arg2*Ny_s(i)*uy(nde(i))*area;
                sigxy(nde(l)) = sigxy(nde(l)) + arg3*Nx_s(i)*uy(nde(i))*area;
                sigxy(nde(l)) = sigxy(nde(l)) + arg3*Ny_s(i)*ux(nde(i))*area;
                sigy(nde(l)) =  sigy(nde(l)) + arg1*arg2*Nx_s(i)*ux(nde(i))*area;
                sigy(nde(l)) =  sigy(nde(l)) + arg1*Ny_s(i)*uy(nde(i))*area;
                sigz(nde(l)) = sigz(nde(l)) + arg1*arg2*(Nx_s(i)*ux(nde(i))+Ny_s(i)*uy(nde(i)))*area;
                
            end
            
            area_nd(nde(l)) = area_nd(nde(l)) + area;
            
        end
        
        
        
        
    end
    
    for n=1:nodes
        sigx(n) = sigx(n)/area_nd(n);
        sigxy(n) = sigxy(n)/area_nd(n);
        sigy(n) = sigy(n)/area_nd(n);
        sigz(n) = sigz(n)/area_nd(n);
    end
    
    
    % Calculate and store the stress changes through time
    for n=1:nodes
        sigx_t(n,it) = sigx(n);
        sigxy_t(n,it) = sigxy(n);
        sigy_t(n,it) = sigy(n);
        sigz_t(n,it) = sigz(n);
        sig_kk_t(n,it) =(sigx_t(n,it) + sigy_t(n,it) + sigz_t(n,it))/3;
       
        % Failure Criteria adopted from Ge et al (2009)
        % using sigx and sigy and sigxy to calculate maximum and minimum
        % principle stress(sig_1 and sig_3).
        % then use sig_1 and sig_3 to calculate the normal stress and shear
        % stress on a fault plane that has a 30 degrees with the vertical
        % direction
        
        % Principle stress in 3D
        % https://en.wikiversity.org/wiki/Principal_stresses
        %if (it>20) % This prevent variable "input" becoming NaN for ice free period
            I_1(n) = sigx(n) + sigy(n) + sigz(n);
            I_2(n) = sigx(n)*sigy(n) + sigy(n)*sigz(n) + sigz(n)*sigx(n) - sigxy(n)^2;
            I_3(n) = sigx(n)*sigy(n)*sigz(n) - sigz(n)*sigxy(n)^2;
            input(n) = (2*I_1(n)^3 - 9*I_1(n)*I_2(n) + 27*I_3(n))/(2*(I_1(n)^2-3*I_2(n))^(3/2));
            fi(n) = (1/3)*acos(input(n));
            % max principal stress
            sig_1(n,it) = I_1(n)/3+(2/3)*sqrt(I_1(n)^2-3*I_2(n))*cos(fi(n));
            % intermediate stress
            sig_2(n,it) = I_1(n)/3+(2/3)*sqrt(I_1(n)^2-3*I_2(n))*cos(fi(n)-(2*pi)/3);
            % min principal stress
            sig_3(n,it) = I_1(n)/3+(2/3)*sqrt(I_1(n)^2-3*I_2(n))*cos(fi(n)-(4*pi)/3);
        %else
%             sig_1(n,it) = 0.0;
%             sig_2(n,it) = 0.0;
%             sig_3(n,it) = 0.0;
%        end
        
        % max principle stress
        %sig_1(n,it) = 0.5*(sigx(n)+sigy(n))+sqrt((0.5*(sigx(n)-sigy(n)))^2+sigxy(n)^2);
        
        % min principle stress
        %sig_3(n,it) = 0.5*(sigx(n)+sigy(n))-sqrt((0.5*(sigx(n)-sigy(n)))^2+sigxy(n)^2);
        
        % change in fault stability margin (dFSM) (see Wu and Hasegawa,
        % 1996)
        mu = 0.6;
        beta = sin (atan(mu))/2/mu;% for optimally oriented faults
        % change in Fault Stability Margin (dFSM)
        dFSM(n,it) = mu*beta*(sig_1(n,it)+sig_3(n,it)) + 0.5*(sig_3(n,it)-sig_1(n,it));
        
        % Coulomb Stress (see Ge et al. (2009)), assuming a reverse
        % faulting stress state (theta = pi/6, theta is the angle between
        % the fault plane and sig_1 (within the ice sheet max extension sig_1 is 
        % essentially sig_Hmax, and essentailly horizontal , beyond the ice sheet
        % max extension, sig_1 is essentially sig_v, and essentially vertical);
        sign_p(n,it) = 0.5*(sig_1(n,it)+sig_3(n,it))-0.5*(sig_1(n,it)-sig_3(n,it))*cos(pi/3);
        tau_p(n,it) = 0.5*(sig_1(n,it)-sig_3(n,it))*sin(pi/3);
        fc(n,it) = tau_p(n,it) -0.6*sign_p(n,it);
      end    

%--------------------------------------------------------------------------    
    if (it == 1)
        for n = 1:nodes
            dsig_dt(n,1) = (sig_kk_t(n,1) - 0)/delt;
        end
    end
    
    if (it>1)
        for n = 1:nodes
            dsig_dt(n,it) = (sig_kk_t(n,it) - sig_kk_t(n,it-1))/delt;
        end
    end
    
end


%--------------------------------------------------------------------------
% The following codes generate the output file for tecplot for the whole
% grid
%--------------------------------------------------------------------------

tecout_whole

%--------------------------------------------------------------------------
% The following codes generate the output variable of temperature for 
% tecplot for the upper 20 km grid
%--------------------------------------------------------------------------

tecout_temp

%--------------------------------------------------------------------------
% The following codes generate the output file for tecplot for the mean 
% normal stress of the whole grid
%--------------------------------------------------------------------------
tecout_rift2d

%--------------------------------------------------------------------------

