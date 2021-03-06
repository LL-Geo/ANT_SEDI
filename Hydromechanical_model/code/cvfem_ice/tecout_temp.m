
% Read in the temperature profiles from each column

sizeA = [76 Inf];
for i =  1:55
fid = fopen(['temp_subsurface_col',num2str(i),'.txt'],'r');
T(:,:,i) = fscanf(fid,'%f',sizeA);
fclose(fid);
end

% Create the distance array 

for i = 1:1381
    dist(i) = 0.13e6/130*(i-1);
end

for i = 1:76
    depth(i) = 10*(i-1);
end

t = zeros(76,1000,1601);

% Interpret the temperature in between adjacent columns linearly
for i = 1:50
  delt_dis = dist(1+(i-1)*26) - dist(27+(i-1)*26);%?
    
  for j = 1:76
      
      for k = 1:1000
          delt_T = T(j,ceil(k/10),i) - T(j,ceil(k/10),i+1);
          slop =  delt_T/delt_dis;
          for n = (2+26*(i-1)):(26+26*(i-1))
              t(j,k,n) = T(j,ceil(k/10),i) + slop*(dist(n)-dist(1+(i-1)*26));
          end
      end
      
  end
end

for i = 51:54
    
    delt_dis1 = dist(1301+20*(i-51)) - dist(1321+20*(i-51));
    
  for j = 1:76
      
      for k = 1:1000
          delt_T = T(j,ceil(k/10),i) - T(j,ceil(k/10),i+1);
          slop =  delt_T/delt_dis;
          for n = (1302+20*(i-51)):(1320+20*(i-51))
              t(j,k,n) = T(j,ceil(k/10),i) + slop*(dist(n)-dist(1301+20*(i-51)));
          end
      end
      
  end
  
end
  
%   t(:,:,1) = T(:,:,1);
%   t(:,:,27) = T(:,:,2);
%   t(:,:,53) = T(:,:,3);
%   t(:,:,79) = T(:,:,4);
%   t(:,:,105) = T(:,:,5);
%   t(:,:,131) = T(:,:,6);
%   t(:,:,157) = T(:,:,7);
%   t(:,:,183) = T(:,:,8);
%   t(:,:,209) = T(:,:,9);
%   t(:,:,235) = T(:,:,10);
%   t(:,:,261) = T(:,:,11);
%   t(:,:,287) = T(:,:,12);
%   
%   t(:,:,313) = T(:,:,13);
%   t(:,:,339) = T(:,:,14);
%   t(:,:,365) = T(:,:,15);
%   t(:,:,391) = T(:,:,16);
%   t(:,:,417) = T(:,:,17);
%   t(:,:,443) = T(:,:,18);
%   t(:,:,469) = T(:,:,19);
%   t(:,:,495) = T(:,:,20);
%   t(:,:,521) = T(:,:,21);
%   t(:,:,547) = T(:,:,22);
%   t(:,:,573) = T(:,:,23);
%   t(:,:,599) = T(:,:,24);
%   
%   t(:,:,625) = T(:,:,25);
%   t(:,:,651) = T(:,:,26);
%   t(:,:,677) = T(:,:,27);
%   t(:,:,703) = T(:,:,28);
%   t(:,:,729) = T(:,:,29);
%   t(:,:,755) = T(:,:,30);
%   t(:,:,781) = T(:,:,31);
%   t(:,:,807) = T(:,:,32);
%   t(:,:,833) = T(:,:,33);
%   t(:,:,859) = T(:,:,34);
%   t(:,:,885) = T(:,:,35);
%   t(:,:,911) = T(:,:,36);
%   
%   t(:,:,937) = T(:,:,37);
%   t(:,:,963) = T(:,:,38);
%   t(:,:,989) = T(:,:,39);
%   t(:,:,1015) = T(:,:,40);
%   t(:,:,1041) = T(:,:,41);
%   t(:,:,1067) = T(:,:,42);
%   t(:,:,1093) = T(:,:,43);
%   t(:,:,1119) = T(:,:,44);
%   t(:,:,1145) = T(:,:,45);
%   t(:,:,1171) = T(:,:,46);
%   t(:,:,1197) = T(:,:,47);
%   t(:,:,1223) = T(:,:,48);
%   
%   t(:,:,1249) = T(:,:,49);
%   t(:,:,1275) = T(:,:,50);
%   t(:,:,1301) = T(:,:,51);
%   t(:,:,1321) = T(:,:,52);
%   t(:,:,1341) = T(:,:,53);
%   t(:,:,1361) = T(:,:,54);
%   t(:,:,1381) = T(:,:,55);

t(:,:,1) = t(:,:,1+1);
t(:,:,27) = t(:,:,27+1);
t(:,:,53) = t(:,:,53+1);
t(:,:,79) = t(:,:,79+1);
t(:,:,105) = t(:,:,105+1); 
t(:,:,131) = t(:,:,131+1); 
t(:,:,157) = t(:,:,157+1); 
t(:,:,183) = t(:,:,183+1); 
t(:,:,209) = t(:,:,209+1); 
t(:,:,235) = t(:,:,235+1); 
t(:,:,261) = t(:,:,261+1); 
t(:,:,287) = t(:,:,287+1); 

t(:,:,313) = t(:,:,313+1); 
t(:,:,339) = t(:,:,339+1); 
t(:,:,365) = t(:,:,365+1); 
t(:,:,391) = t(:,:,391+1); 
t(:,:,417) = t(:,:,417+1); 
t(:,:,443) = t(:,:,443+1); 
t(:,:,469) = t(:,:,469+1); 
t(:,:,495) = t(:,:,495+1); 
t(:,:,521) = t(:,:,521+1); 
t(:,:,547) = t(:,:,547+1); 
t(:,:,573) = t(:,:,573+1); 
t(:,:,599) = t(:,:,599+1); 

t(:,:,625) = t(:,:,625+1); 
t(:,:,651) = t(:,:,651+1); 
t(:,:,677) = t(:,:,677+1); 
t(:,:,703) = t(:,:,703+1); 
t(:,:,729) = t(:,:,729+1); 
t(:,:,755) = t(:,:,755+1); 
t(:,:,781) = t(:,:,781+1); 
t(:,:,807) = t(:,:,807+1); 
t(:,:,833) = t(:,:,833+1); 
t(:,:,859) = t(:,:,859+1); 
t(:,:,885) = t(:,:,885+1); 
t(:,:,911) = t(:,:,911+1); 

t(:,:,937) = t(:,:,937+1); 
t(:,:,963) = t(:,:,963+1); 
t(:,:,989) = t(:,:,989+1); 
t(:,:,1015) =t(:,:,1015+1);
t(:,:,1041) =t(:,:,1041+1);
t(:,:,1067) =t(:,:,1067+1);
t(:,:,1093) =t(:,:,1093+1);
t(:,:,1119) =t(:,:,1119+1);
t(:,:,1145) =t(:,:,1145+1);
t(:,:,1171) =t(:,:,1171+1);
t(:,:,1197) =t(:,:,1197+1);
t(:,:,1223) =t(:,:,1223+1);

t(:,:,1249) =t(:,:,1249+1);
t(:,:,1275) =t(:,:,1275+1);
t(:,:,1301) =t(:,:,1301+1);
t(:,:,1321) =t(:,:,1321+1);
t(:,:,1341) =t(:,:,1341+1);
t(:,:,1361) =t(:,:,1361+1);
t(:,:,1381) =t(:,:,1381+1);
  t_min = zeros(1601,ntime);
 
  for i = 1:1601
      for j = 1:1000
          t_min(i,j) = min(t(:,j,i));
      end
  end
  
  temp = zeros(nodes,ntime); 
%--------------------------------------------------------------------------
  ndarray = zeros(3,161);

for j = 1:size(ndarray,2);
    ndarray(1,j) = 2249+j-1;
end
for i = 1:size(ndarray,1);
    ndarray(i,:)  = ndarray(1,:)+281*(i-1);
end

ndarray = ndarray';
ndarray = ndarray(:);
%--------------------------------------------------------------------------
  
for it = 1:ntime
    
        for i = 1:161
            temp(ndarray(323+i-1),it) = t_min(1+10*(i-1),it);
        end
   
end


% for it = 1:ntime
%     
%         for i = 1:161
%             temp(ndarray(323+i-1),it) = 0;
%         end
%    
% end
