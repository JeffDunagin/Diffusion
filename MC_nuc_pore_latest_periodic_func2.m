function [ xp,t, avg_na ] = MC_nuc_pore_latest_periodic_func2( lp, Lc,kap,kdp,alpha1)
close all;
a = 1; %Spacing between links
% L=a*(Nlink-1);% Length of pore
% lp = 20*a; %length of plate/cargo
w = a; %Vertical distance between lower wall and plate.
% pt = 0.1*w; %Plate thickness for plotting
% Lc = 50*a; %Contour Length of the polymer chain
rm = 0*4*a; % Minimum distance for repulsion forces between plate and tethers

Dp = 0; %Plate Diffusivity
Dw = 250; %Wall attachments Diffusivity
mu = 1;%1/Dp; %Fluid drag coefficient 
epso = 1E-10; %Energy parameter of the WCA potential

N = 100; %Number of monomer units in the tether
nb2 = Lc^2/3/N; %Natural polymer length Nb^2/3- to be fine tuned
klink = 1/nb2;  %stiffness of link - 3/2Nb^2
kdw = 0*30;  % rate of detachment at zero force
kaw = 0*30;  % rate of attachment
% fatw0 = kaw/(kdw+kaw);  %fraction of attached links (starting from equilibrium)
fatw0 = 0;  %fraction of attached links (starting from equilibrium)
% kdp = 10;  % rate of detachment at zero force
% kap = 10;  % rate of attachment
fatp0 = kap/(kdp+kap);  %fraction of attached links (starting from equilibrium)
%fatp0 = 0;  %fraction of attached links (starting from equilibrium)
dt1 = 1/max([kdp,kap,kdw,kaw]); % length of time step- fix to include spring
dt2 = mu/klink; %Characteristic time scale for harmonic well
% dt3 = rm^2/(2*Dw);%Should it be a function of epsilon?
% dt_vect = [dt1 dt2 dt3];
dt_vect = [dt1 dt2];
dt = min(dt_vect)/500;
dt_vect(isnan(dt_vect))=0;
ttot = 20*max(dt_vect); %Long total 
nstep = round(ttot/dt); %number of time steps
% nstep = 1E5;
% nstep = 1500; %number of time steps
% alpha1 = 1; % Force sensitivity of link detachment from plate
alpha2 = 0*2; % Force sensitivity of link detachment from wall
f_max = 2*sqrt(2*Dp/dt)*mu;

%Initialization
Nat = zeros(nstep,2); % Number of attached links
% dF = zeros(nstep,1);  % Force increments
% Pa = zeros(nstep,1);  % Probabilty of attachment as a function of time
% Pd = zeros(nstep,1);  %Probability of detachment as a function of time
% fat = zeros(nstep,1); % Fraction of attached links
% fdt = zeros(nstep,1); %Fraction of detached links
% deltaU = zeros(nstep,1); %Displacement increments
% U = zeros(nstep,1); % Total sliding displacement at each time step (assumed equal for each link)
xp = zeros(nstep,1); % Position of plate
rep_range = 2*(sqrt(Lc^2-w^2)+rm)+lp; %Defining max length range for consideration
L = rep_range;

xp(1)=0;
% temp = xp(1)+ -L/2:a:L/2;
temp = -L/2 :a:L/2;
Nlink = length(temp);
matlink = zeros(Nlink,4); %Matrix of links
%matlink(:,1) - wall state ( attached = 1, detached = 0) 
%matlink(:,2) - plate state ( attached = 1, detached = 0) 
%matlink(:,3) - Position of link on the lower wall
%matlink(:,4) - Displacement of link
matlink(:,3) = temp';
pos_noise = 0*a/3;
rnum = normrnd(0,pos_noise,Nlink,1);
matlink(:,3)=matlink(:,3)+rnum;



for ilink = 1:Nlink % Matlink at time 0
   matlink(ilink,4) = normrnd(0,sqrt(nb2),1); %strain
   rn = rand;
   if abs(matlink(ilink,3)+matlink(ilink,4)-xp(1))<lp/2
       if rn < fatp0
          matlink(ilink,2) = 1;    % bond is attached to plate
       end
   else
       if rn < fatw0
          matlink(ilink,1) = 1;    % bond is attached to wall
       end 
   end
end
matlink(matlink(:,1)+matlink(:,2)<1,4)=0;

t=zeros(nstep,1);
Fa=zeros(nstep,1);
F_repuls=zeros(nstep,1);
stretch=zeros(Nlink,nstep);
t(1)=0;
% k1=0;
for it = 1:nstep-1
    
    tether_repls = zeros(Nlink,1);
%     tether_pos = zeros(Nlink,1);
    Fp = 0;
%     relv_link = find(abs(matlink(:,3)-xp(it)) < lp/2+rep_range);
    relv_link = 1:Nlink;
    %update life of attached bond to plate
%     for ilink = 1:Nlink
    for ilink = relv_link(1):relv_link(end)
       rn = rand;
       %Force before for detachment sensitivity
       tether_length = sqrt( (matlink(ilink,4))^2 + w^2);
       force_tether = klink*matlink(ilink,4)*(1 + 3/5 *(tether_length/Lc)^2 + 99/175 *(tether_length/Lc)^4 + 513/875 *(tether_length/Lc)^6); %From series expansion of inverse Langevin function
       kdp1 = kdp*exp(alpha1*force_tether); %Force dependent detachment
       kdw1 = kdw*exp(alpha2*force_tether); %Force dependent detachment
       if matlink(ilink,2)==1
           if rn < kdp1*dt
               matlink(ilink,2) = 0;   % bond is detached from plate
               matlink(ilink,4) = 0;
           end
       elseif matlink(ilink,1)==1
           if rn < kdw1*dt
               matlink(ilink,1) = 0;   % bond is detached from wall
               matlink(ilink,4) = 0;
           end
       else
           disp = normrnd(0,sqrt(nb2),1);
           if abs(matlink(ilink,3)+disp-xp(it))<lp/2
               if rn < kap*dt
                   matlink(ilink,2) = 1;    % bond is attached to plate
                   matlink(ilink,4) = disp;
               end
           else
               if rn < kaw*dt
                   matlink(ilink,1) = 1;    % bond is attached to wall
                   matlink(ilink,4) = disp;
               end
           end
       end
       
       %Force after detachment/attachment events
%        tether_length = sqrt( (matlink(ilink,4))^2 + w^2);
%        force_tether = klink*matlink(ilink,4)*(1 + 3/5 *(tether_length/Lc)^2 + 99/175 *(tether_length/Lc)^4 + 513/875 *(tether_length/Lc)^6); %From series expansion of inverse Langevin function
       force_tether = klink*matlink(ilink,4); %Linear Spring
       link_tip_pos = matlink(ilink,3)+matlink(ilink,4);
       F_bound =  matlink(ilink,2)*force_tether;
       
       %Plate-tether repulsion forces
       x_edge_front = matlink(ilink,1)*(link_tip_pos-(xp(it)+lp/2));
       x_edge_back = matlink(ilink,1)*(link_tip_pos-(xp(it)-lp/2));
       if x_edge_front>0 && abs(x_edge_front)< rm
           f_edge_front =  abs (epso*rm^6/abs(x_edge_front)^7 * (1 - rm^6/abs(x_edge_front)^6 ));%Magnitude of repulsive force from behind the plate
           if f_edge_front > f_max
               f_edge_front=f_max;
               display('error- large force');
           end
           tether_repls(ilink) = f_edge_front*dt/mu;
       elseif x_edge_back<0 && abs(x_edge_back)< rm
           f_edge_back =  abs(epso*rm^6/abs(x_edge_back)^7 * (1 - rm^6/abs(x_edge_back)^6 )); %Magnitude of repulsive force from behind the plate
           if f_edge_back > f_max
               f_edge_back=f_max;
               display('error- large force');
           end
           tether_repls(ilink) = -f_edge_back*dt/mu;
       else 
           f_edge_front = 0; f_edge_back = 0;
       end
       Fp = Fp - F_bound - f_edge_front + f_edge_back;
       
    end    
    
    F_repuls(it) = sum(tether_repls*mu/dt);
    Fa(it)=Fp; %Force from attached links in time
    stretch(:,it) = matlink(:,4);
    Nat(it,1) = sum(matlink(:,1)); %number attached to wall
    Nat(it,2) = sum(matlink(:,2));%number attached to plate
    xd = normrnd(0,sqrt(2*Dp*dt),1); % Diffusive displacement (white noise)
%    Fp=0;
    delx = Fp*dt/mu +xd;
    xp(it+1) = xp(it)+ delx; %Change in displacement of plate center
    
    matlink(:,4)=matlink(:,4)+matlink(:,2)*delx + matlink(:,1).*(normrnd(0,sqrt(2*Dw*dt),Nlink,1) + tether_repls); %Modified strains for plate and wall tethers
    
    %After change in plate position- Update matlink to include new tethers and delete old tethers
    if xp(it+1)+L/2 >= matlink(end,3)+2*pos_noise %Moves to the right too much
        nlinks_chnge = round( ((xp(it+1)+L/2) - (matlink(end,3)))/a ) ;
        matlink_tmp = zeros(nlinks_chnge,4);
        temp = matlink(end,3)+ (a:a:nlinks_chnge*a);
        rows = Nlink-nlinks_chnge+1:Nlink;
        matlink(1:Nlink-nlinks_chnge,:) = matlink(nlinks_chnge+1:Nlink,:);
        matlink_tmp(:,3) = temp';
        rnum = normrnd(0,pos_noise,abs(nlinks_chnge),1);
        matlink_tmp(:,3)=matlink_tmp(:,3)+rnum;
        for ilink=1:nlinks_chnge %Define new bonds
            matlink_tmp(ilink,4) = normrnd(0,sqrt(nb2),1); %strain
            rn = rand;
            if abs(matlink_tmp(ilink,3)+matlink_tmp(ilink,4)-xp(it+1))<lp/2
                if rn < fatp0
                    matlink_tmp(ilink,2) = 1;    % bond is attached to plate
                end
            else
                if rn < fatw0
                    matlink_tmp(ilink,1) = 1;    % bond is attached to wall
                end
            end
        end
        matlink_tmp(matlink_tmp(:,1)+matlink_tmp(:,2)<1,4)=0;
        matlink(rows,:) = matlink_tmp;
    elseif  xp(it+1)-L/2 <= matlink(1,3)-2*pos_noise %Moves to the left too much
        nlinks_chnge = round ( (matlink(1,3) - (xp(it+1)-L/2))/a ) ;
        matlink_tmp = zeros(nlinks_chnge,4);
        temp = matlink(1,3)- (a:a:nlinks_chnge*a);
        rows = 1:nlinks_chnge;
        matlink(nlinks_chnge+1:Nlink,:) = matlink(1:Nlink-nlinks_chnge,:);
        matlink_tmp(:,3) = temp';
        rnum = normrnd(0,pos_noise,abs(nlinks_chnge),1);
        matlink_tmp(:,3)=matlink_tmp(:,3)+rnum;
        for ilink=1:nlinks_chnge %Define new bonds
            matlink_tmp(ilink,4) = normrnd(0,sqrt(nb2),1); %strain
            rn = rand;
            if abs(matlink_tmp(ilink,3)+matlink_tmp(ilink,4)-xp(it+1))<lp/2
                if rn < fatp0
                    matlink_tmp(ilink,2) = 1;    % bond is attached to plate
                end
            else
                if rn < fatw0
                    matlink_tmp(ilink,1) = 1;    % bond is attached to wall
                end
            end
        end
        matlink_tmp(matlink_tmp(:,1)+matlink_tmp(:,2)<1,4)=0;
        matlink(rows,:) = matlink_tmp;
    end
    
    %check=sum(matlink(:,1).*matlink(:,2));
    t(it+1)= t(it)+dt;
% %% Animation of moving plate through a nuclear pore made of a nuclear brush
%     if mod(it,100) ==0
%         k1=k1+1;
%         multip_plate = 1;
%         multip_tether = 1;
%          h = figure;
%     %      end_link_pos = find(abs(matlink(:,3)- xp(it+1))< lp/2+5*a);
% %         end_link_pos = find(abs(matlink(:,3)- xp(1))< lp/2+5*a);
%         end_link_pos = relv_link;
%          for i=1: length(end_link_pos)
%              y=end_link_pos(i);
%              hold on; 
%     %          h1(i)=plot(matlink(y,3)+[0 5*matlink(y,4)],[0 w],'color','b','linewidth',1);
%              if matlink(y,2)== 1 %attached to plate
% %                  h1(i)=plot(matlink(y,3)+[0 multip_tether *matlink(y,4)],[0 w],'color','g','linewidth',1);
% %                  h2(i)=plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[2*w+pt w+pt],'color','g','linewidth',1);
%                  plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[0 w],'color','g','linewidth',1);
%                  plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[2*w+pt w+pt],'color','g','linewidth',1);
%              elseif matlink(y,1)== 1 %attached to wall
% %                  h1(i)=plot(matlink(y,3)+[0 multip_tether *matlink(y,4)],[0 w],'color','b','linewidth',1);
% %                  h2(i) =plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[2*w+pt w+pt],'color','b','linewidth',1);
%                    plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[0 w],'color','b','linewidth',1);
%                    plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[2*w+pt w+pt],'color','b','linewidth',1);
%              else % detached
% %                  h1(i)=plot(matlink(y,3)+[0 multip_tether *matlink(y,4)],[0 w],'color','r','linewidth',1);
% %                  h2(i)=plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[2*w+pt w+pt],'color','r','linewidth',1);
%                  plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[0 w],'color','r','linewidth',1);
%                  plot(matlink(y,3)+[0 multip_tether*matlink(y,4)],[2*w+pt w+pt],'color','r','linewidth',1);
%              end
%          end
%          xprofile = [multip_plate *xp(it+1)-lp/2 multip_plate*xp(it+1)+lp/2 multip_plate*xp(it+1)+lp/2 multip_plate*xp(it+1)-lp/2 multip_plate*xp(it+1)-lp/2];
%          yprofile = [w w w+pt w+pt w];
%          hold on; 
% %          h3=plot(xprofile,yprofile,'color','k','linewidth',2);
%          plot(xprofile,yprofile,'color','k','linewidth',2);
%          xlim([xp(1)-lp/2-1*a xp(1)+lp/2+1*a])
%          ylim([0 w+pt+w+pt]);
% %          pause(0.1)
% %          delete (h1);
% %          delete (h2);
% %          delete (h3);
%          M(k1)=getframe(h);
%          close all;
%     end
end
% movie2gif(M, 'plate_move_alpha_00_test.gif', 'LoopCount', 0, 'DelayTime', 0)
% % Plots w.r.t time
figure;
hold on
subplot(3,1,1) %Plate position
plot(t(1:it),xp(1:it)); %ylim([ -0.2 0.2]);
hold on
subplot(3,1,2)%Number of attached tethers to plate
plot(t(1:it),Nat(1:it,2))
% hold on
% subplot(5,1,3) %Number of attached tethers to wall
% plot(t(1:it),Nat(1:it,1))
hold on
subplot(3,1,3) %Mechanical Force on plate from tethers
plot(t(1:it),Fa(1:it,1))
hold on
% subplot(5,1,5) %Repsulsion Force on plate from tethers
% plot(t(1:it),F_repuls(1:it,1))

beg = round(nstep/10); %number 10 chosend for equilibrium
avg_na = mean(Nat(beg:end,2));

end

