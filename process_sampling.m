clear;clc;close all;
delete(gcp('nocreate'));
% parpool(25) %for parallel computing

n_samples = 100; %Number of repeats
LP = linspace(1,15,15); %Testing differnt plate length values
K=0.25; %Spring constant value to be tried
LC = sqrt(300/K); %Length of tether corresponding to spring constant with N=100

%Initialize variables used for storing
avg_mat = cell(1,length(LP)); 
avg_tk = cell(1,length(LP));
asymp_val = zeros(length(LP),1);
avg_NA = zeros(length(LP),1);
for k = 1:length(LP)
    lpp = LP(k);
    msd_smpl_i = cell(1,n_samples);
    for i=1:n_samples %Repeat for different samples
    %     [xp,t]  = MC_nuc_pore_latest_periodic_func( plate length, length of chain,ka_plate,kd_plate,force-dependent parameter);
        [xp,t_k,av_na]  = MC_nuc_pore_latest_periodic_func(lpp, LC, 1, 10, 0); %Run a Monte Carlo time simulation for given parameters
        X = zeros(1,1,length(xp)); X(1,1,:)=xp; 
        [msd,dtime]=computeMSD(X, length(t_k), 0, 2); %Calculate mean squared distance(msd) of the plate
        msd_smpl_i{i}  = msd(:,1);
    end
    nstep = length(t_k);
    avg = zeros(nstep-1,1);
    msd_smpl = cell2mat(msd_smpl_i);
    for j=1:nstep-1
        avg(j) = mean(msd_smpl(j,:));
    end
    avg_mat{k} = avg;
    avg_tk{k} = t_k(1:nstep-1);
    avg_NA(k) = av_na;
end
save LP_high
