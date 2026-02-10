clc;clear;close all

wH = 3; 
phi = 0.32; 
t_max = 6*24*60*60; 
L = 100;
G = 0.1;
dx = 1e-4; 


phpa_set = [2, 4, 6,8,10,12]; 
kh_set = logspace(-7, -2, 30); 
wH_set = [3];
for phpa = phpa_set
    for kh = kh_set
        for wH = wH_set
            P0 = 1000*9.8*wH; 
            gradP_p0= -1;
            gradP = P0/L; 
            % rheology
            mu_w = 0.001; 
            load(sprintf('PHPA%02d_para.mat', phpa));
            % geo
            k = kh*(1e-3/1000/9.8); 
            R_eff = 0.18*sqrt(k)/(phi^2.3*(1-phi)^2.1);

            counter = 1;
            t = 0;
            ui = zeros(1, 10);
            Sat = zeros(1, 10);
            gradP_p = zeros(1, 10);
            mu_eff = zeros(1, 10);
            stat = zeros(1, 10);

            while t<t_max
                Sat(counter) = counter*dx/L;
                [gradP_p(counter), mu_eff(counter), stat(counter)] = calculateGradPp(mu_w, mu_0, mu_inf, m, n, G, R_eff, gradP, Sat(counter), gradP_p0);
                ui(counter) = gradP_p(counter)*k/mu_eff(counter)/phi; %interstitial velocity
                dt = dx/ui(counter);
                t = t+dt;
                counter=counter+1;
            end
            dt = dx./(0.5*(ui(1:end-1)+ui(2:end)))/60/60/24; %second->day
            t = [0, cumsum(dt)];
            xp = [0, 0.5*(Sat(1:end-1)+Sat(2:end))*L];
            slope = diff(log(xp))./diff(log(t));
            save(sprintf('./Data/Trench_PHPA%02d_k%1.1e_wH%1.1f.mat', phpa, kh, wH), 't', 'xp', 'ui', 'gradP_p','slope', 'stat');
        end
    end
end

