clc;clear;close all


wH = 3; 
phi = 0.32; 
P0 = 1000*9.8*wH; 
t_max = 6*24*60*60; 
L = 100;
G = 0.1;
dx = 1e-4; 
r0 = 0.5; 
rw = -1;
gradP_p0 = -1;
gradP = P0/L; 

phpa_set = [2, 4, 6, 8, 10, 12]; 
kh_set = logspace(-7, -2, 30); 
for phpa = phpa_set
    for kh = kh_set
        
        mu_w = 0.001; 
        load(sprintf('PHPA%02d_para.mat', phpa));
        
        k = kh*(1e-3/1000/9.8); % permeability m2
        R_eff = 0.18*sqrt(k)/(phi^2.3*(1-phi)^2.1);

        
        gamma_inf = (mu_inf/m)^(1/(n-1));
        gamma_0 = (mu_0/m)^(1/(n-1));
        gradP = linspace(0.1*2*mu_0*gamma_0/R_eff, 10*2*mu_inf*gamma_inf/R_eff, 10000);
        for i = 1:length(gradP)
            [mu_eff_(i), u_(i), ~, ~, ~] = calculateEffectiveViscosity(mu_0, mu_inf, m, n, G, R_eff, gradP(i));
        end

        counter = 1;
        t = 0;
        C0 = 10;
        k_ = k/phi; 
        ui = zeros(1, 10);
        rp = zeros(1, 10);
        C = zeros(1, 10);
        stat = zeros(1, 10);
        while t<t_max
            rp(counter) = r0+counter*dx;
            C(counter) = calculateGradPp_cylinder_upper(mu_w, mu_0, mu_inf, m, n, G, R_eff, mu_eff_, u_, k_, phi, P0, r0, rp(counter), rw, C0);
            C0 = C(counter);
            ui(counter) = C(counter)/rp(counter);
            [~, ~, ~, stat(counter)] = calculateEffectiveViscosity_u(mu_0, mu_inf, m, n, G, R_eff, mu_eff_, u_, ui(counter));
            dt = dx/ui(counter);
            t = t+dt;
            counter=counter+1;
        end
        dt = dx./(0.5*(ui(1:end-1)+ui(2:end)))/60/60/24;
        t = [0, cumsum(dt)];
        xp = [r0, 0.5*(rp(1:end-1)+rp(2:end))]-r0;
        slope = diff(log(xp))./diff(log(t));
        save(sprintf('./Data/BoreHole_PHPA%02d_k%1.1e.mat', phpa, kh), 't', 'xp', 'ui','slope', 'stat');
        fprintf('./Data/BoreHole_PHPA%02d_k%1.1e.mat\n', phpa, kh)
    end
end

