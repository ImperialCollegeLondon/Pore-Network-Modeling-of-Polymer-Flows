classdef PowerLawModel < handle
    %% rheology model Q = g(gradP, G, Rt)
    % power law model with two bounds
    properties
        mu_0, mu_inf, m, n; % model parameters
        gamma_0, gamma_inf;
    end

    methods
        function obj = PowerLawModel(mu_0, mu_inf, m, n)
            obj.mu_0 = mu_0;
            obj.mu_inf = mu_inf;
            obj.m = m;
            obj.n = n;
            obj.gamma_0 = (mu_0/m)^(1/(n-1));
            obj.gamma_inf = (mu_inf/m)^(1/(n-1));
        end

        %% public methods
        function mu_app = calculatemu_app(obj, gamma)
            % calculate apparent mu
            if gamma<obj.gamma_0
                mu_app = obj.mu_0;
            elseif gamma<obj.gamma_inf
                mu_app = obj.m*gamma.^(obj.n-1);
            else
                mu_app = obj.mu_inf;
            end
        end

        function [q, g, stat] = calculateQ(obj, gradP, G, Rt)
            % calculate throat flow rate (q_ij)
            sgn = sign(gradP);
            gradP = abs(gradP);
            tau_w = Rt/2*gradP;
            if G>=1.0/4.0/pi            % circular
                a = 0.25; b = 0.75;
                A = pi*Rt*Rt; C = 2*pi*Rt; Dh = 2.0*Rt; 
            elseif G>sqrt(3)/36         % square
                a = 0.2121; b = 0.6766;
                A = 4*Rt*Rt; C = 8*Rt; Dh = 4*A/C; 
            else                        % isosceles triangle
                a = 12.8035*G^2 + 0.9881*G + 0.1103;
                b = -344.3033*G^3 + 26.7163*G^2 + -0.2891*G + 0.6363;
                Dh = 2*Rt;
                A = Rt*Rt/(4*G);
            end
            k = Dh^2/(32*(a+b));

            if abs(tau_w) < obj.mu_0*obj.gamma_0
                q = Dh^2*gradP*A/(32*(b + a)*obj.mu_0);
                g = Dh^2*A/(32*(b + a)*obj.mu_0);
                stat = 1;
            elseif abs(tau_w) <= obj.mu_inf*obj.gamma_inf
                q = -Dh*(4^(b/a)*(Dh*gradP).^(-b/a)*obj.mu_0^((b*obj.n + a)/(a*(-1 + obj.n)))*a*(-1 + obj.n)*obj.m^((-a - b)/(a*(-1 + obj.n))) - obj.n*(Dh*gradP).^(1/obj.n)*obj.m^(-1/obj.n)*4^(-1/obj.n)*(b + a))*A/(8*(b + a)*(b*obj.n + a));
                g = q./gradP;
                stat = 2;
            else
                q = -Dh*(4^(b/a)*(Dh*gradP).^(-b/a)*a*obj.mu_inf*(-1 + obj.n)*(obj.mu_0^((b*obj.n + a)/(a*(-1 + obj.n))) - obj.mu_inf^((b*obj.n + a)/(a*(-1 + obj.n))))*obj.m^((-a - b)/(a*(-1 + obj.n))) - Dh*gradP*(b*obj.n + a)/4)*A/(8*(b*obj.n + a)*(b + a)*obj.mu_inf);
                g = q./gradP;
                stat = 3;
            end
            q = q*sgn;
        end
    end
end