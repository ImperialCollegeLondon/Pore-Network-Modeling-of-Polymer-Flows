classdef PoreNetwork < handle
    %% pore networks: including pores and throats

    properties
        L_stream; 

        throatConn;
        throatProp;
        nodeProp;

        CModel; 

        nPore; 
        nThroat; 

        P_vec;
        Q_vec; 
        q_vec; 
        k_vec;
        dP_vec;
        stat_vec; 

        K; 
    end

    properties(SetAccess=protected)
        NumIter_total; 
        ResError_plot;
        IncError_plot;

    end

    methods
        function obj = PoreNetwork(L, throatConn_,throatProp_, nodeProp_, CModel_)
            obj.L_stream = L;
            obj.nodeProp = nodeProp_;
            obj.throatConn = throatConn_;
            obj.throatProp = throatProp_;
            obj.nThroat = height(obj.throatConn);
            obj.nPore = height(obj.nodeProp);
            obj.P_vec = zeros(obj.nPore+2, 1, 'double');
            obj.Q_vec = zeros(obj.nPore+2, 1, 'double');
            obj.q_vec = zeros(obj.nThroat, 1, 'double');
            obj.k_vec = zeros(obj.nThroat, 1, 'double');
            obj.dP_vec = zeros(obj.nThroat, 1, 'double');
            obj.stat_vec = zeros(obj.nThroat, 1, 'double');

            obj.K = sparse(obj.nPore+2, obj.nPore+2);

            obj.CModel = CModel_;

        end

        %% potected methods
        function assemblyConductanceMatrix(obj)
            obj.K = sparse(obj.nPore+2, obj.nPore+2);
            for i = 1:obj.nThroat
                pore_i = obj.throatConn.poreI(i)+2;
                pore_j = obj.throatConn.poreJ(i)+2;
                obj.K(pore_i, pore_i) = obj.K(pore_i, pore_i) + obj.k_vec(i);
                obj.K(pore_i, pore_j) = obj.K(pore_i, pore_j) - obj.k_vec(i);
                obj.K(pore_j, pore_j) = obj.K(pore_j, pore_j) + obj.k_vec(i);
                obj.K(pore_j, pore_i) = obj.K(pore_j, pore_i) - obj.k_vec(i);
            end

        end

        function updateThroats(obj)
            obj.dP_vec = (obj.P_vec(obj.throatConn.poreI+2) - obj.P_vec(obj.throatConn.poreJ+2));
            gradP = obj.dP_vec./obj.throatConn.lenTot;
            for i = 1:obj.nThroat
                [obj.q_vec(i), g, obj.stat_vec(i)] = obj.CModel.calculateQ(gradP(i), obj.throatConn.shapeFact(i), obj.throatConn.radius(i));
                obj.k_vec(i) = g/obj.throatConn.lenTot(i);
            end
        end

        function calculateUnbalancedFlowRateAtPore(obj)
            obj.Q_vec = zeros(obj.nPore+2, 1, 'double');
            for i = 1:obj.nThroat
                pore_i = obj.throatConn.poreI(i)+2;
                pore_j = obj.throatConn.poreJ(i)+2;
                obj.Q_vec(pore_i) = obj.Q_vec(pore_i) + obj.q_vec(i);
                obj.Q_vec(pore_j) = obj.Q_vec(pore_j) - obj.q_vec(i);
            end
        end

        %% public methods
        function [Q, dP, ratio_stat2] = SolvePQcurve(obj, dP_min, dP_max, num, maxIter, tol, outPath)
            dP = linspace(dP_min, dP_max, num)';
            ratio_stat2 = zeros(2, 1);
            ratio_stat3 = zeros(2, 1);
            Q = zeros(size(dP)); Q(1) = obj.LinearSolver(dP_min);
            obj.OutputNetwork(outPath, 1);

            figure(1); grid on
            C = colororder;
            obj.ResError_plot = animatedline(Color=C(2,:));
            obj.IncError_plot = animatedline(Color=C(3,:));
            xlabel("iteration")
            ylabel("err")

            obj.NumIter_total = 0;
            for i = 2:length(dP)
                [Q(i), iter] = obj.IterativeSolver(dP(i), maxIter, tol);
                obj.NumIter_total = obj.NumIter_total+iter;
                ratio_stat2(i) = sum(obj.stat_vec==2)/obj.nThroat;
                ratio_stat3(i) = sum(obj.stat_vec==3)/obj.nThroat;
                
                obj.OutputNetwork(outPath, i);
            end
        end

        function Q = LinearSolver(obj, dP)
            Res = ones(obj.nPore+2, 1, 'double'); 
            Q_target = zeros(obj.nPore+2, 1, 'double');
            obj.updateThroats();

            obj.assemblyConductanceMatrix();
            K_ = obj.K;
            Q_target(1) = obj.K(1,1)*1e8*dP;
            Q_target(2) = 0.0;
            K_(1,1) = obj.K(1,1)*1e8;
            K_(2,2) = obj.K(2,2)*1e8;

            obj.P_vec = K_\Q_target;
            obj.updateThroats();
            obj.calculateUnbalancedFlowRateAtPore();
            Q = obj.K*obj.P_vec;
            Q = Q(1);
        end

        function [Q, iter] = IterativeSolver(obj, dP, maxIter, tol)
            Res = ones(obj.nPore+2, 1, 'double'); 
            Q_target = zeros(obj.nPore+2, 1, 'double');
            obj.P_vec(1) = dP;
            obj.updateThroats();
            obj.calculateUnbalancedFlowRateAtPore();
            Res = Q_target - obj.Q_vec; Res(1) = 0; Res(2) = (0);
            Q = obj.Q_vec(1);
            % iterating
            iter = 0;
            while iter<maxIter
                err_res = max(abs(Res(3:end)/Q));
                addpoints(obj.ResError_plot,obj.NumIter_total+iter,err_res);
                drawnow
                if err_res > tol
                    obj.assemblyConductanceMatrix();
                    obj.K(1,2:end) = 0; obj.K(2:end,1) = 0;
                    obj.K(2,3:end) = 0; obj.K(3:end,2) = 0;
                    inc = obj.K\Res;
                    err_inc = max(abs(inc(3:end)./dP));
                    addpoints(obj.IncError_plot,obj.NumIter_total+iter,err_inc);
                    drawnow
                    if err_inc < tol
                        break;
                    end
                    lamda = 0.5;
                    obj.P_vec(3:end) = obj.P_vec(3:end) + lamda*inc(3:end);

                    obj.updateThroats();
                    obj.calculateUnbalancedFlowRateAtPore();
                    Res = Q_target - obj.Q_vec; Res(1) = 0; Res(2) = (0);
                    iter = iter+1;
                    Q = obj.Q_vec(1);
                else
                    break;
                end

            end
        end

        function OutputNetwork(obj, Path, t)
            poreFileName = sprintf('%spore_%i.csv', Path, t);
            dlmwrite(poreFileName, obj.P_vec, ',');
            throatFileName = sprintf('%sthroat_%i.csv', Path, t);
            dlmwrite(throatFileName, [obj.q_vec, obj.dP_vec, obj.k_vec, obj.stat_vec], ',');
        end

    end
end