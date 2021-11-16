% This MATLAB class was developed to generate simulation results to "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads". Tomas Ramirez, Carlos Mosquera, Nader Alagha, arXiv, 2021
% License: This code is licensed under the GPLv3 license. If you in any way use this code for research that results in publications, please cite our  work as described above.
classdef ResourceAssignmentOneDimesion
    % This class enables the simulation of the different resource assignment strategies for a one-dimensional satellite scenario with
    % a two-colour scheme. The following strategies are simulated, with a two-step optimization strategy:
    % - POW: Flexible power allocation. Fixed user assignment to the beams
    % - BW: Flexible bandwidth allocation. Fixed user assignment to the beams
    % - BW-POW: Flexible power and bandwidth allocation. Fixed user assignment to the beams
    % - MAP: Fixed resources per beam. Flexible beam user assignment.
    % - BW-MAP: Flexible bandwidth allocation. Flexible beam user assignment.
    
    properties
        K % Number of beams
        M % Number of carriers per color
        N % Number of users
        Delta_W; % Carrier Bandwidth ( normalized to the total bandwidth)
        
        Nuser_beam %  Vector with the number of user per beam
        user_beams % Vector that indicates the dominant beam for each user
        users_b_index % Cell with the user indexes for each beam
        Gamma_wo_P  % Matrix with the channel values
        snr_car % Carrier SNR for uniform carrier power allocation.
        Max_P_sat % Overall power constraint
        Ref_Pb_Sat % Value for uniform power allocation
        Max_Pamp % Maximum amplifier power
        
        
        Req_beam  % Vector with the user requested traffic
        Req_user % Vector with the beam requested traffic
        
        Max_time_opt_second % Auxiliar variable for the second step process. Sets the maximum time for the optimization.
        effective_snr_beam % Auxiliar vector for the beam resource assignment. Effective SNR per beam
        tol % Auxiliar variable. Tolerance value
        
        SE_av  % Auxilixar variable.  log2(1+ SNR_eff_beam) , with SNR_eff_beam  the average SNR for a uniform distribution of user within a beam (  approx 13.5 dB )
        User_per_car  % Auxiliar variable.Estimation of the number of user which can be served by one carrier
        ub_B % Auxiliar variable.Upper bound for the bandwidth allocation in the genetic algorithm
        aux_Cn % Auxiliar variable for the genetic algorithm. Matrix with SNR values.
        
    end
    
     
    
    
    methods
        function obj = ResourceAssignmentOneDimesion(arg)
            % Construct function. Load the scenario data from the argument into the class object
            obj.K=arg.K;
            obj.M= arg.M;
            obj.N=arg.N;
            obj.Delta_W=arg.Delta_W;
            
            obj.Nuser_beam=  arg.Nuser_beam;
            obj.user_beams= arg.user_beams;
            obj.users_b_index= arg.users_b_index;
            obj.Gamma_wo_P=  arg.Gamma_wo_P;
            obj.snr_car=  arg.snr_car;
            obj.Max_P_sat=  arg.Max_P_sat;
            obj.Ref_Pb_Sat= arg.Ref_Pb_Sat;
            obj.Max_Pamp= arg.Max_Pamp;
            
            obj.Req_user=arg.Req_user;
            obj.Req_beam = arg.Req_beam;
            obj.Max_time_opt_second=arg.Max_time_opt_second;
            obj.tol=1e-5;
            
            obj.effective_snr_beam= ComputEffectiveSNR(obj);
            
            obj.SE_av=arg.SE_av;
            obj.User_per_car= round(1./(obj.Req_user(1)/(obj.SE_av/(2*obj.M))));
            obj.ub_B =   min(ceil(obj.Nuser_beam/obj.User_per_car),2*obj.M) ;
            obj.aux_Cn=arg.aux_Cn;
            
        end
        
        
        
        
        
        function [R_off,M_beam] = FlexibleBandwidth(obj)
            % The function implementes the two-step optimization process for flexible bandwidth and fixed beam user mapping
            % as it is described in "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021
            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       M_beam: Vector with the number of allocated carriers per beam
            
            
            %%%%%%%%%%%% FIRST STEP: Beam resource allocation
            
            % Define convex optimization problem with CVX
            cvx_solver default
            cvx_begin quiet
            variable w_beam(obj.K)
            expression R_off_beam(obj.K)
            expression obj_fun
            
            % Constraint to avoid the use of the same portion of spectrum in adjacent beams
            C_OverB_s=zeros(obj.K-1,obj.K);
            v_OverB_s= ones(obj.K-1,1);
            for aux_C=1:obj.K-1
                C_OverB_s(aux_C,:)= [  zeros(1, aux_C-1)  1  1 zeros(1,obj.K-aux_C-1)];
            end
            
            % Offered traffic per beam
            for index=1:obj.K
                if obj.effective_snr_beam(index)~=0
                    R_off_beam(index)=w_beam(index)* log( 1+ obj.effective_snr_beam(index)/obj.K   )/log(2);
                else
                    R_off_beam(index)=0;
                end
            end
            
            % Objective function, Quadratic unment demand (QU)
            obj_fun=  (obj.Req_beam.'-R_off_beam).^2; %QU
            % Minimization problem
            minimize( sum(obj_fun) );
            
            % Apply optimization constraints
            subject to
            w_beam>=0
            C_OverB_s*w_beam<=v_OverB_s
            cvx_end
            
            
            % Obtain discrete carrier allocation for each beam from the continuos bandwidth from the first step( w_beam).
            M_beam=BandwidthDisc(obj,w_beam,1);
            
            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(size(obj.Req_user)); % Variable with the offerd traffic to each user
            
            for ind_beam=1:obj.K   % Beam selection
                N_user_beam=obj.Nuser_beam(ind_beam); % Number of users  in the selected beam
                M_b=M_beam(ind_beam); % Number of carriers in the selected beam
                users_index=obj.users_b_index{ind_beam}; % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end
            
            
        end
        
        
        
        
        function [R_off,P_flex] = FlexiblePower(obj)
            % The function implementes the two-step optimization process for flexible power and fixed beam user mapping
            % as it is described in "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021
            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       P_flex: Vector with the allocated power per beam.
            
            
            %%%%%%%%%%%% FIRST STEP: Beam resource allocation
            
            % Optimization of power per HPA
            Nvar=obj.K/2; % One amplifier per pair of beams.
            p0=zeros(1,Nvar);
            ub=obj.Max_Pamp/obj.Max_P_sat*ones(1,Nvar); % Upper bound for the power allocation
            lb=zeros(1,Nvar); % Lower bound for the power allocation
            
            % Settings for the optimization
            alg='sqp';
            options_FlexP_Beam = optimoptions('fmincon','OptimalityTolerance',1e-10,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',4e5,'Display','off','Diagnostics','off','Algorithm',alg);
            A=ones(1,Nvar);
            b=1;
            
            % Optimization ( with fmincon)
            [P_flex_aux,~,~,~] =  fmincon( @(p) Flex_Power_Beam(p,obj.Req_beam,obj.K,obj.effective_snr_beam),p0,A,b,[],[],lb,ub,[],options_FlexP_Beam);
            
            % Transalte power per HPA to power per beam
            P_flex=obj.Max_P_sat*repelem(P_flex_aux,2)/2; % Power per beam ( Same power for each beam sharing the HPA)
            
            
            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(size(obj.Req_user)); % Variable with the offerd traffic to each user
            
            for index_beam=1:obj.K % Index of the selected beam
                N_user_beam=obj.Nuser_beam(index_beam); % Number of users  in the selected beam
                users_index=obj.users_b_index{index_beam};  % Auxiliar variable
                M_b=obj.M; % Fixed number of carriers per beam
                Req_user_aux=obj.Req_user(users_index).'; % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,index_beam,users_index,Req_user_aux, P_flex(index_beam), 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end
            
            
        end
        
        
        
        
        
        
        function [R_off,Nuser_beam_UM_FixedRes]=FixResFlexibleMapping(obj)
            % The function implementes the two-step optimization process for fixed beam resources and flexible beam user mapping
            % as it is described in "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021
            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       Nuser_beam_UM_FixedRes: Vector that contains the association between users and beams.
            
            
            
            % Relaxed problem to obtain the beam-user assignment
            cvx_solver default
            cvx_begin quiet
            variable w_beam_user(obj.K,obj.N)
            expression Roff_sim2(obj.K,obj.N)
            expression obj_fun
            for index_beam=1:obj.K
                for index_user=1:obj.N
                    Roff_sim2( index_beam,index_user)=w_beam_user(index_beam,index_user)*obj.Delta_W* log(1 +  obj.snr_car(index_user,index_beam)   )/log(2);
                end
            end
            obj_fun= ( obj.Req_user.' - sum(Roff_sim2,1)).^2;
            minimize( sum(obj_fun) );
            subject to
            0<=w_beam_user<=1
            sum(w_beam_user,1)<= 1 %  Constraint of one carrier per user
            sum(w_beam_user,2)<=obj.M % Constraint of fixed bandwidth ( M carriers) per beam
            cvx_end
            
            % Number of user per beam
            Assig_UM_FlexB=zeros(1,obj.N);
            for ind_n=1:obj.N
                [~,indmax]=max( Roff_sim2(:,ind_n) );
                Assig_UM_FlexB(ind_n)= indmax;
            end
            
            % Obtain association between users and beams.
            Nuser_beam_UM_FixedRes=hist(Assig_UM_FlexB,1:obj.K);
            
            
            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            
            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user
            
            for ind_beam=1:obj.K % Index of the selected beam
                N_user_beam=Nuser_beam_UM_FixedRes(ind_beam); % Number of user in the selected beam
                M_b=obj.M; % Fixed carrier allocation per beam
                users_index=find(Assig_UM_FlexB==ind_beam); % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end
            
            
        end
        
        
        function  [R_off,Nuser_beam_UM_FlexB,M_beam]=FlexBandwidthFlexMapping(obj)
            % The function implementes the two-step optimization process for flexible bandwidth and beam user mapping
            % as it is described in "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021
            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       Nuser_beam_UM_FlexB: Vector that contains the association between users and beams
            %       M_beam: Vector with the number of allocated carriers per beam
            
            
            % Relaxed problem to obtain the beam-user assignment and bandwidth allocation
            C_OverB_sim=zeros(obj.K-1,obj.K);
            v_OverB_sim= 2*obj.M*ones(obj.K-1,1);
            for aux_C=1:obj.K-1
                C_OverB_sim(aux_C,:)= [  zeros(1,aux_C-1)  ones(1,2 )   zeros( 1,obj.K-aux_C-1)];
                
            end
            cvx_solver default
            cvx_begin quiet
            variable w_sim(obj.K,obj.N)
            expression Roff_sim(obj.K,obj.N)
            expression obj_fun
            for index_beam=1:obj.K
                for index_user=1:obj.N
                    
                    Roff_sim( index_beam,index_user)=w_sim(index_beam,index_user)*obj.Delta_W* log(1 +  obj.snr_car(index_user,index_beam)   )/log(2);
                end
            end
            obj_fun= ( obj.Req_user.' - sum(Roff_sim,1)).^2;
            minimize( sum(obj_fun) );
            subject to
            0<=w_sim<=1
            sum(w_sim,1)<= 1 %  Constraint of one carrier per user
            C_OverB_sim*sum(w_sim,2)<=v_OverB_sim % Constraint to avoid the reuse of bandwith in adjacent beams
            cvx_end
            
            
            % Extract beam user mapping from the results of the relaxed problem. If no rates are assigned to a user, the user is assigned to its dominant beam
            Assig_UM_FixedRes=zeros(1,obj.N);
            for ind_n=1:obj.N
                
                [vmax,indmax]=max( Roff_sim(:,ind_n) );
                if vmax<obj.tol
                    Assig_UM_FixedRes(ind_n)=obj.user_beams(ind_n);
                else
                    Assig_UM_FixedRes(ind_n)= indmax;
                end
            end
            
            % Obtain discrete carrier allocation for each beam from the continuos bandwidth from the first step( w_beam).
            w_beam=sum(w_sim,2)+0;
            M_beam=BandwidthDisc(obj,w_beam,2);
            
            % Obtain beam-user mapping
            Nuser_beam_UM_FlexB=hist(Assig_UM_FixedRes,1:obj.K);
            
            
            
            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user
            
            for ind_beam=1:obj.K % Index of the selected beam
                N_user_beam=Nuser_beam_UM_FlexB(ind_beam); % Number of users  in the selected beam
                M_b=M_beam(ind_beam); % Number of carriers  in the selected beam
                users_index=find(Assig_UM_FixedRes==ind_beam); % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % Requested uset traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end
        end
        
        
        
        function  [R_off,M_GA,P_GA]=FlexBandwidthPower(obj,P_TwoStep,M_beam,optionsGA)
            % The function implementes the two-step optimization process for flexible bandwidth and power allocation  as it
            % is described in "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021
            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       M_GA: Vector with the number of allocated carriers per beam
            %       M_beam: Vector with the allocated power per beam.
            
            
            
            %%%%%%%%%%%% FIRST STEP: Beam resource allocation
            
            
            % Genetic algorithm settings
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            aux_P=1:obj.K/2;  % Auxiliar variable
            aux_B=obj.K/2+1:obj.K+obj.K/2; % Auxiliar variable
            lb_P = zeros(1,obj.K/2); % Lower bound for the power allocation
            ub_P =obj.Max_P_sat/(obj.K*obj.M)*2*ones(1,obj.K/2);  % Upper bound for the power allocation
            lb_B = zeros(1,obj.K);  % Lower bound for the bandwidth allocation
            lb= [lb_P lb_B ];  % Lower bound for the genetic algorithm
            ub= [ub_P obj.ub_B]; % Upper bound for the genetic algorithm
            aux_pflex_car= P_TwoStep(1:2:obj.K)/obj.M; % Auxiliar variable. Power allocation for the fixed bandwidth and flexible power assignment.
            un_case=  [ obj.Max_P_sat/(obj.K*obj.M)*ones(1,obj.K/2) obj.M*ones(1,obj.K) ; % Initial individuals for the genetic algorithm. Solutions from other resource assigments
                aux_pflex_car obj.M*ones(1,obj.K);                                        % are employed as individuals.
                obj.Max_P_sat/(obj.K*obj.M)*ones(1,obj.K/2)   M_beam' ];
            options = optimoptions('ga','Display','none','InitialPopulationMatrix',un_case,'MaxGenerations',optionsGA.Max_Gen,'UseParallel', true,'CreationFcn','gacreationuniform','EliteCount',optionsGA.Elite_num,'MutationFcn', {@mutationadaptfeasible,optionsGA.pmut},'PopulationSize',optionsGA.PoP,'SelectionFcn', {@selectiontournament,5});
            
            
            % Genetic algorithm
            [x_ga,~,~,~,~,~] = ga(@(x) ga_function(x,obj),length(lb),A,b,Aeq,beq,lb,ub,[],4:9,options);
            
            
            %  Process the results from the GA from MATLAB ( Cheking of the values is needed due to the GA implementation in MATLAB)
            % Check Bandwidth
            [B] = Check_B(x_ga(aux_B), obj);
            % Check Power
            P=repelem(x_ga(aux_P),2).*B*2*obj.M;
            P_amp=zeros(1,obj.K/2);
            for index=1:obj.K/2
                P_amp(index)=sum(P( 2*(index-1)+1:2*index));
                
                if    P_amp(index)> obj.Max_Pamp
                    k=  obj.Max_Pamp/P_amp(index);
                    P( 2*(index-1)+1:2*index)= k*P( 2*(index-1)+1:2*index);
                end
            end
            if sum(P_amp)>obj.Max_P_sat
                k=  obj.Max_P_sat/sum(P_amp);
                P= k*P;
            end
            
            % Final resource assignment
            P_GA=P;
            M_GA=B*2*obj.M;
            
            
            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user
            for index_beam=1:obj.K
                N_user_beam=obj.Nuser_beam(index_beam); % Number of users  in the selected beam
                users_index=obj.users_b_index{index_beam}; % Auxiliar variable
                Req_user_aux=obj.Req_user(users_index).'; % Requested user traffic in the selected beam
                M_b=M_GA(index_beam); % Number of carriers in the selected beam
                %%%% Opt problem with explicit carriers
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,index_beam,users_index,Req_user_aux, P_GA(index_beam), M_b*obj.Delta_W);
                    R_off(users_index)=R_off_users.';
                end
            end
            
            
        end
        
        
        
    end
    
    
    %---------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------------------------------------------
    
    
    methods (Access=protected )
        
        
        function effective_snr=ComputEffectiveSNR (obj)
            % Compute effective snr per beam
            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       effective_snr: Vector with effective snr per beam
            
            effective_snr=zeros(1,obj.K);
            for index=1:obj.K
                users_index=obj.users_b_index{index};
                if ~isempty(users_index)
                    aux=obj.Gamma_wo_P( users_index,index)*obj.Max_P_sat/0.5  ;
                    effective_snr(index)=  geo_mean(aux) ;
                end
            end
            
        end
        
        
        
        function M_beam= BandwidthDisc(obj,W_cont,mode)
            % Obtain discrete carrier allocation for each beam from the obtained results( W_cont).
            % If two beams reuse the same carrier after the discretization, the conflict is resolve by allocating the carrier to the beam with higher traffic demand.
            % INPUT:
            %       obj: The class object that contains the scenario information
            %       W_cont: Continous bandwidth allocation
            %       mode: Value 1 when discretizing the bandwidh from the solution with flexible bandwidth allocation. Value 2 when discretizing the bandwidh from
            %             the solution with flexible bandwidth and beam-user mapping.
            % OUTPUT:
            %       M_beam: Discrete bandwidth. Vector with the number of carriers per beam
            
            
            % Process the input
            switch mode
                case 1 % Flexible Bandwidth
                    M_beam_aux= W_cont*obj.M/0.5;
                    M_beam_aux( obj.Nuser_beam==0)=0;
                case 2 % Flexible Bandwidth. Flexibe Mapping
                    M_beam_aux=W_cont;
            end
            
            
            % First step of discretization
            M_beam=floor(M_beam_aux); % Perform a floor operation as first step for the discretization.
            Number_carrier_conflict=  obj.M*obj.K-sum(M_beam); % Check if more carriers can be allocated or the first step violates the allowed numbers of carriers.
            [~,indS]=sort(M_beam_aux-M_beam,'descend'); % Order the beams regarding the unsastiefied bandwidth
            
            
            
            % Second step of the discretization. Carriers area allocated following the previous sorting order. The carrier allocation takes into account the bandwidth
            % of the adjacent beam.
            aux_cont=0; % Auxiliar variable.
            for index=1:obj.K
                if obj.Nuser_beam(indS(index))>0
                    if aux_cont==Number_carrier_conflict
                        break
                    end
                    switch  indS(index)
                        case 1
                            if  (M_beam(indS(index)+1)+ M_beam(indS(index))+1)<=2*obj.M
                                M_beam(indS(index))= M_beam(indS(index))+1;
                                aux_cont=aux_cont+1;
                            end
                        case obj.K
                            if  (M_beam(indS(index)-1)+ M_beam(indS(index))+1)<=2*obj.M
                                M_beam(indS(index))= M_beam(indS(index))+1;
                                aux_cont=aux_cont+1;
                            end
                        otherwise
                            if  (M_beam(indS(index)-1)+ M_beam(indS(index))+1 )<=2*obj.M && (M_beam(indS(index)+1)+ M_beam(indS(index))+1)<=2*obj.M
                                M_beam(indS(index))= M_beam(indS(index))+1;
                                aux_cont=aux_cont+1;
                            end
                    end
                end
                
            end
            
            
            
        end
        
        
        function  R_off_users=SecondStep(obj,N,M,ind_beam,users_index,Req_user_aux,P,B)
            %  Second Step of the two-step optimization process( User-Carrier assignment within the beams), as it
            % is described in "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021
            % INPUT:
            %       obj: The class object that contains the scenario information
            %       N: Number of users served by the selected beam
            %       M: Number of allocated carriers to the selected beam
            %       users_index: Auxiliar variable. It associates the users to the beam
            %       Req_user_aux: Vector with the user demanded rates
            %       P: Vector with the power per beam
            %       B: Vector with the bandwidth per beam ( Normalized to the total avaliable bandwidth)
            
            % OUTPUT:
            %       R_off_users: Offered user rate to the user withint the selected beam
            
            
            % Mixed Binary Quadratic Program (MBQP) problem, solved with CVX and MOSEK
            cvx_solver Mosek
            cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', obj.Max_time_opt_second,'MSK_DPAR_LOWER_OBJ_CUT',-1.0e30);
            cvx_begin quiet
            variable w(M,N)
            variable u(M,N) binary
            expression Roff_cvx(M,N)
            expression obj_fun
            for index_car=1:M
                for index_user=1:N
                    Roff_cvx(index_car,index_user)=w(index_car,index_user)*obj.Delta_W* log(1 +  obj.Gamma_wo_P(users_index(index_user),ind_beam)*P/B    )/log(2);
                end
            end
            obj_fun= ( Req_user_aux  - sum(Roff_cvx,1)).^2; % QU
            minimize( sum(obj_fun) );
            subject to
            u>=w
            sum(u,1)<= 1 % Constraint of one carrier per user
            sum(w,2)<=1 % Constraint of overall carrier time allocation
            0<=w<=1
            cvx_end
            % Store the obtained values
            R_off_users=sum(Roff_cvx,1);
            
        end
        
    end
    
end

%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------
% Auxiliar functions
%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------------


function [fval] = Flex_Power_Beam(x,Req_b,K,effective_snr_beam)
    % Objective function for the flexible power allocation, for the optimization process with fmincon.
    % INPUT:
    %       x: Solution input ( Power allocation in this case)
    %       Req_b: Requested traffic per beam
    %       K: Number of beams
    %       effective_snr_beam: Vector with the effective snr per beam
    % OUTPUT:
    %       fval: Evaluation of the objective function ( Quadratic Unment demand)


    % Process the input
    x=repelem(x,2)/2;
    h=zeros(1,K);

    % Measure the quadratic unmnet per beam
    Roff=zeros(1,K);
    for index=1:K
        if effective_snr_beam(index)~=0
            Roff(index)= 0.5*log2(  1+ effective_snr_beam(index)*x(index)   );
        else
            Roff(index)=0;
        end

        h(index)= ( Req_b(index)-Roff(index))^2;
    end

    % Aggregate the quadratic unmnet per beam
    fval= sum(h);



end

 


function [f] = ga_function(x,obj )
            % Objective function for genetic algorithm.
            % INPUT: 
            %       x: Solution input ( Power allocation in this case)
            %       obj: Struct with the simulated scenario information
            % OUTPUT: 
            %       f: Evaluation of the objective function ( Quadratic Unment demand)
            


    tol=1e-8;
    nvars=length(x);
    % Process the input
    aux_P=1:obj.K/2;
    aux_B=obj.K/2+1:nvars;

    % Cheking of the values is needed due to the GA implementation
    % Verify Bandwidth
    [B] = Check_B(x(aux_B),obj);
    % Carrier power to beam power
    P=repelem(x(aux_P),2).*B*2*obj.M;
    P_amp=zeros(1,obj.K/2);
    for index=1:obj.K/2
        P_amp(index)=sum(P( 2*(index-1)+1:2*index));
        if    P_amp(index)> obj.Max_Pamp
            k=  obj.Max_Pamp/P_amp(index);
            P( 2*(index-1)+1:2*index)= k*P( 2*(index-1)+1:2*index);
        end
    end
    if sum(P_amp)>obj.Max_P_sat
        k=  obj.Max_P_sat/sum(P_amp);
        P= k*P;
    end

 % Obtain offered traffic per beam
    R_off_beam=zeros(1,obj.K);
    for index=1:obj.K

        users_index=find(  obj.user_beams==index);
        if ~isempty(users_index)
            if B(index)~=0
                aux=1+  obj.aux_Cn( users_index,index)*P(index)/(B(index));
            else
                aux=1;
            end
            snr_eff2_b=  geo_mean(aux) ;
        else
            snr_eff2_b=1;
        end
        if B(index)~=0
            R_off_beam(index)=B(index)* log(  snr_eff2_b   )/log(2);
        else
            R_off_beam(index)=0;
        end

    end

 % Measure the Quadratic Unmet demand
    f= sum ( (R_off_beam-obj.Req_beam).^2);


    if f<tol
        f=0;
    end

 

end



function [B] = Check_B(B_in,obj)
    % Verification and correction of the allocated bandwidth. This is described int Appendix B from as "Flexible User Mapping 
    % for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021. Is based on the process described in " A Genetic
    % Algorithm for Joint Power and Bandwidth Allocation in Multibeam Satellite Systems", IEE Aeroespace Conference, 2019
    % INPUT:
    %       B_in: Allocated bandwidth
    %       obj: Class Object
    % OUTPUT:
    %       B:  Allocated bandwidth after verification


    % Process the input
    B=B_in/(2*obj.M);

    % Selects randomly a processing order
    coin_flip=rand(1);

    if coin_flip>0.5
        order_lim=1:length(B)-1;
        step=1;
    else
        order_lim=length(B):-1:2;
        step=-1;
    end


    % Ensure no bandwidth overlaping in neighbour beams
    for ind=order_lim

        if (B(ind)+B(ind+step))>1

            B(ind)=1-B(ind+step);
        end

    end


    % Ensure that the all avaliable bandwidth is employed
    [~,order_req_Beam]=sort(obj.Req_beam,'descend');

    for ind_center_b=order_req_Beam

        ind_left=ind_center_b-1;
        ind_right=ind_center_b+1;

        B_center=B(ind_center_b);

        if ind_left>0
            B_left=B(ind_left);
        else
            B_left=0;
        end

        if ind_right<=length(B)
            B_right=B(ind_right);
        else
            B_right=0;
        end


        B_un= 1- B_center -max(B_right,B_left);

        if  (B_center+B_un)> obj.ub_B(ind_center_b)/(2*obj.M)

            B(ind_center_b)= obj.ub_B(ind_center_b)/(2*obj.M);

        else
            B(ind_center_b)= B_center+B_un;
        end



    end




end



