function sol = ChimeraModel(params,domain,ins)
    domain = domain(2);         % Simulation domain
    nmodels = ins.nmodels;
    M = zeros(23*nmodels, domain + 1);  % Simulation matrix
    npars = length(params)/nmodels;
    j = 0:nmodels - 1;
    
    %Configurable initial conditions
    M(4+j*23,1) = params(2 + j*npars);       % E_3
    M(14+j*23,1) = params(4 + j*npars);      % J_1
    M(16+j*23,1) = params(5 + j*npars);      % J_L
    M(17+j*23,1) = params(6 + j*npars);      % Dv
    M(20+j*23,1) = params(7 + j*npars);      % R_j
    M(21+j*23,1) = params(8 + j*npars);      % RJ_L
    
    %Parameters
    beta_L = params(9 + j*npars);       % beta_L
    beta_T = params(10 + j*npars);      % beta_T 
    beta_P = params(11 + j*npars);      % beta_P
    phi_EP = params(12 + j*npars);      % phi_{EP}
    lambda_fq = params(13 + j*npars);   % lambda_fq_{fq}
    vartheta_E = params(14 + j*npars);  % varvartheta_E_E
    gamma_L = params(15 + j*npars);     % gamma_L_L
    k_L  = params(16 + j*npars);        % k_L
    k_P  = params(17 + j*npars);        % k_P
    phi_T   = params(18 + j*npars);     % phi_T
    lambda_qf = params(19 + j*npars);   % lambda_{qf}
    psi_e = params(20 + j*npars);       % psi_e
    phi_PH = params(21 + j*npars);      % phi_{PH}
    delta = params(22 + j*npars);       % delta
    
    %Inputs
    m = params(23 + j*npars);           %Number of inmigrants
    
    %Modifying parameters
    eta_L = params(24 + j*npars);       % eta_L
    vartheta_P = params(25 + j*npars);  % vartheta_P
    eta = params(26 + j*npars);         % gamma_H
    z = params(27 + j*npars);           % z
    phi_PL = params(28 + j*npars);      % phi_{PL}
    eta_vartheta = params(29 + j*npars);% eta_vartheta
    nons = params(30 + j*npars);        % initial proportion of inmune population
    
    %Difussion system parameters
    a_L = params(31 + j*npars);         % a_L
    b_L = params(32 + j*npars);         % b_L
    a_mu = params(33 + j*npars);        % a_mu
    b_mu = params(34 + j*npars);        % b_mu
    mu = params(35 + j*npars);          % mu
    a_H = params(36 + j*npars);         % a_H
    b_H = params(37 + j*npars);         % b_H
    nu = params(38 + j*npars);          % nu
    lon = round(params(39 + j*npars));  % time
    
    lambda_H = lambda_fq.^eta_L;
    alpha_H = lambda_qf.^(1./eta_L);
    vartheta_E_L = vartheta_P;
    vartheta_E_H = vartheta_E_L.^eta_vartheta;

    %Parameters for building the difussion matrices
    ML = cell(1,nmodels);
    ML_j = cell(1,nmodels);
    ML_gamma_L = cell(1,nmodels);
    MJL = cell(1,nmodels);
    MJL_gamma_L = cell(1,nmodels);
    MJH = cell(1,nmodels);
    MJH_mu = cell(1,nmodels);
    MJH_gamma_L = cell(1,nmodels);
    
    for ii = 1:nmodels
        dgamma_L = distri_est(a_L(ii), b_L(ii), lon(ii));
        dmu = distri_mu(a_mu(ii), b_mu(ii), lon(ii));
        dgamma_L2 = distri_est(a_H(ii), b_H(ii), lon(ii));
        dgamma_L = dgamma_L * gamma_L(ii);
        dgamma_L2 = dgamma_L2 * eta(ii);
        dmu = dmu * mu(ii);
        
        len = length(dgamma_L);
        
        [ML{ii}, ML_j{ii}, ML_gamma_L{ii}] = L_DS(dgamma_L, vartheta_E_H(ii), lambda_H(ii), alpha_H(ii));
        [MJL{ii}, MJL_gamma_L{ii}] = JL_DS(dgamma_L);
        [MJH{ii}, MJH_mu{ii}, MJH_gamma_L{ii}] = JH_DS(dmu, dgamma_L2);
    end
    
    L = zeros(nmodels, len * 2);
    JH = zeros(nmodels, len);
    JL = JH;
    L_gamma_L=zeros(nmodels,1);
    L_j = JH;
    JI_gamma_L = L_gamma_L;
    JI_mu = L_gamma_L;
    JA_gamma_L = L_gamma_L;
    

    JH(:,1) = params(4 + j*npars);    
    JL(:,1) = params(5 + j*npars);
    
    M(1+j*23, 1) = params(1 + j*npars) .* (1 - nons);      % N_0
    M(9+j*23, 1) = params(3 + j*npars) .* delta;           % A1
    M(11+j*23, 1) = params(3 + j*npars) .* (1 - delta);    % A1
    M(18+j*23, 1) = params(1 + j*npars) .* nons + params(7 + j*npars) + params(8 + j*npars); % Recovered or inmune initial population
    
    k_L = k_P .* k_L;
    M(19+j*23, 1) = dot((M(9+j*23, 1)+M(11+j*23, 1)),k_P) + dot(M(7+j*23, 1),k_L);       % T

    for i=1:domain
       
        % Population definition    
        E_f = sum(M([3+j*23,4+j*23], i));                     % Free exposed
        N_f = sum(M([1+j*23,18+j*23], i)) + E_f;     % Total population in free circulation
        N_t = N_f;
        tota = zeros(nmodels,1);
        H_f = tota;
        L_f = tota;
        C = 0;
        
        for k = 0:nmodels-1
            N_q = sum(M([2,5:6,8,10,12,13:16,22:23]+k*23, i));   % Total population in quarantine
            H_f(k+1) = sum(M(7+k*23, i));                                  % Free infected
            L_f(k+1) = sum(M([9+k*23,11+k*23], i));                  % Free asymptomatic
            N_t = N_t+N_q; % Total population 
            tota(k+1) = sum(M(7+k*23:16+k*23, i)); % Total infectious
            C = C + sum(M([1:18,22:23]+k*23,1))*k_L(k+1);          % Max virus concentration
        end       
        
        % Contact probabilities definition
        N_t = N_t+sum(H_f+L_f); % Total population 
        Ttota = sum(tota);
        
        aleph = 1 + nu*Ttota/(sum(M(1+j*23, i)) + Ttota);
        probAux = M(1+j*23, i) > 1;
        
        
        red = (M(1+j*23, i)/N_t).^aleph;
        PA = dot(z,L_f)*red;
        PI = dot(H_f,z)*red;

        % Probability of infection with pre-symptomatic individuals
        probL = real(1 - ((M(1+j*23, i) - 1)./M(1+j*23, i)).^PA).*probAux;     % 
        % Probability of infection with symptomatic individuals
        probH = real(1 - ((M(1+j*23, i) - 1)./M(1+j*23, i)).^PI).*probAux;

        
        % Probability of infection with the environmental repository
        Phi_T = min(((M(19+j*23, i)/C).^aleph), 1);
      
        % Susceptible population in free circulation
        M(1+j*23, i + 1) = (1 - psi_e).* m + lambda_qf.* M(2+j*23, i) + ... 
                   (1 - beta_P.* probL).*...
                   (1 - beta_L.* probH).* (1 - Phi_T.* beta_T).*...
                   (1 - lambda_fq).* M(1+j*23, i);                               % S_f
               
        % Susceptible population in quarantine
        M(2+j*23, i + 1) = (1 - (beta_P).* probL).* (1 - beta_L.* probH).* (1 - Phi_T.* beta_T).*...
                    (lambda_fq).* M(1+j*23, i) + (1 - lambda_qf).* M(2+j*23, i);      % S_q            
               
        % Exposed population in free circulation
        M(3+j*23, i + 1) = (beta_P.* probL + (1 - (beta_P.* probL)).* (beta_L.* probH)...
                    + (1 - (beta_P.* probL)).* (1 - (beta_L.* probH)).*...
                    Phi_T.* beta_T).* (1 - lambda_fq).* M(1+j*23, i);            % E_f^1

        M(4+j*23, i + 1) =  psi_e.* m + lambda_qf.* M(5+j*23, i) + (1 - lambda_fq).* M(3+j*23, i)...
                     + exp( - phi_EP).* ((lambda_qf).* M(6+j*23, i) + ...
                     (1 - lambda_fq).* M(4+j*23,i));                             % E_f^2
        
        %Exposed in quarantine
        M(5+j*23, i + 1) = (beta_P.* probL + (1 - (beta_P.* probL)).* (beta_L.* probH)...
                    + (1 - (beta_P.* probL)).* (1 - (beta_L.* probH)).* ...
                    Phi_T.* beta_T).* (lambda_fq).* M(1+j*23, i);                %E_q^1
                        
        M(6+j*23, i + 1) = (1 - lambda_qf).* M(5+j*23, i) + (lambda_fq).* M(3+j*23, i) +...
                  exp(-phi_EP).* ((1 - lambda_qf).* M(6+j*23, i) + ...
                  (lambda_fq).* M(4+j*23, i));                                   %E_q^2
        
        % Low symptomatic population (Presymptomatic + Symptomatic)
        M(11+j*23, i + 1) = (1 - delta).* (1 - vartheta_E) ...
                    .* (1 - exp( - phi_EP)).* ((1 - lambda_fq).* M(4+j*23, i)...
                    + (lambda_qf).* M(6+j*23, i)) + (1-vartheta_E_L).* (exp( - phi_PL))...
                    .* ((1 - lambda_fq).* M(11+j*23, i) + lambda_qf.* M(12+j*23, i)); %A1
        
        M(12+j*23, i + 1) = (1 - delta).* (1 - vartheta_E) ...
            .* (1 - exp( - phi_EP)).* ((lambda_fq).* M(4+j*23, i)...
            + (1 - lambda_qf).* M(6+j*23, i)) + (1 - vartheta_E_L)...
            .* (exp( - phi_PL)).* (lambda_fq.* M(11+j*23, i)...
            + (1 - lambda_qf).* M(12+j*23, i));%QA_1       
        
        % Identified Low symptomatic population (as difussion system)
        for k = 1:nmodels
            L_gamma_L(k) = sum(L(k,:) * ML_gamma_L{k}) + sum(L(k,end-1:end));
            aux = L(k,:) * ML_j{k};
            L_j(k,:) = aux(1:2:end) + aux(2:2:end);
            L(k,:) = L(k,:) * ML{k};
            JI_gamma_L(k) = sum(JH(k,:) * MJH_gamma_L{k});
            JI_mu(k) = sum(JH(k,:) * MJH_mu{k}) + JH(k,end);
            JH(k,:) = JH(k,:) * MJH{k};
            JA_gamma_L(k) = sum(JL(k,:) * MJL_gamma_L{k}) + JL(k,end);
            JL(k,:) = JL(k,:) * MJL{k};
        end
        
        L(:,1:2) = [(1 - lambda_fq) .* M(11+j*23, i) + lambda_qf .* M(12+j*23, i),lambda_fq .* M(11+j*23, i)...
            + (1 - lambda_qf) .* M(12+j*23, i)] .* (1 - vartheta_E_L) .* (1 - exp( - phi_PL));
        
        M(7+j*23, i + 1) = sum(L(:,1:2:end),2);
        M(8+j*23, i + 1) = sum(L(:,2:2:end),2);
            
        % High symptomatic population (Presymptomatic + Symptomatic)
        M(9+j*23, i + 1) = delta.* (1 - vartheta_E).* (1 - exp( - phi_EP))...
            .* ((1 - lambda_fq).* M(4+j*23, i)...
            + (lambda_qf).* M(6+j*23, i)) + (1 - vartheta_E_L).* (exp( - phi_PH))...
            .* ((1 - lambda_fq).* M(9+j*23, i)...
            + lambda_qf.* M(10+j*23, i)); %I1
        
        M(10+j*23, i + 1) = delta.* (1 - vartheta_E).* (1 - exp( - phi_EP))...
            .* ((lambda_fq).* M(4+j*23, i) + (1 - lambda_qf)...
            .* M(6+j*23,i)) + (1 - vartheta_E_L).* (exp( - phi_PH))...
            .* (lambda_fq.* M(9+j*23, i) + (1 - lambda_qf).* M(10+j*23, i));%QI1

        %Detecciones JH  
        M(13+j*23, i + 1) = vartheta_E_L .* (exp( - phi_PH)) .* (M(9+j*23, i)...
                           + M(10+j*23, i)) + M(13+j*23, i) .* exp( - phi_PH);%JI_1
        
        M(22+j*23, i + 1) = (M(6+j*23, i) + M(4+j*23, i)) .* vartheta_E...
                           .* (1 - exp( - phi_EP)) .* delta + M(22+j*23, i) .* exp( - phi_PH);
        
        JH(:,1) = (M(13+j*23, i) + M(22+j*23, i)) .* (1 - exp( - phi_PH))...
              + (1 - exp( - phi_PH)) .* (M(9+j*23, i) + M(10+j*23, i));%J2
        
        M(14+j*23, i + 1) = sum(JH,2);%JI_2 en adelante
        
        %Detecciones de JL
        M(15+j*23, i + 1) = vartheta_E_L .* (exp( - phi_PL)) .* (M(11+j*23, i)...
                           + M(12+j*23, i)) + M(15+j*23, i) .* exp( - phi_PL);%JA_1
        
        M(23+j*23, i + 1) = (M(6+j*23, i) + M(4+j*23, i)) .* vartheta_E...
                           .* (1 - exp( - phi_EP)) .* (1 - delta) + M(23+j*23, i) .* exp( - phi_PL);
              
        JL(:,1) = (M(15+j*23, i) + M(23+j*23, i)) .* (1 - exp( - phi_PL))...
              + vartheta_E_L .* (1 - exp( - phi_PL)) .* (M(11+j*23, i) + M(12+j*23, i));%J_A2
        JL = JL + L_j;
        M(16+j*23, i + 1) = sum(JL,2);%JA_2 en adelante
        
        
        %D, R y V
        M(17+j*23, i + 1) = M(17+j*23, i) + JI_mu;%D
        M(18+j*23, i + 1) = M(18+j*23, i) + JA_gamma_L + JI_gamma_L + L_gamma_L; %R
        
        M(19+j*23,i+1) = (1 - min(M(19+j*23, i)/C, 1))* (dot(k_L,H_f) + dot(k_P,L_f))...
                    + exp(-phi_T).* M(19+j*23, i); %T  
        M(20+j*23, i + 1) = M(20+j*23, i) + JI_gamma_L; %RI_j
        M(21+j*23, i + 1) = M(21+j*23, i) + JA_gamma_L; %RA_j
        
%        for k=0:nmodels-1
%             tester(k+1,i)=sum(M([1:18,22:23]+k*23,i));
%        end
    end
%     tester  
    
    outs = zeros(20*nmodels,domain+1);

    for i = 1:nmodels
        k = i-1;
        outs(20*(k)+1:20*i,:) = [M([1,2]+k*23,:); sum(M([3,4]+k*23,:)); sum(M([5,6]+k*23,:));...
            sum(M([9,11]+k*23,:)); sum(M([10,12]+k*23,:));M([7,8]+k*23,:); sum(M([13:16,22:23]+k*23,:));...
            M((17:21)+k*23,:);sum(M([20,21]+k*23,:)); sum(M([7:16,22:23]+k*23,:)); M(14+k*23,:);...
            sum(M([13,15:16,22:23]+k*23,:));sum(M((13:16)+k*23,:));sum(M((22:23)+k*23,:))];
        
    end
    
    sol = struct();
    sol.x = 0:domain;
    sol.y = outs;      
end
