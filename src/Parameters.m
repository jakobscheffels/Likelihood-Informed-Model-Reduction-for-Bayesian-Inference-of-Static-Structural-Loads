function Parameters(inputName, inputVal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Input: inputName is cell containing variable names to be modified
    %%%            inputVal is list containing values of variables
    %%%
    %%%     Output: Stores modified parameters in Parameter.mat to be loaded
    %%%
    %%%     Valid variable names (default value):
    %%%       mu_q: mean of distributed load (1)
    %%%       sigma_q: standard deviation of distributed load (0.2)
    %%%       L: length of bar element (2)
    %%%       theta: correlation length (2)
    %%%       D: rigidity of bar (4e8)
    %%%       E: Youngs modulus of beam (10000)
    %%%       nele: number of finite elements used (100)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scale of fluctuation of q (=2*correlation length)
    theta = 2;
    % number of elements
    nele = 100;
    nnode = nele+1;
    % length of domain
    L = 2;

    D = 4e8;
    E = 10000;
    BC_dofs = 1;

    beam_bool = false;
    fixed = {'true'};
    % distributed load
    mu_q = 4e6;
    sigma_q = 0.3*mu_q;
    

    moments = false;
    tunnel = false;
    k_vector = zeros(nnode,1);
    zeta = 1/7;
    t = 0.35;
    % set Parameters
    if nargin > 0
        for i=1:size(inputName,2)
            input = inputVal(i);
            if strcmp(inputName(i),'theta')
                theta=input{1};      
            elseif strcmp(inputName(i),'nele')
                nele=input{1};
                nnode=nele+1;
            elseif strcmp(inputName(i),'L')
                L=input{1};
            elseif strcmp(inputName(i),'D')
                D=input{1};
            elseif strcmp(inputName(i),'E')
                E=input{1};
            elseif strcmp(inputName(i),'I')
                I=input{1};
            elseif strcmp(inputName(i),'mu_q')
                mu_q=input{1};
            elseif strcmp(inputName(i),'sigma_q')
                sigma_q=input{1};
            elseif strcmp(inputName(i),'beam_bool')
                beam_bool=input{1};
            elseif strcmp(inputName(i),'fixed')
                fixed=input{1};
            elseif strcmp(inputName(i),'BC_dofs')
                BC_dofs=sort(input{1});
            elseif strcmp(inputName(i),'moments')
                moments=input{1};
            elseif strcmp(inputName(i),'tunnel')
                tunnel=input{1};
            elseif strcmp(inputName(i),'k_vector')
                k_vector=input{1};
            elseif strcmp(inputName(i),'zeta')
                zeta=input{1};
            elseif strcmp(inputName(i),'t')
                t=input{1};
            else
                inputStr = inputName(i);
                inputStr = inputStr{:};
                text = append('Invalid variable name: ',inputStr,'! Using default parameter instead!');
                warning(text)
            end
        end
    end
    
    if ~beam_bool
        moments=false;
    end

    l=L/nele;
    x_dofs=0:l:L;
    
    if beam_bool
        index_disp = 1:2*nnode;
    else
        index_disp = 1:nnode;
    end

    L_mat = loadMapper(beam_bool,nnode,nele,l);

    % Store parameters in 'Parameters.mat'
    save('Parameters','BC_dofs','D', ...
        'E','l','L','L_mat','mu_q','nele','nnode','sigma_q','theta','tunnel',...
        'k_vector','x_dofs','beam_bool','fixed','index_disp','moments',...
        'zeta','t')
end