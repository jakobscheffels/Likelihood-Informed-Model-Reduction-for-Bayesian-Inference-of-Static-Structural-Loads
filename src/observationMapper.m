function B_obs = observationMapper(nm_pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Function to assemble mapping matrices
%%%
%%%     Input:  nm_pos      Position indices of observable nodes
%%%
%%%     Output: B_obs       Operator to extract observable nodes from full
%%%                         state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('Parameters.mat','beam_bool','nnode')
    
    nm = length(nm_pos);

    if beam_bool
        B_obs = zeros(nm,2*nnode);
        j=1;
        for i=1:nm
            B_obs(j,2*(nm_pos(i)-1)+1)=1;
            j=j+1;
        end
    else
        B_obs = zeros(nm,nnode);
        j=1;
        for i=1:nm
            B_obs(j,nm_pos(i))=1;
            j=j+1;
        end        
    end

end