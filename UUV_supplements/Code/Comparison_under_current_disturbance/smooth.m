function [K_co_state_smooth] = smooth(K_co_state)
    K_co_smooth =  [1;1;1;1;1;1;1;1;1;1;1;1];
    K_co_state_smooth =[K_co_smooth];
    length = size(K_co_state,2);
    for i=2:length
        K_co_smooth = K_co_smooth * 0.95 + K_co_state(:,i).*0.05;
        K_co_state_smooth = [K_co_state_smooth, K_co_smooth];
    end
end