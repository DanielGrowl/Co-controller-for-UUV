function [Error_state_RMSE] = RMSE(Error_state)
    Error_state_RMSE =[];
    Error_total = 0;
    length = size(Error_state,2);
    start = 800;
    for i=start:length
        Error_total = Error_total + Error_state(:,i).^2;
        Error_RMSE = Error_total./i;
        Error_state_RMSE = [Error_state_RMSE, Error_RMSE];
    end
end