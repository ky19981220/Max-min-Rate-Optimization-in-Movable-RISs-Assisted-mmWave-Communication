function [stopnow, reason] = mystopfun(problem, x, info, last)
    stopnow = (last >= 5 && info(last-4).cost - info(last).cost < 1e-10);
    reason = 'two successive iterations combined have decreased the cost by less than 1e-10';
    % if xyz % decide if solver should stop based on inputs
    %     stopnow = true;
    %     reason = 'This optional message explains why we stopped.';
    % end
end