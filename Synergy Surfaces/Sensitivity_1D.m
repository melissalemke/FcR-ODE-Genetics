% Arnold Lab, University of Michigan
% Robby Theisen, PhD Student
% Last edit: August 10, 2021
function [yout] = Sensitivity_1D(params, paramnames, dparams, p1)
%% Iterate
yout = zeros(1, length(dparams));
for i = 1:length(dparams)
    newP = params;
    newP(p1) = dparams(i) * newP(p1);
    [yend, ~] = Simulate(newP, paramnames, [], []);
    yout(i) = yend(33);
end


end