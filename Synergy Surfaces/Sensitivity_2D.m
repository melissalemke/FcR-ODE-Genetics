% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Last edit: March 11th, 2020

function [YOutput, filettl] = Sensitivity_2D_1v3(params, dParams_p1, dParams_p2,...
    paramnames,ycol,p1,p2)
% If not inside a PBS job, use 4 processors
if isempty(getenv('PBS_NP'))
    NP = 4; % May want to change this back if it takes too long - talk to Melissa
else
    NP = str2double(getenv('PBS_NP'));
end

myPool = parpool('local',NP);

% Create Output Matrices
    [YOutput] = iterate([1 2],dParams_p1,dParams_p2,params,paramnames,ycol,p1,p2);
    clear Outputs;
    filettl = strrep(strcat('8_9_2021','_','2D_',paramnames(p1),'_vs_',paramnames(p2),'_Surface.mat'),' ','_');

delete(myPool)
%% Iterate Through Parameter Changes
      function [YOutput] = iterate(jloop,dParams_p1,dParams_p2,params,paramnames,ycol,p1,p2)
        YOutput = zeros(length(dParams_p1)*length(dParams_p2),ycol);
        lendParams = length(dParams_p1);
        lenjloop = length(jloop);
       steadystate_total = 0;
        parfor r = 1:lendParams*lendParams
            i = ceil(r/lendParams);
            if rem(r,lendParams)
                j=rem(r,lendParams);
            else
                j=lendParams;
            end
            newP = params;

             newP(p1) = dParams_p1(i)*params(p1);
            newP(p2) = dParams_p2(j)*params(p2); %#ok<PFBNS>
            [yout steadystate] = Simulate(newP,paramnames, [],[]);
            YOutput(r,:) = yout;
            steadystate_total = steadystate_total+steadystate;
            if steadystate>0
              disp([paramnames(p1)+" at "+num2str(dParams_p1(i))+"X "+paramnames(p2)+" at "+num2str(dParams_p2(j))+"X did not run to steady state"])
else
disp([paramnames(p1)+" at "+num2str(dParams_p1(i))+"X "+paramnames(p2)+" at "+num2str(dParams_p2(j))+"X ran to steady state"])
                   end  

        end
        if steadystate_total >0
    disp(["Some simulations did not reach steady state"])
else
    disp(["All simulations at steady state"])
end
    end
end            
