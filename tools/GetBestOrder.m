function [ odorOrder, ORNOrder ] = GetBestOrder( data3D )
% use simulated annealing algorithm to find the order of ORN and odor to
% gather the large responses close to the diagonal. 

%% design a punishment matrix
[odorNum, ornNum, zDim] = size(data3D);

pM = zeros(odorNum, ornNum);
for i = 1: odorNum
    for j = 1: ornNum
        pM(i, j) = abs(ornNum*i - odorNum*j)/sqrt(odorNum^2 + ornNum^2); % defined using the distance to the diagonal
    end
end
pM = repmat(pM, [1 1 zDim]);

%% define the initial order
odorOrderInit = 1:odorNum;
ORNOrderInit = 1:ornNum;

%% simulated annealing
parent = cell(1,2);     %initialize the permutation
parent{1} = odorOrderInit;  parent{2} = ORNOrderInit;

% main settings
Tinit = 1;      % initial temp
minT = 1e-10;	% stopping temp
cool = @(T) (.975*T);	% annealing schedule
minF = 0;
max_consec_rejections = 50000;
max_try = 50000;
max_success = 5000;
k = 1;	% boltzmann constant

% counters etc
itry = 0;
success = 0;
finished = 0;
consec = 0;
T = Tinit;
initenergy = CostFun(parent, pM, data3D, zDim);
oldenergy = initenergy;
total = 0;
% disp('----------Simulated Annealing Process----------');
% fprintf(1,'\n  T = %7.5f, loss = %10.5f\n',T,oldenergy);

while ~finished
    itry = itry+1; % just an iteration counter
    current = parent; 
    
    % Stop / decrement T criteria
    if itry >= max_try || success >= max_success
        if T < minT || consec >= max_consec_rejections
            finished = 1;
            total = total + itry;
            break;
        else
            T = cool(T);  % decrease T according to cooling schedule
%             fprintf(1,'  T = %7.5f, loss = %10.5f\n',T,oldenergy);
            total = total + itry;
            itry = 1;
            success = 1;
        end
    end
    
    newparam = getNewSeq( current );
    newenergy = CostFun(newparam, pM, data3D,zDim);
    
    if (newenergy < minF)
        parent = newparam; 
        oldenergy = newenergy;
        break
    end
    
    if (oldenergy-newenergy > 1e-6)
        parent = newparam;
        oldenergy = newenergy;
        success = success+1;
        consec = 0;
    else
        if (rand < exp( (oldenergy-newenergy)/(k*T) ))
            parent = newparam;
            oldenergy = newenergy;
            success = success+1;
        else
            consec = consec+1;
        end
    end
end

minimum = parent;
odorOrder = minimum{1};  
ORNOrder = minimum{2};
fval = oldenergy;

disp('----------Search Odor and ORN Orders----------');
disp('Simulated Annealing Results:');
fprintf('%30s %.0f \n', 'Initial temperature:',Tinit);
fprintf('%30s %.6f \n', 'Final temperature:',T);
fprintf('%30s %.0f \n', 'Consecutive rejections:',consec);
fprintf('%30s %.0f \n', 'Number of function calls:',total);
fprintf('%30s %.2f \n', 'Total final loss:',fval);
fprintf('%30s ','ORN order: '); disp( num2str(ORNOrder));
fprintf('%30s ','Odor order: '); disp( num2str(odorOrder));
end

%define new sequence generator
function newSeq = getNewSeq(oldSeq )
    %first randomly decide change row or col
    index = round(rand) + 1;

    len = length(oldSeq{index});
    %get two random ids to swap
    swap_id =  unidrnd(len,2,1); %generate two random ids 
    while swap_id(1) == swap_id(2) %make sure the two ids are unequal
        swap_id =  unidrnd(len,2,1);
    end
    
    %swap the two elements in the sequence
    newSeq = oldSeq;
    e1 = oldSeq{index}(swap_id(1));
    e2 = oldSeq{index}(swap_id(2));
    newSeq{index}(swap_id(1)) = e2;
    newSeq{index}(swap_id(2)) = e1;
end

%define cost function
function cost = CostFun(seq, pM, dataM, zDim)
    if zDim == 1
        dataMNew = dataM(seq{1}, seq{2} );
        cost = sum(sum(dataMNew(:,:).*pM(:,:)))/sum(dataMNew(:));
    else
        dataMNew = dataM(seq{1}, seq{2}, :);
        cost = 0;
        for i = 1:zDim
            cost = cost + sum(sum(dataMNew(:,:,i).*pM(:,:,i)))/sum(sum(dataMNew(:,:,i)));
        end
    end
end