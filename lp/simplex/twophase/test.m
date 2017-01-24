warning('off', 'all')

%% Test correctness
disp('Test correctness')
nTest = 1000;

diff = inf*ones(nTest,1);

%opt = optimoptions('linprog', 'display', 'off');

opt = optimset('Display', 'off');

dwork      = zeros((10+10)*(10+4)+10,1);
iwork      = zeros((10+10)*3+10,1);

i = 1;
while i <= nTest
    size = 6;
    m    = 2;
    c    = randn(size,1);
    Aeq  = randn(m,size);
    beq  = randn(m,1);
    lb   = randn(size,1)-2;
    ub   = randn(size,1)+2;

    if any(lb>ub)
        i = i-1;
        continue
    end

    [x2, fval2, inform2] = linprog(c,[],[],Aeq,beq,lb,ub, opt);
    
    if inform2 ~= 1
        i = i-1;
        continue
    end
    
    [x1,fval1,inform1] = mxbdlpext(c,Aeq,beq,lb,ub,dwork,iwork);

    diff(i) = x1'*c-x2'*c;

    if diff(i) > 0.00001
        save 'data.mat' c Aeq beq lb ub
        error('diff(fval) is too large')
    end
    
    i = i+1;
end
% save 'data.mat' c Aeq beq lb ub
disp('correct!')

%% Test speed
nTest = 1000;
size  = 10;
m     =  5;
c_list     = cell(nTest,1);
Aeq_list   = cell(nTest,1);
beq_list   = cell(nTest,1);
lb_list    = cell(nTest,1);
ub_list    = cell(nTest,1);

% generate data
disp('generate data for speed testing')
i = 1;
while i <= nTest
    c    = randn(size,1);
    Aeq  = randn(m,size);
    beq  = randn(m,1);
    lb   = randn(size,1)-2;
    ub   = randn(size,1)+2;

    if any(lb>ub)
        i = i-1;
        continue
    end

    [x2, fval2, inform2] = linprog(c,[],[],Aeq,beq,lb,ub, opt);
    
    if inform2 ~= 1
        i = i-1;
        continue
    end
    
    c_list{i}   =   c;
    Aeq_list{i} = Aeq;
    beq_list{i} = beq;
    lb_list{i}  =  lb;
    ub_list{i}  =  ub;
    
    i = i+1;
end

disp('Time used by Linprog:')
tic
for i = 1:nTest
    [x2, fval2, inform2] = linprog(c_list{i},[],[],Aeq_list{i},beq_list{i},lb_list{i},ub_list{i}, opt);
end
toc

disp('Time used by mxbdlp:')
tic
for i = 1:nTest
    %c = c_list{i};
    %Aeq = Aeq_list{i};
    %beq = beq_list{i};
    %lb  = lb_list{i};
    %ub  = ub_list{i};
    %save 'tmp.mat' i c Aeq beq lb ub
    %pause(0.00001)
    %i
    %save 'data.mat' c Aeq beq lb ub
    x = mxbdlpext(c_list{i},Aeq_list{i},beq_list{i},lb_list{i},ub_list{i}, ...
               dwork, iwork);
    %[x1,fval1,inform1] = mxbdlp(c,Aeq,beq,lb,ub);
end
toc

%save 'speeddata.mat' c_list Aeq_list beq_list lb_list ub_list