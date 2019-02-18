function [TrainedSeqs] = seqwin_targetfile(SubCode,what,varargin)

% TrainedSeqs = seqwin_targetfile(SubCode , 'Behav&Memo','baseDir','...');

% Creates the tatgert files for the memory trainig blocks as well as the behavioral blocks.
% saves them in two separate folders, in the baseDir directory  
% Returns TrainedSeqs the set of 4 trained sequences per subject
% Neda Kordjazi - July 2018
%% variables and input arguments
% hand = 2; %right
TrainedSeqs = [];
MaxPress = 9;
baseDir = '/Users/nkordjazi/Documents/SeqWin/TargetFiles';
c = 1;

while(c<=length(varargin))
    switch(varargin{c})
        case {'hand'}
            %right (2 , default) or left
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'TrainedSeqs'}
            % [4x9]  the 4 trained sequences can be provided to the program
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'MaxPress'}
            % maximum number of presses in the sequence
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'baseDir'}
            % the folder where the target files are going to be saved.
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end


%%
cd(baseDir)
rng('shuffle');
switch what
    case 'Behavioral'
        %% check if the trained sequecnes are provided, else generate them
        reGen = 'n';
        if isempty(TrainedSeqs)
            while reGen == 'n'
                TrainedSeqs(1:4,1) = [1:4]'; % make sure each sequecne starts with a different number
                for seqNum = 1:4
                    seq = [TrainedSeqs(seqNum,1) , sample_wor([1:5],1,MaxPress-1)];
                    tempSeq = diff(seq);
                    while sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                        seq = [TrainedSeqs(seqNum,1) , sample_wor([1:5],1,MaxPress-1)];
                        tempSeq = diff(seq);
                    end
                    TrainedSeqs(seqNum , 1:MaxPress) = seq;
                end
                TrainedSeqs
                reGen = input('Happy with the sequences? (y/n).... (n will regenerate sequences)', 's');
            end
        end
        Fname  = [SubCode,'_', what, '_tgtFiles'];
        mkdir(Fname)
        OrderFields = {'seqNumb','press1','press2','press3','press4','press5','press6','press7','press8','press9','hand','cueP','iti','sounds' , 'Window'};
        %% creat random blocks
        clear RndmSeq
        Window = [1:8];
        Trials = 1:5; % Trials per window which is 40 altogether
        for BN = 1:20
            RndmSeq = [];
            for e = 1:length(Window)
                A.seqNumb(1:length(Trials),1) = 0;
                A.Window(1:length(Trials),:) = Window(e);
                for i = 1:length(Trials)
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    while sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                        seq = sample_wor([1:5],1,MaxPress);
                        tempSeq = diff(seq);
                    end
                    A.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
                    A.iti(1:i,:) = 1000;
                    A.hand(i,:) = 2;
                    A.sounds(i,:) = 1;
                    for press= 1:MaxPress
                        comnd  = [' A.press' , num2str(press) , '(i,1) = seq(press);'];
                        eval(comnd);
                    end
                end
                RndmSeq = addstruct(RndmSeq , A);
                clear A
            end
            idxAll = randperm(length(RndmSeq.Window));
            for kk = 1:length(OrderFields)
                statement = ['RndmSeq.' , OrderFields{kk},' = ','RndmSeq.' , OrderFields{kk},'(idxAll);'];
                eval(statement);
            end
            name = [Fname ,'/'  , SubCode, '_RAND' , '_B' , num2str(BN) , '.tgt'];
            dsave(name,orderfields(RndmSeq,OrderFields));
        end
        %% creat Trained blocks
        clear TrndSeq
        Window = [1:8];
        Trials = 1; % Trials per window, per sequence number which is 4*8*2 altogether
        for BN = 1:20
            TrndSeq = [];
            for e = 1:length(Window)
                for i = 1:length(Trials)
                    for seqNum = 1:size(TrainedSeqs , 1)
                        A.Window= Window(e);
                        A.seqNumb = seqNum;
                        seq = TrainedSeqs(seqNum , :);
                        A.cueP {1}= char(regexprep(cellstr(num2str(seq)),'\s',''));
                        A.iti = 1000;
                        A.hand = 2;
                        A.sounds = 1;
                        for press= 1:MaxPress
                            comnd  = [' A.press' , num2str(press) , ' = seq(press);'];
                            eval(comnd);
                        end
                        TrndSeq = addstruct(TrndSeq , A);
                        clear A
                    end
                end
            end
            idxAll = randperm(length(TrndSeq.Window));
            for kk = 1:length(OrderFields)
                statement = ['TrndSeq.' , OrderFields{kk},' = ','TrndSeq.' , OrderFields{kk},'(idxAll);'];
                eval(statement);
            end
            name = [Fname ,'/'  , SubCode, '_TRND' , '_B' , num2str(BN) , '.tgt'];
            dsave(name,orderfields(TrndSeq,OrderFields));
            clear x
        end
    case 'Memory'
        %% check if the trained sequecnes are provided, else generate them
        reGen = 'n';
        if isempty(TrainedSeqs)
            while reGen == 'n'
                TrainedSeqs(1:4,1) = [1:4]'; % make sure each sequecne starts with a different number
                for seqNum = 1:4
                    seq = [TrainedSeqs(seqNum,1) , sample_wor([1:5],1,MaxPress-1)];
                    tempSeq = diff(seq);
                    while sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                        seq = [TrainedSeqs(seqNum,1) , sample_wor([1:5],1,MaxPress-1)];
                        tempSeq = diff(seq);
                    end
                    TrainedSeqs(seqNum , 1:MaxPress) = seq;
                end
                TrainedSeqs
                reGen = input('Happy with the sequences? (y/n).... (n will regenerate sequences)', 's');
            end
        end
        Fname  = [SubCode,'_', what, '_tgtFiles'];
        mkdir(Fname)
        OrderFields = {'seqNumb'    'press1'	'press2'	'press3'	'press4'	'press5'	'press6'	'press7'	'press8'	'press9'	'cueI'	'cueP'	'cuePPost'	'iti'	'sounds'	'trialtype'};
        %% creat Memory blocks with one sequecne per block
        
        CUEI = {'Memorize_the_Sequence_and_Press_Space_to_Continue', 'Recall_and_Type_the_Sequence.'}; % instructions for trial types 0 and 1
        ITI = [100 , 8000]; %inter-trial interval for trial types 0 and 1 resectively
        Trials = 1:30; % Trials per sequence, per tirlaType (0,1 for 'memorize' and 'recall' respectively) which is 4*2*5 altogether
        
        for seqNum = 1:size(TrainedSeqs , 1)
            for BN = 1:5
                trn = 1;
                tro = 1;
                TrndSeq = [];
                trialtype= [0;ones(max(Trials)-1 , 1)];
                for i = 1:length(Trials)
                    A.seqNumb = seqNum;
                    seq = TrainedSeqs(seqNum , :);
                    if trialtype(i) == 0
                        A.cueP{1}     = char(regexprep(cellstr(num2str(seq)),'\s',''));
                    else
                        A.cueP{1}     = char(regexprep(cellstr('*********'),'\s',''));
                    end
                    A.cuePPost{1} = char(regexprep(cellstr(num2str(seq)),'\s',''));
                    A.cueI{1}= CUEI{trialtype(i)+1};
                    A.iti = ITI(trialtype(i)+1);
                    A.sounds = 1;
                    for press= 1:MaxPress
                        comnd  = [' A.press' , num2str(press) , ' = seq(press);'];
                        eval(comnd);
                    end
                    A.trNums = trn;
                    A.trOrder = tro;
                    TrndSeq = addstruct(TrndSeq , A);
                    trn = trn+1;
                    tro = tro+1;
                end
                
                TrndSeq.trialtype= trialtype;
                TrndSeq = rmfield(TrndSeq , {'trNums' , 'trOrder'});
                
                name = [Fname ,'/'  , SubCode, '_MEMO',num2str(seqNum) , '_B' , num2str(BN) , '.tgt'];
                dsave(name,orderfields(TrndSeq,OrderFields));
                clear x
            end
        end
        %% creat Memory blocks with alternate 'Memorize' and 'Recall' trials
        clear TrndSeq
        CUEI = {'Memorize_the_Sequence_and_Press_Space_to_Continue...', 'Recall_and_Type_the_Previous_Sequence.'}; % instructions for trial types 0 and 1
        ITI = [500 , 5000]; %inter-trial interval for trial types 0 and 1 resectively
        Trials = 1:5; % Trials per sequence, per tirlaType (0,1 for 'memorize' and 'recall' respectively) which is 4*2*5 altogether
        
        for BN = 1:20
            trn = 1;
            tro = 1;
            TrndSeq = [];
            for i = 1:length(Trials)
                for seqNum = 1:size(TrainedSeqs , 1)
                    for e = 0:1 % represents trial type (0,1)
                        A.trialtype= e;
                        A.seqNumb = seqNum;
                        seq = TrainedSeqs(seqNum , :);
                        if e == 0
                            A.cueP{1}     = char(regexprep(cellstr(num2str(seq)),'\s',''));
                        else
                            A.cueP{1}     = char(regexprep(cellstr('*********'),'\s',''));
                        end
                        A.cuePPost{1} = char(regexprep(cellstr(num2str(seq)),'\s',''));
                        A.cueI{1}= CUEI{e+1};
                        A.iti = ITI(e+1);
                        A.sounds = 1;
                        for press= 1:MaxPress
                            comnd  = [' A.press' , num2str(press) , ' = seq(press);'];
                            eval(comnd);
                        end
                        A.trNums = trn;
                        A.trOrder = tro;
                        TrndSeq = addstruct(TrndSeq , A);
                        trn = trn+1;
                    end
                    tro = tro+1;
                end
            end
            idxAll = randperm(tro-1);
            T = [];
            for id = 1:length(idxAll)
                T = addstruct(T , getrow(TrndSeq , TrndSeq.trOrder == idxAll(id)));
            end
            TrndSeq = rmfield(T , {'trNums' , 'trOrder'});

            name = [Fname ,'/'  , SubCode, '_MEMO_All_' , '_B' , num2str(BN) , '.tgt'];
            dsave(name,orderfields(TrndSeq,OrderFields));
            clear x
        end
        
        %% creat Memory blocks with just 'Recall' trials
        clear TrndSeq
        ITI = [100 , 5000]; %inter-trial interval for trial types 0 and 1 resectively
        
        Trials = 1:8; % Trials per sequence, per tirlaType (0,1 for 'memorize' and 'recall' respectively) which is 4*2*5 altogether
        for BN = 1:10
            trn = 1;
            tro = 1;
            TrndSeq = [];
            for i = 1:length(Trials)
                for seqNum = 1:size(TrainedSeqs , 1)
                    for e = 1 % just recall
                        A.trialtype= e;
                        A.seqNumb = seqNum;
                        seq = TrainedSeqs(seqNum , :);
                        cueItemp = [num2str(seq(1)) , '********'];
                        A.cueP{1}     = char(regexprep(cellstr(cueItemp),'\s',''));
                        A.cuePPost{1} = char(regexprep(cellstr(num2str(seq)),'\s',''));
                        A.cueI{1}= 'Type_the_Sequence_Starting_with:';
                        A.iti = ITI(e+1);
                        A.sounds = 1;
                        for press= 1:MaxPress
                            comnd  = [' A.press' , num2str(press) , ' = seq(press);'];
                            eval(comnd);
                        end
                        A.trNums = trn;
                        A.trOrder = tro;
                        TrndSeq = addstruct(TrndSeq , A);
                        trn = trn+1;
                    end
                    tro = tro+1;
                end
            end
            idxAll = randperm(trn-1);
            T = [];
            for id = 1:length(idxAll)
                T = addstruct(T , getrow(TrndSeq , TrndSeq.trOrder == idxAll(id)));
            end
            TrndSeq = rmfield(T , {'trNums' , 'trOrder'});

            name = [Fname ,'/'  , SubCode, '_MEMO_Rec_' , '_B' , num2str(BN) , '.tgt'];
            dsave(name,orderfields(TrndSeq,OrderFields));
            clear x
        end
    case 'Behav&Memo'
        %% check if the trained sequecnes are provided, else generate them
        reGen = 'n';
        if isempty(TrainedSeqs)
            while reGen == 'n'
                TrainedSeqs(1:4,1) = [1:4]'; % make sure each sequecne starts with a different number
                for seqNum = 1:4
                    seq = [TrainedSeqs(seqNum,1) , sample_wor([1:5],1,MaxPress-1)];
                    tempSeq = diff(seq);
                    while sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                        seq = [TrainedSeqs(seqNum,1) , sample_wor([1:5],1,MaxPress-1)];
                        tempSeq = diff(seq);
                    end
                    TrainedSeqs(seqNum , 1:MaxPress) = seq;
                end
                TrainedSeqs
                reGen = input('Happy with the sequences? (y/n).... (n will regenerate sequences)', 's');
            end
        end
        TrainedSeqs = seqwin_targetfile(SubCode , 'Behavioral', 'TrainedSeqs' , TrainedSeqs,'baseDir' , baseDir); 
        TrainedSeqs = seqwin_targetfile(SubCode , 'Memory' , 'TrainedSeqs' , TrainedSeqs,'baseDir' , baseDir);
end
