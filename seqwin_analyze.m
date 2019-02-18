function [varargout] = seqwin_analyze(what, varargin)
%% function [varargout] = seqwin_analyze(what, varargin)
% SequenceWindow experiment: analysis of behavioral data
%
% usage|example calls:
%
%                       seqwin_analyze('all_subj');                                 %pre-analysis: run the subject routine for all_subj
%                       seqwin_analyze('all_subj', {'s01'});                        %pre-analysis: run the subject routine for selected subjects
%                       [all_data] = seqwin_analyze('all_data');                    %pre-analysis: create .mat file containing data from subjects in subj
%
%                       [D] = seqwin_analyze('train');                              %group          analysis of training effects
%                       [D] = seqwin_analyze('train', 's01');                       %single subject analysis of training effects
%
%                       [D] = seqwin_analyze('seqWin');                             %group          analysis of sequence window interaction with training
%                       [D] = seqwin_analyze('seqWin', 's01');                      %single subject analysis of sequence window interaction with training
%
%                       [D] = seqwin_analyze('IPI');                                %group          analysis of inter-press intervals (IPIs)
%                       [D] = seqwin_analyze('IPI', 's01');                         %single subject analysis of inter-press intervals (IPIs)
%
% --
% gariani@uwo.ca - 2018.09.19

%% paths
pathToData = '../../../../data/SeqWindow';
pathToAnalyze = '../../../../data/SeqWindow/analyze';
if ~exist(pathToAnalyze, 'dir'); mkdir(pathToAnalyze); end % if it doesn't exist already, create analyze folder

%% globals

% subjects
% pilot: 'sl', 's99'
% incomplete:
% high-error:
subj = {'sl', 's99'};
ns = numel(subj);
subvec = zeros(1,ns);
for i = 1:ns
    if numel(subj{i})<3
        subvec(1,i) = 98;
    else
        subvec(1,i) = str2double(subj{i}(2:3));
    end
end

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;

% plot defaults
fs = 16; %default axes font size for all figures
fsm_l = 1.5; %default font size multiplier for axes labels
fsm_t = 1.75; %default title font size for all figures
lw = 3; %4; %default line width for all figures
ms = 10; %12; %default marker size for all figures
ws = 1:8; %[1,2,3,4,7,8]; %default window size subset (full range = 1:8)
maxWin = 0; %default ceiling level for window size (0 = no ceiling)

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,lightgray,green,lightgreen,black,silver,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
trsty = style.custom({orange, purple}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');

% legends
trleg = {'Random', 'Trained'};

%% types of analysis
switch (what)
    case 'all_subj' % pre-analysis: run the subject routine for all_subj
        if nargin>1; subj = varargin{1}; end
        for s = 1:numel(subj)
            fprintf(1, '\nsubject: %s\n', subj{s});
            datafilename = fullfile(pathToData, sprintf('Win_%s.dat', subj{s})); %input
            outfilename  = fullfile(pathToData, sprintf('seqwin_%s.mat', subj{s})); %output
            D.SN = []; % create SN info field
            D = dload(datafilename); %load dataset for this subj
            D.SN = ones(numel(D.TN), 1) * subvec(s); % fill in SN info field
            D.seqNum = D.seqNumb; D = rmfield(D, 'seqNumb'); % remove old seqNum field
            D.seqWin = D.Window;  D = rmfield(D, 'Window');  % remove old seqNum field
            D.day = ceil( ( (1:numel(D.TN)) / 360) )';
            D.train = double(D.seqNum>0);
            % add BN info separately per trained / untrained sequences
            D.BN_train = zeros(numel(D.TN), 1);
            c1 = 1; c2 = 1;
            for b = 1:max(D.BN)
                if any(D.BN==b & D.train==1)
                    D.BN_train(D.BN==b & D.train==1) = c1; c1 = c1 + 1;
                elseif any(D.BN==b & D.train==0)
                    D.BN_train(D.BN==b & D.train==0) = c2; c2 = c2 + 1;
                end
            end
            D = orderfields(D, [43,1,48,2,46,44,45,47,3:42]); % reorder fields
            save(outfilename, '-struct', 'D');
        end
        
    case 'all_data' % pre-analysis: create .mat file containing data from all subjects
        all_data = [];
        for s = 1:ns
            fprintf('\n%s\n\n', subj{s});
            D = load(fullfile(pathToData, sprintf('seqwin_%s.mat', subj{s}))); % load data structure for each subject
            all_data = addstruct(all_data, D); % append data structures from each subject
        end
        save( fullfile( pathToAnalyze, 'seqwin_all_data.mat'), '-struct', 'all_data'); % save all_data.mat file
        varargout = {all_data};
        
    case 'train' % analysis of training effects
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('seqwin_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'seqwin_all_data.mat'));
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Training - subj %s', varargin{1})); else; figure('Name',sprintf('Training - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % ET
        T = tapply(D, {'SN', 'BN_train', 'day', 'train'}, ...
            {D.MT, 'nanmean', 'name', 'ET'}, ...
            'subset', D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.BN_train, T.normET, 'length');
        
        subplot(2,2,1); title('Execution time'); hold on;
        plt.line([T.day T.BN_train], T.normET, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Block number'); ylabel('ET (ms)'); axis square; ylim([1300 4200]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        
        % stats
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % RT
        T = tapply(D,{'SN', 'BN_train', 'day', 'train'}, ...
            {D.pressTime0, 'nanmean', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.BN_train, T.normRT, 'length');
        
        subplot(2,2,2); title('Reaction time'); hold on;
        plt.line([T.day T.BN_train], T.normRT, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Block number'); ylabel('RT (ms)');  axis square; ylim([100 3000]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        
        % stats
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % ACC
        T = tapply(D,{'SN', 'BN_train', 'day', 'train'}, ...
            {(1-D.isError)*100, 'nanmean', 'name', 'ACC'});
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.BN_train, T.normACC, 'length');
        
        subplot(2,2,3); title('Accuracy'); hold on;
        plt.line([T.day T.BN_train], T.normACC, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Block number'); ylabel('ACC (%)');  axis square; ylim([50 100]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        
        % stats
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % PTS
        T = tapply(D,{'SN', 'BN_train', 'day', 'train'}, ...
            {D.points, 'nanmean', 'name', 'PTS'});
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'PTS'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.BN_train, T.normPTS, 'length');
        
        subplot(2,2,4); title('Points'); hold on;
        plt.line([T.day T.BN_train], T.normPTS, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Block number'); ylabel('PTS');  axis square; ylim([1 8]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'seqWin' % analysis of sequence window interaction with training
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('seqwin_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'seqwin_all_data.mat'));
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Learning - subj %s', varargin{1})); else; figure('Name',sprintf('Learning - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % select window sizes of interest
        D = getrow(D, ismember(D.seqWin, ws));
        
        % put a ceiling to window size
        if maxWin > 0; D.seqWin(D.seqWin >= maxWin) = maxWin; end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % ET
        T = tapply(D, {'SN', 'day', 'train', 'seqWin'}, ...
            {D.MT, 'nanmean', 'name', 'ET'}, ...
            'subset', D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ET'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.seqWin, T.normET, 'length');
        
        subplot(2,2,1); title('Execution time'); hold on;
        plt.line([T.day T.seqWin], T.normET, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Window size'); ylabel('ET (ms)'); axis square; ylim([1500 5000]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        for l = 1:numel(ax.XTickLabel); if strcmp(ax.XTickLabel(l), num2str(maxWin)); ax.XTickLabel{l} = sprintf('%s+', ax.XTickLabel{l}); end; end
        
        % stats
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % RT
        T = tapply(D, {'SN', 'day', 'train', 'seqWin'}, ...
            {D.pressTime0, 'nanmean', 'name', 'RT'}, ...
            'subset', D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'RT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.seqWin, T.normRT, 'length');
        
        subplot(2,2,2); title('Reaction time'); hold on;
        plt.line([T.day T.seqWin], T.normRT, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Window size'); ylabel('RT (ms)'); axis square; ylim([100 3000]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        for l = 1:numel(ax.XTickLabel); if strcmp(ax.XTickLabel(l), num2str(maxWin)); ax.XTickLabel{l} = sprintf('%s+', ax.XTickLabel{l}); end; end
        
        % stats
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % ACC
        T = tapply(D,{'SN', 'day', 'train', 'seqWin'}, ...
            {(1-D.isError)*100, 'nanmean', 'name', 'ACC'});
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.seqWin, T.normACC, 'length');
        
        subplot(2,2,3); title('Accuracy'); hold on;
        plt.line([T.day T.seqWin], T.normACC, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Window size'); ylabel('ACC (%)');  axis square; ylim([50 100]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        for l = 1:numel(ax.XTickLabel); if strcmp(ax.XTickLabel(l), num2str(maxWin)); ax.XTickLabel{l} = sprintf('%s+', ax.XTickLabel{l}); end; end
        
        % stats
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % PTS
        T = tapply(D,{'SN', 'day', 'train', 'seqWin'}, ...
            {D.points, 'nanmean', 'name', 'PTS'});
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'PTS'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], T.seqWin, T.normPTS, 'length');
        
        subplot(2,2,4); title('Points'); hold on;
        plt.line([T.day T.seqWin], T.normPTS, 'split',T.train, 'style',trsty, 'leg',trleg);
        xlabel('Window size'); ylabel('PTS');  axis square; ylim([1 8]);
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        for l = 1:numel(ax.XTickLabel); if strcmp(ax.XTickLabel(l), num2str(maxWin)); ax.XTickLabel{l} = sprintf('%s+', ax.XTickLabel{l}); end; end
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'IPI' % analysis of inter-press intervals (IPIs)
        if nargin>1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('seqwin_%s.mat', subj)));
        else % load group data
            D = load( fullfile(pathToAnalyze, 'seqwin_all_data.mat'));
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Learning - subj %s', varargin{1})); else; figure('Name',sprintf('Learning - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % select window sizes of interest
        D = getrow(D, ismember(D.seqWin, ws));
        
        % put a ceiling to window size
        if maxWin > 0; D.seqWin(D.seqWin >= maxWin) = maxWin; end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % IPI
        D.IPI=diff([D.pressTime0, D.pressTime1, D.pressTime2, D.pressTime3, D.pressTime4, D.pressTime5, D.pressTime6, D.pressTime7, D.pressTime8], 1, 2);
        D.IPI_1 = D.IPI(:,1); D.IPI_2 = D.IPI(:,2); D.IPI_3 = D.IPI(:,3); D.IPI_4 = D.IPI(:,4); D.IPI_5 = D.IPI(:,5); D.IPI_6 = D.IPI(:,6); D.IPI_7 = D.IPI(:,7); D.IPI_8 = D.IPI(:,8);
        
        % create summary table for IPIs
        T = tapply(D, {'SN', 'day', 'train', 'seqWin'}, ...
            {D.IPI_1,'nanmean','name','IPI1'}, ...
            {D.IPI_2,'nanmean','name','IPI2'}, ...
            {D.IPI_3,'nanmean','name','IPI3'}, ...
            {D.IPI_4,'nanmean','name','IPI4'}, ...
            {D.IPI_5,'nanmean','name','IPI5'}, ...
            {D.IPI_6,'nanmean','name','IPI6'}, ...
            {D.IPI_7,'nanmean','name','IPI7'}, ...
            {D.IPI_8,'nanmean','name','IPI8'}, ...
            'subset', D.isError==0);
        for i = 1:size(D.IPI, 2)
            T.IPI(:, i) = eval(sprintf('T.IPI%d', i));
            T = rmfield(T, sprintf('IPI%d', i));
            T.IPInum(:, i) = ones( size(T.SN, 1), 1) * i;
            T.SN(:, i) = T.SN(:, 1);
            T.day(:, i) = T.day(:, 1);
            T.train(:, i) = T.train(:, 1);
            T.seqWin(:, i) = T.seqWin(:, 1);
        end
        T.IPI = reshape(T.IPI, size(T.IPI,1) * size(T.IPI,2), 1);
        T.IPInum = reshape(T.IPInum, size(T.IPInum,1) * size(T.IPInum,2), 1);
        T.SN = reshape(T.SN, size(T.SN,1) * size(T.SN,2), 1);
        T.day = reshape(T.day, size(T.day,1) * size(T.day,2), 1);
        T.train = reshape(T.train, size(T.train,1) * size(T.train,2), 1);
        T.seqWin = reshape(T.seqWin, size(T.seqWin,1) * size(T.seqWin,2), 1);        
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.day T.train], [T.seqWin T.IPInum], T.normIPI, 'length');
        
        plt.line([T.day T.IPInum], T.normIPI, 'split',T.train, 'style',trsty, 'leg',trleg, 'subset',T.seqWin==8 & T.train==0);
        xlabel('IPI number'); ylabel('IPI duration (ms)'); ylim([0 800]); axis square;
        ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        
        %         subplot(2,2,1); title('Day 1'); hold on;
        %         plt.line([T.seqWin T.IPInum], T.normIPI, 'split',T.train, 'style',trsty, 'leg',trleg, 'subset',T.day==1);
        %         xlabel('IPI number'); ylabel('IPI duration (ms)'); ylim([0 800]); %axis square;
        %         ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        %
        %         subplot(2,2,2); title('Day 2'); hold on;
        %         plt.line([T.seqWin T.IPInum], T.normIPI, 'split',T.train, 'style',trsty, 'leg',trleg, 'subset',T.day==2);
        %         xlabel('IPI number'); ylabel('IPI duration (ms)'); ylim([0 800]); %axis square;
        %         ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        %
        %         subplot(2,2,3); title('Day 3'); hold on;
        %         plt.line([T.seqWin T.IPInum], T.normIPI, 'split',T.train, 'style',trsty, 'leg',trleg, 'subset',T.day==3);
        %         xlabel('IPI number'); ylabel('IPI duration (ms)'); ylim([0 800]); %axis square;
        %         ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        %
        %         subplot(2,2,4); title('Day 4'); hold on;
        %         plt.line([T.seqWin T.IPInum], T.normIPI, 'split',T.train, 'style',trsty, 'leg',trleg, 'subset',T.day==4);
        %         xlabel('IPI number'); ylabel('IPI duration (ms)'); ylim([0 800]); %axis square;
        %         ax = gca; ax.FontSize = fs; ax.LabelFontSizeMultiplier = fsm_l; ax.TitleFontSizeMultiplier = fsm_t;
        
        % stats
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    otherwise
        error('no such case!')
end

end