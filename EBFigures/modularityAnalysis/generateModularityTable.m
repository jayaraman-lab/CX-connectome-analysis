function generateModularityTable()

% GENMODULARITYTABLE computes the EB-wedge-specific modularity of ring
% neuron inputs to EPGs for different classes of ring neurons, and 
% tabulates and plots the results 


%define options for constructing modularity
erTypes     = {'ER1_a','ER1_b','ER2_a','ER2_b','ER2_c','ER2_d','ER4d','ER4m'};
weightTypes = {'ROIweight','weightRelative','ROIweight','weightRelative'};
prepost     = {'pre','pre','post','post'};
titles      = {'PRE_ROIweight', 'PRE_weightRelative', 'POST_ROIweight', 'POST_weightRelative'};
nShuffles   = 1000;
simMeas     = 'corr';

%specify which condition to plot:
% 1: PRE_ROIweight
% 2: PRE_weightRelative
% 3: POST_ROIweight
% 4: POST_weightRelative
plotInd = 2; 

%set figure specifications
figure('Color','w');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [.2, 0.1, 0.5, .9]);

pvals = nan(numel(erTypes),numel(weightTypes));
for i=1:numel(erTypes)
    for j=1:numel(weightTypes)
        
        [Q,Qshuffle,pvals(i,j),A] = computeModularity('erType',erTypes{i},'prePost',prepost{j},'weightType',weightTypes{j},'simMeas',simMeas,'nShuffles',nShuffles);
        if j==plotInd
            subplot(4,4,i);hold on;
            hist(Qshuffle);plot([Q,Q],[0,300],'-r','linewidth',2);ylim([0,300])
            if pvals(i,j)<1/nShuffles
                title([erTypes{i},': p < ',num2str(1/nShuffles)])
            else
                title([erTypes{i},': p = ',num2str(pvals(i,j))])
            end
            set(gca,'fontsize',16);
            
            subplot(4,4,i+8);
            if strcmp(simMeas,'corr')
                imagesc(2*A-1);set(gca,'ydir','normal');
                caxis([-1,1]);colormap(redblue);axis off
            else
                imagesc(A);set(gca,'ydir','normal');
                caxis([0,1]);colormap(redblue);axis off
            end
        end
    end
end

%write results to table
filename = ['ModularityPvals_',simMeas,'.csv'];
T = array2table(pvals,'VariableNames',titles,'RowNames',erTypes);
writetable(T,filename,'WriteRowNames',true)

end

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 
end
