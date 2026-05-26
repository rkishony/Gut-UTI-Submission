clear
close all
clc

% paths
expname = 'IsolateResB';
folders = get_folders();
data_fldr   = fullfile(folders.source_data, 'Phenotyping', expname);
output_path = folders.figures;
examples_filename = fullfile(folders.metadata, 'example_isolate_nums.csv');

% get example choice
chosen_examples = readcell(examples_filename);
fig_name = 'fig_example_MIC_heatmap_UTI_FOBT_isolates';

% load well labels
sheetname = 'SourcePlate';
fname = fullfile(data_fldr, 'PlateLayout.xlsx');
labels =  readcell(fname, 'Sheet', sheetname);

% load MIC data for all antibiotics
distinguishing_conc = readtable(fullfile(folders.metadata, 'distinguishing_conc.xlsx'));
Abs = distinguishing_conc.Properties.VariableNames; 
breakpoints = distinguishing_conc{1,:};
breakpoint_conc_order=[3 4 1 1 3];
clear distinguishing_conc
% sheetname = 'MICs'; %%% TG
for a = 1:numel(Abs)
    filename = 'results'; %Abs{a};
    fname = fullfile(data_fldr, [filename '.xlsx']);
    % meas = readmatrix(fname, 'Sheet', sheetname); %%% TG
    meas = readmatrix(fname); %%% TG
    measurements(:,a) = meas(:);
end

chosen_inds = zeros(size(chosen_examples));
for i = 1:length(chosen_examples)
    chosen_inds(i) = find(arrayfun(@(c) strcmp(labels{c}, chosen_examples{i}), 1:numel(labels)));
end

%% plot

figure(1)
clf
set(gcf,'color','white')
%cmap=gray(256);%colormap('gray');
%cmap=redbluecmap;
cmap=turbo;
cmap=cmap*0.9;
%cmap=cmap(end:-1:80,:);
%cmap=cmap([1 2 3 7:end-1],:);
colormap(cmap);
toplot=(measurements(chosen_inds,:)-1)-repmat(breakpoint_conc_order,[length(chosen_inds) 1]);
toplot(measurements(chosen_inds,:)<2)=-3;
toplot(measurements(chosen_inds,:)>6)=3;
%imagesc(measurements(chosen_inds,:))
imagesc(toplot)
set(gca,'Units','pixels')
set(gcf,'Units','pixels')
W = 50;
nisolates = length(chosen_inds);
nabs = length(Abs);
set(gcf,'position',1.4*[3*W W (nabs+2)*W nisolates*W])
set(gca,'innerposition',[2.5*W W*0.75 nabs*W nisolates*W],'fontsize',6)
set(gca,'TickLabelInterpreter','none')
yticksnow=1:nisolates;
yticklabels(labels(chosen_inds))
xticksnow=1:length(Abs);
xticks(xticksnow);
xticklabels(Abs)
c=colorbar;
c.Limits=[-3.5 3.5];
c.Position=[0.75 0.38 0.03 0.3];
clim([-3.5 3.5]);
ticksnow=[-3:3];
for i=1:nisolates
    line([xticksnow(1)-0.5 xticksnow(end)+0.5], [yticksnow(i)+0.5 yticksnow(i)+0.5],'color','black' )
end
for i=1:length(Abs)
    line([xticksnow(i)+0.5 xticksnow(i)+0.5], [yticksnow(1)-0.5 yticksnow(end)+0.5],'color','black' )
end
ticklabelsnow=arrayfun(@(x) sprintf('2^{%i}',x),ticksnow,'UniformOutput',false);
ticklabelsnow{1}='Below detection';
ticklabelsnow{4}='Breakpoint';
ticklabelsnow{end}='Above detection';
set(c,'Ticks',ticksnow,'TickLabels',ticklabelsnow,'fontsize',6);
pbaspect([1 2.5 1])
%set(c,'Ticks',[-1:4],'TickLabels',{'Sensitive', 'Breakpoint','2XBreakpoint','4XBreakpoint', '8XBreakpoint'})

set(gcf,'PaperPosition',[2 2 3 3])
print_figure(gcf, fullfile(folders.figures, fig_name))
