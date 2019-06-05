function spyc_grid(sA,cmap, tick_minor, tick_major);

%SPYC Visualize sparsity pattern with color-coded scale.
%   SPYC(S) plots the color-coded sparsity pattern of the matrix S.
%
%   SPYC(S,CMAP) plots the sparsity pattern of the matrix S USING 
%                    COLORMAP CMAP.
%
%   SPYC(S,CMAP,PB) allows turning off the display of a colorbar by passing
%                   flag PB=0
%
%   written by Try Hard
%   $Revision: 0.0.0.2 $  $Date: 2013/08/24 11:11:11 $

if nargin<1 || nargin>4 && ~isempty(cmap)
    error( 'spyc:InvalidNumArg', 'spyc_grid takes one to five inputs')
end

if isempty(sA)
    error( 'spyc:InvalidArg', 'sparse matrix is empty')
end

if nargin>1 && ~isempty(cmap)
    % colorspy does not check whether your colormap is valid!
    if ~isnumeric(cmap)
        cmap=colormap(cmap);        
    end
else
    cmap=flipud(colormap('autumn'));
end

pb=1;
    
indx=find(sA);
[Nx, Ny]=size(sA);
sA=full(sA(indx));
ns = length(indx);
[ix, iy]=ind2sub([Nx Ny],indx);

imap = round((sA-min(sA))*63/(max(sA)-min(sA)))+1;

hold on
colormap(cmap)
scatter(iy,ix,[],imap,'Marker','.','SizeData',200)
set(gca,'ydir','reverse')
axis equal;
xlabel(['nz = ' num2str(ns)])
axis([0.5 Nx+0.5 0.5 Ny+0.5])
hold on
minors = linspace(0, Nx, Nx/tick_minor+1);
minors = minors + 0.5;
majors = linspace(0, Nx, Nx/tick_major+1);
majors = majors + 0.5;
for minorTick = minors
    yline(minorTick,'Color', [0.5 0.5 0.5], 'LineWidth',1);
    xline(minorTick,'Color', [0.5 0.5 0.5], 'LineWidth',1);
end
for majorTick = majors
    yline(majorTick,'g', 'LineWidth',3);
    xline(majorTick,'g', 'LineWidth',3);
end
hold off
box on

if pb
    colorbar
end



