function [] = exportFigPaper()

% Export all open figures to eps to prep for editing in Illustrator

for i = 1:length(get(0,'children'))
    filename = sprintf('fig%i',i);
    export_fig filename -eps -transparent -preserve_size
end