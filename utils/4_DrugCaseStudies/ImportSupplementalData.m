function [outputs,sheets] = ImportSupplementalData(dbdir,snmdir)

% load full sheets: LV control ECM array, LV infarct ECM array, SEM for
% both
study = split(dbdir,"_");
studytitle = study(3);
sheets = sheetnames(dbdir);
for sheet = 1:length(sheets)
    fprintf("Importing database: %s %s...\n",studytitle,sheets(sheet));
    opts = detectImportOptions(dbdir,"Sheet",sheets(sheet));
    dat{sheet} = readtable(dbdir,opts);
end

% load snm: species sheet
opts = detectImportOptions(snmdir,"Sheet","species");
snm = readtable(snmdir,opts);

% filter rows for experimental groups (based on study)
if studytitle == "Ramirez"
    treatments = ["Day0" "Saline" "Valsartan"];
    for i = 1:length(dat)
        single = dat{i};
        single = single(ismember(single.Treatment,treatments),:);
        filt{i} = single;
    end
end

% rename columns for snm outputs
if studytitle == "Ramirez"
    for i = 1:length(filt)
        single = filt{i};
        single = renameCols(single,snm);
        outputs{i} = single;
    end
else
    outputs = dat;
end



% utility functions
function renamed = renameCols(db,snm)
% filters table columns for those matching list of gene names
names = [];
for gene = 2:length(db.Properties.VariableNames)
    row = find(contains(snm.geneName,db.Properties.VariableNames(gene),"IgnoreCase",true));
    if isempty(row)==0 && length(row)==1
        names = [names row];                % save snm species
    elseif isempty(row)==0 && length(row)>1
        names = [names row(end)];           % for FN1 gene
    end
end
renamed = db;
renamed.Properties.VariableNames(2:end) = snm.ID(names);   % rename to snm species



