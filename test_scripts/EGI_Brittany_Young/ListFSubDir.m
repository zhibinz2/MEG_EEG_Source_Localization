%%% Alex Chalard Jan 2019
%%% Return all Subdir path

function [ListF,Folder] = ListFSubDir(Path,Tag,varargin)

if nargin == 3 
    Filter = varargin{1};
    rawListF = subdir(Path);
    c = 1;
    ListF  = {};
    Folder = {};
    for n = 1 : numel(rawListF) 
        if contains(rawListF(n).name,Tag) && ~contains(rawListF(n).name,Filter)
            ListF{c}    = rawListF(n).name;
            Folder{c}   = rawListF(n).folder;
            c=c+1;
        end
    end
    
    
else
    
    
    rawListF = subdir(Path);
    c = 1;
    ListF  = {};
    Folder = {};
    for n = 1 : numel(rawListF)
        if contains(rawListF(n).name,Tag)
            ListF{c}    = rawListF(n).name;
            Folder{c}   = rawListF(n).folder;
            c=c+1;
        end
    end
    
end

end
