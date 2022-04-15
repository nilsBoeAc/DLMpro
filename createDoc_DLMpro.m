
% automatically creates doc using m2docgen
%
% create options struct for m2doc
%% Options:
%   toolboxName - string : "Name_of_the_toolbox"
%       Distinct name that will be shown in the documentation
%   delOld - boolean: true
%       If documentation folder opts.outputFolder already exist, delete it 
%       first.
%   mFolder - string array : ["absolute_path_to_scripts"]
%       The path specified in this variable (and subfolders) will be 
%       searched for .m and .mlx files to convert to html. Multiple folders
%       are possible. 
%   outputFolder - string array : ["absolute_path_to_output_folder"]
%       The path specified in this variable will contain the converted html 
%       files, a subfolder with the .css files, and the toc xml file. 
%   excludeFolder - string array : ["folder_names_to_exclude"]
%       If the path of an m file contains these words, they will be ignored
%       and not be converted to html.
%   excludeFile - string array : ["file_name_to_exclude"]
%       If the file name contains one of these words, they will be ignored
%       and not converted to html.
%   htmlMetaFolder - string : ["relative_folder_name"]
%       This folder will contain the css-files and images of the html files
%       and will be a subfolder of opts.outputFolder.
%   htmlTemplate - string:  ["relative_folder_name"]
%       Define the folder containing the html template files that will
%       define the structure and look of the exported documents. If the
%       path is not valid, then the default template is used.
%   startpage - string array: ["name_of_landing_page_html_file_name"]
%       The very first toc-element will be displayed when opening the html
%       documentation. Specify instead here a html site as landing page.
%       Create this page by writing an m/mlx-file with the specified name
%       The file will be converted to HTML and used as landing page.
%   toc - cell: 
%       The html documentation requires an xml file (helptoc.xml) that
%       structures the documentation. If this variable is empty, then the
%       original folder structure from opts.mFolder will be used.
%       Alternatively, a custom structure can be defined:
%       First cell column:  Names displayed in toc
%       Second cell column: Folder of origin
%       Third cell column:  cell that can define a substructure
%       Example: opts.toc = {"MyToolbox", "/", {}};
%           - All files from the root directory will be inside "MyToolbox"
%       opts.toc{1,3} = {"Vehicles", ["cars" "rockets"], {}};
%           - All files whose last folder is either "cars" or "rocket" will
%           be found under a new sub-toc element instead of the root dir:
%           Mytoolbox->Vehicles
%   verbose - boolean: false
%       If true, then more intermediate steps will be documented in the
%       command window.

clearvars;

cd(fileparts(which(mfilename)));
cPath = fileparts(which(mfilename));

mF = fullfile(cPath);
oF = mF+"\doc";
htmlF = "m2doc-standard";

opts = struct(  'toolboxName',      "DLMpro", ...
                'delOld',           true, ...
                'mFolder',          mF, ...
                'outputFolder',     oF, ...
                'excludeFolder',	[], ...
                'excludeFile',      ["DLMPro_install","DLMPro_uninstall","createDoc_DLMpro"], ...
                'htmlMetaFolder',   "ressources", ...
                'htmlTemplate',     htmlF, ...
                'startPage',        "DLMPro", ...
                'toc',              [], ...
                'verbose',          false);
             
opts.toc = {"DLMpro","DLMPro",{},...
            "Examples","Exampels",{}};

% make sure to have added m2doc to the matlab path
res = m2docgen(opts);