
% installing DLMpro
%
%% Description
% The Toolbox will be installed by adding the folder to the matlab search path
% 
% Hint: you can fast install the toolbox by marking this file in the MATLAB
% "Current Folder Browser" and pressing F9. So you do not need to open the
% script. 
%
%% Disclaimer
%   Copyright (c) 2021 Nils Böhnisch, Marc Bangel.
%
%   This file is part of DLMPro.
%
%   DLMPro is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version. Also see the file "License"

cd(fileparts(which(mfilename)));
addpath(genpath('.'));
disp("DLMpro is Installed / Added!");
% eof;