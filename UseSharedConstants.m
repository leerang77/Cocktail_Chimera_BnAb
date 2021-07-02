classdef UseSharedConstants < handle
%% This class allows access to the shared constants from anywhere
    properties (Constant) %class property
       Constants = SharedConstants; 
    end
end
