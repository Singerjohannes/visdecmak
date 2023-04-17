% output = transres_other_average(decoding_out, varargin)
% 
% This function averages other output that was generated across decoding
% steps (in the field decoding_out.opt), but works only if this other
% output is numerical
%
% To use it, use
%
%   cfg.results.output = {'other_average'}
%
% Martin, 2017-06-30
function output = transres_mean_decision_values(decoding_out, varargin)

currdim = ndims(decoding_out(1).decision_values)+1;

% average across the optional output
output = {mean(abs(cat(currdim, decoding_out.decision_values)),currdim)};