function figclosekey(h, key)
%FIGCLOSEKEY Add a hotkey for closing the figure.
% Usage:
%   figclosekey(h, key)
% 
% Args:
%   h: figure handle (default: gcf)
%   key: hotkey to use (default: 'q')
% 
% See also: event, addlistener

if nargin < 2 || isempty(key); key = 'q'; end
if nargin == 1 && ischar(h); [h,key] = swap(h,key); end
if nargin < 1 || isempty(h); h = gcf(); end

set(h, 'KeyPressFcn',@(h,evt)KeyPressFcn_cb(h,evt,key));

end

function KeyPressFcn_cb(h,evt,key)
    if strcmp(evt.Key,key)
        delete(h)
    end
end