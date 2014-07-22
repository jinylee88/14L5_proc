%COPYFIG Creates exact duplicate of figure
%
%h2 = copyfig(h1) - copies h1 to new figure
%h2 = copyfig(h1,h2) - copies h1 to h2 (clears h2 if it exists)
%h2 = copyfig(...,'PropertyName','PropertyValue') applies the new
%properties to h2
%
%Inputs: h1 - handle of figure to be copied. Defaults to the current figure
%Outputs: h2 - handle of new figure
%
%Modification History
% Created by Peter Hollender 2013/01/09

function h2 = copyfig(varargin)
if nargin == 0 || ischar(varargin{1})
    h1 = gcf;
else
    h1 = varargin{1};
end

if nargin>1 && mod(nargin,2)==0 && ~ischar(varargin{1})
    h2 = varargin{2};
    if ishandle(h2);close(h2);end
    figure(h2);
else
    h2=figure;
end

DefaultProperties = get(h1,'default');
if ~isempty(DefaultProperties)
    PropertyFieldNames = fieldnames(DefaultProperties);
    for i = 1:length(PropertyFieldNames)
        set(h2,PropertyFieldNames{i},DefaultProperties.(PropertyFieldNames{i}));
    end
end
ExcludeList = {'BeingDeleted','CurrentCharacter','CurrentAxes','CurrentObject','Children','Type'};
if strcmpi(get(h2,'WindowStyle'),'Docked')
    ExcludeList = [ExcludeList,'Position'];
end
Properties = get(h1);
PropertyFieldNames = fieldnames(Properties);
for i = 1:length(PropertyFieldNames)
    if ~any(strcmp(PropertyFieldNames{i},ExcludeList))
        set(h2,PropertyFieldNames{i},Properties.(PropertyFieldNames{i}));
        drawnow;
    end
end

if nargin > 1;
    if ischar(varargin{1})
        NewFields = varargin(1:2:end);
        NewValues = varargin(2:2:end);
    elseif ischar(varargin{2})
        NewFields = varargin(2:2:end);
        NewValues = varargin(3:2:end);
    elseif nargin > 2 && ischar(varargin{3})
        NewFields = varargin(3:2:end);
        NewValues = varargin(4:2:end);
    else
        NewFields = {};
        NewValues = {};
    end
    for i = 1:length(NewFields)
        set(h2,NewFields{i},NewValues{i})
    end
end

copyobj(get(h1,'children'),h2);
