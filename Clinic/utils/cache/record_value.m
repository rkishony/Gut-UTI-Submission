function out = record_value(description, label, val, fmt, verbose)
persistent RECORDED_NUMS

if nargin==1
    RECORDED_NUMS = description;
end

if isempty(RECORDED_NUMS)
    RECORDED_NUMS = table('Size', [0, 4], ...
        'VariableTypes', {'string', 'string', 'string', 'double'}, ...
        'VariableNames', {'Description', 'Label', 'sValue', 'Value'});
end

if nargin>1
    if isempty(label)
        label = inputname(3);
    end
    if nargin<4
        fmt = [];
    end
    if ischar(fmt)
        fmt = string(fmt);
    end
    if nargin<5
        verbose = true;
    end

    is_val_str = ischar(val) || isstring(val);
    if is_val_str
        val = string(val);
    end

    if ~isscalar(val) && ~is_val_str
        for i = 1:numel(val)
            v = val(i);
            snum = sprintf('%d', i);
            if isscalar(fmt) || isempty(fmt)
                f = fmt;
            else
                f = fmt(i);
            end
            record_value(char(description), [char(label) '(' snum ')'], v, f, verbose)
        end
    elseif isstruct(val)
        flds = fields(val);
        for i = 1:numel(flds)
            fld = flds{i};
            v = val.(fld);
            if ~isstruct(fmt)
                f = fmt;
            else
                fmt_flds = fields(fmt);
                if ismember(fld, fmt_flds)
                    f = fmt.(fld);
                else
                    f = [];
                end
            end
            record_value(char(description), [char(label), '.' fld], v, f, verbose)
        end
    else
        description = string(description);
        label = string(label);
        if label~="N/A" && any(label==RECORDED_NUMS.Label)
            % error('Label %s already recorded.', label)
        end
        sval = format_value(fmt, val);
        if is_val_str
            val = nan;
        end
        RECORDED_NUMS(end + 1, :) = {description, label, string(sval), val};
        if verbose
            fprintf('*** %-35s %45s: %7s\n', extractBefore(description,min(strlength(description)+1,35)), ...
                label, sval)
        end
    end
end
if nargout
    out = RECORDED_NUMS;
end
end

function out = round_to_sig_digits(val, n)
if val == 0
    out = 0;
else
    out = round(val, n - floor(log10(abs(val))) - 1);
end
end

function sval = format_value(fmt, val)
if ischar(val) || isstring(val)
    sval = string(val);
    return
end

if ~isempty(fmt)
    sval = sprintf(fmt, val);
elseif val==fix(val) || isinf(val)
    sval = num2str(val);
else
    val = round_to_sig_digits(val, 2);
    sval = sprintf('%g', val);
end

sval = regexprep(sval, 'e([+-])0*(\d+)', 'e$1$2');  %1e06 -> 1e6

end
