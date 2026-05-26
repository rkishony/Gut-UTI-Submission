function out = safe_index(VEC, IND)
% SAFE_INDEX Indexes VEC at positions IND, returns NaN where IND == 0

  out = NaN(size(IND));
    mask = IND ~= 0;
      out(mask) = VEC(IND(mask));
end
