function isMeasured = is_measured(T)

isMeasured = ~T.from_lawn & ~T.is_neg_ctrl & ~T.rep & T.FOBT_ID>0;

end