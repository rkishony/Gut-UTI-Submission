function [num_puc, num_measured, num_res] = count_UTI_resistances(...
    v123count, DayEdges, time_frame, time_method)
% time_method: 'all', 'first', 'first_measured'

iT1 = find(time_frame(1)-0.5 == DayEdges);
iT2 = find(time_frame(2)+0.5 == DayEdges) - 1;
v123count = v123count(:,:,iT1:iT2,:,:);
switch time_method
    case 'all'
        num_puc = sum(v123count(:,:,:,:),3:4);
        num_measured = sum(v123count(:,:,:,2:3),3:4);
        num_res = sum(v123count(:,:,:,3),3:4);
    case {'first', 'first_measured'}
        V4 = convert_v123count_to_v4(v123count, [], time_method, 3);
        V4 = V4(:,1);  % Looking forward from the first timebin
        num_puc = double(V4~=1);
        num_measured = double(V4==3 | V4==4);
        num_res = double(V4==4);
    otherwise
        error('Unknown time-method %s', time_method)
end

end