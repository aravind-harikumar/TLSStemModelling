
function b = createlastemplate(a,inp)

% inp = rand(1000,3);

% load('matlab11a.mat');
% fn = fieldnames(a.record);

b.header = a.header;
b.header.n_point_records = size(inp,1);
b.header.n_points_by_return(1) = size(inp,1);
b.header.max_x = max(inp(:,1));
b.header.min_x = min(inp(:,1));
b.header.max_y = max(inp(:,2));
b.header.min_y = min(inp(:,2));
b.header.max_z = max(inp(:,3));
b.header.min_z = min(inp(:,3));

valss = randperm(size(a.record.x,1),size(inp,1));

b.record.x = inp(:,1);
b.record.y = inp(:,2);
b.record.z = inp(:,3);
b.record.intensity = a.record.intensity(valss);
b.record.return_number = a.record.return_number(valss);
b.record.number_of_returns = a.record.number_of_returns(valss);
b.record.scan_direction_flag = a.record.scan_direction_flag(valss);
b.record.flightline_edge_flag = a.record.flightline_edge_flag(valss);
b.record.classification = a.record.classification(valss);
b.record.classification_synthetic = a.record.classification_synthetic(valss);
b.record.classification_keypoint = a.record.classification_keypoint(valss);
b.record.classification_withheld = a.record.classification_withheld(valss);
b.record.scan_angle = a.record.scan_angle(valss);
b.record.user_data = a.record.user_data(valss);
b.record.point_source_id = a.record.point_source_id(valss);
% b.record.gps_time = a.record.gps_time(valss);

end
