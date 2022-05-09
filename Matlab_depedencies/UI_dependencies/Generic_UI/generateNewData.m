function data = generateNewData(data,newDataPoint)
data.x= newDataPoint(:,1);
data.y= newDataPoint(:,2);
data.z= newDataPoint(:,3);
intensity = data.get_intensity();
data.intensity = intensity(true(size(newDataPoint(:,1),1),1));
data.bits =[];
classification = data.get_classification();
data.classification = classification(true(size(newDataPoint(:,1),1),1));
user_data = data.get_user_data();
data.user_data = user_data(true(size(newDataPoint(:,1),1),1));
scan_angle = data.get_scan_angle();
data.scan_angle = scan_angle(true(size(newDataPoint(:,1),1),1));
point_source_id = data.get_point_source_id();
data.point_source_id = point_source_id(true(size(newDataPoint(:,1),1),1));
gps_time = data.get_gps_time();
try
    data.gps_time = gps_time(true(size(newDataPoint(:,1),1),1));
catch
    ss = 0;
end
data.read_color();
data.red=[];
data.green=[];
data.blue=[];
data.nir=[];
data.extradata = [];
data.read_point_wave_info();
data.Xt=[];
data.Yt=[];
data.Zt=[];
data.wavedescriptors=[];
data.selection = true(size(newDataPoint(:,1),1),1);
data.wave_return_point=[];
data.extradata=[];

end

