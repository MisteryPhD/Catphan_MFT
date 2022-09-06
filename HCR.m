%% Clear all: workspace, figures and command window
clear all;
close all;
clc;

%% Read\Load the data
    % Load DICOM images
for Isl = 1:93
    data(:,:, Isl) = double(dicomread(['CBCT_Head_20190814\Head_20190814_',num2str(Isl-1),'.dcm']));
end

    % Read \ Parse METADATA - assume all read pictures have the same
    % tech. parameters (read "metadata" for one of them)
dcm_info = dicominfo(['CBCT_Head_20190814\Head_20190814_',num2str(1-1),'.dcm']);


%% Localize "resolution gauge"
% [37:54 - circles (contrast gauge?)]
% [66:68 - plates (resolution gauge?)]
%
% Assume z-location is already known: [66-68] == 67 index
% Assume central symmetry: image center AND phantom axes (center of phantom
% z-slice) are the same | all gauge is equidistant from the center (image
% central point)
%
Ig = 67;
data_slice = data(:,:,Ig);

c_x = 1 + (size(data,1)-1)/2;
c_y = 1 + (size(data,2)-1)/2;


% Take the mean and std values on the "free-of-objects" space of phantom
% (center)
FS_x_ind = (floor(c_x)-20):(ceil(c_x)+20);
FS_y_ind = (floor(c_y)-20):(ceil(c_y)+20);

c_mean = mean( reshape( data_slice(FS_x_ind,FS_y_ind, 1) ,1,[]) );
c_std  = std ( reshape( data_slice(FS_x_ind,FS_y_ind, 1) ,1,[]) );

% Look for gauge "regions" (rectangle that contain gauge) - assume that
% such region starts and ends with the signal value much greater than 
% "free-of-objects" region mean value (10 std more).


    % Scan the image area by "ray"
fig = figure;
    % Define the length of "scanning ray"
R = 100;
    % Define the initial ray direction (straight down)
angle_0 = -90;  % degrees
    % scan all the circle with discrete angle step
angle_step = atan(1/R);
angle = (0:angle_step:360) + angle_0;

endp_b_R_previous = 0;
endp_e_R_previous = 0;

gauge_index = 1;

for Iang = 1:length(angle)
    
        % Estimate the indices of the image pixels that lies on the current
        % ray 
            % Compute ray endpoints (starting on the center and ending
            % defined by its radius and angle)
    x_ray_endp = c_x + [0 R*cos((angle(Iang)/180)*pi)];
    y_ray_endp = c_y + [0 R*sin((angle(Iang)/180)*pi)];
    
            % Compute indices
    if((abs(angle(Iang)) == 90) || abs(angle(Iang)) == 270)
                % vertical line indices computed in one way (without actual
                % rays' line equation)...
                    % (keep the direction of the ray saved in indices)
        if(y_ray_endp(2)>=y_ray_endp(1))
            y_ray = floor(y_ray_endp(1)):1:ceil(y_ray_endp(2));
            x_ray = floor(x_ray_endp(1))*ones(size(y_ray));
        else
            y_ray = flip(floor(y_ray_endp(2)):1:ceil(y_ray_endp(1)));
            x_ray = floor(x_ray_endp(1))*ones(size(y_ray));
        end
    else
                % all others lines indices computed in using actual
                % rays' line equation...
        b = (y_ray_endp(2)-y_ray_endp(1))/(x_ray_endp(2)-x_ray_endp(1));
        a = y_ray_endp(1) - x_ray_endp(1)*b;
        
                   % and try to collect all points on that lies on the ray
                   % (relying on the x-unique indices where more x-unique
                   % indices appear and relying on the y-unique indices 
                   % where more y-unique indices appear.)
        if((angle(Iang)>=-45)&&(angle(Iang)<=45)||...
           (angle(Iang)>=135)&&(angle(Iang)<=225))
                       % (more x-unique indices here)
                       % (keep the direction of the ray saved in indices)
            if(x_ray_endp(2)>=x_ray_endp(1))
                x_ray = floor(x_ray_endp(1)):1:ceil(x_ray_endp(2));
                y_ray = floor(a + b*x_ray);
            else
                x_ray = flip(floor(x_ray_endp(2)):1:ceil(x_ray_endp(1)));
                y_ray = floor(a + b*x_ray);
            end
        else
                       % (and more y-unique indices here)
                       % (keep the direction of the ray saved in indices)
            if(y_ray_endp(2)>=y_ray_endp(1))
                y_ray = floor(y_ray_endp(1)):1:ceil(y_ray_endp(2));
                x_ray = floor((y_ray-a)/b);
            else
                y_ray = flip(floor(y_ray_endp(2)):1:ceil(y_ray_endp(1)));
                x_ray = floor((y_ray-a)/b);
            end
        end
    end
    
    % Get the values of the image pixels that lies on this ray
    ray_data = permute(data_slice(bsxfun(@plus,(y_ray+size(data_slice,1)*(x_ray-1)).',(size(data_slice,1)*size(data_slice,2))*(0:size(data_slice,3)-1).')),[2 1]);
    % Find: plate endpoints and the mean value over the plate (if there 
    % is a plate on the ray path)
        % Begin (b) and end (e) plate point indices 
    endp_b_i = find( (ray_data-c_mean) > (10*c_std),1, 'first');
    endp_e_i = find( (ray_data-c_mean) > (10*c_std),1, 'last');
        % If there was a plate - get the distance to the begin and end
        % (in respect to the center)
    if(~isempty(endp_b_i))
        endp_b_R(Iang) = sqrt( (x_ray(endp_b_i) - x_ray(1))^2 + (y_ray(endp_b_i) - y_ray(1))^2);
        endp_e_R(Iang) = sqrt( (x_ray(endp_e_i) - x_ray(1))^2 + (y_ray(endp_e_i) - y_ray(1))^2);
        
        endp_b_R_previous = endp_b_R(Iang);
        endp_e_R_previous = endp_e_R(Iang);
        
        mean_value(Iang) = mean(ray_data(endp_b_i:endp_e_i));
    else
        % If there was no plate - just note it by zero b_R and e_R...
        endp_b_R(Iang) = 0;
        endp_e_R(Iang) = 0;
        % ... and compute the mean value on the projection of the previous
        % "plate region", just for case (could be used later, but will not)
        endp_b_i_eq = find( (sqrt( (x_ray - x_ray(1)).^2 + (y_ray - y_ray(1)).^2) >= endp_b_R_previous),1, 'first');
        endp_e_i_eq = find( (sqrt( (x_ray - x_ray(1)).^2 + (y_ray - y_ray(1)).^2) >= endp_e_R_previous),1, 'first');
        
        mean_value(Iang) = mean(ray_data(endp_b_i_eq:endp_e_i_eq));
    end
    
% visualize the process of plates/gauge seeking - uncomment to watch   
%      subplot(2,1,1);
%       hold off;
%       imagesc(data_slice(:,:,1));
%       set(gca, 'YDir','normal');
%       axis square
%       hold on;
%       line(c_x+[0 R*cos((angle(Iang)/180)*pi)], c_y+[0 R*sin((angle(Iang)/180)*pi)],'Color','r');
%       %scatter(x_ray, y_ray);
%       subplot(2,1,2);
%       plot( ray_data);
%      
%      pause(0.1)

end

    % Use the data gathered above to finally localize gauge: find
    % rectangle (defined by two angles and two distances ~ polar 
    % coordinates) that contain gauge

        % At first, find image all gauge like regions
        % (it could be full gauge or its' pieces:
        % left_border - is an angle index of the region beginning,
        % right_border - is an angle index of the region ending
        % lengthes - is an angle index "distance" between begin and end)
bin = (endp_e_R~=0);
bin_diff = bin(2:end)-bin(1:(end-1));
left_borders = find(bin_diff==1);
right_borders = find(bin_diff==-1);
lengthes = right_borders - left_borders;

        % At second, refine the series of "gauge like regions" - grouping
        % some regions relying on the fact (in some cases it could fail,
        % BE CAREFULL) that each next gauge should occupy more space than
        % previous.
left_borders_refined = [];
right_borders_refined = [];
lengthes_refined = [];

Il_refined = 1;
left_borders_refined(Il_refined) = left_borders(1);
right_borders_refined(Il_refined) = right_borders(1);
lengthes_refined(Il_refined)= right_borders_refined(Il_refined) - left_borders_refined(Il_refined);

Il_refined = 2;
left_borders_refined(Il_refined) = left_borders(2);
right_borders_refined(Il_refined) = right_borders(2);
lengthes_refined(Il_refined)= right_borders_refined(Il_refined) - left_borders_refined(Il_refined);

    % (Group regions in the way, that each gauge will be no less than 97% 
    % of the previous == 97% or 0.97 empirical here, could not work in other cases
    %. Assume here should be added additional rule on the distances between gauges,
    % but leaving this for future)
for Il = 3:length(lengthes)
    if( lengthes_refined(Il_refined) < 0.97*lengthes_refined(Il_refined-1) ) 
        right_borders_refined(Il_refined) = right_borders(Il);
        lengthes_refined(Il_refined)= right_borders_refined(Il_refined) - left_borders_refined(Il_refined);
    else
        Il_refined = Il_refined + 1;
        left_borders_refined(Il_refined) = left_borders(Il);
        right_borders_refined(Il_refined) = right_borders(Il);
        lengthes_refined(Il_refined)= right_borders_refined(Il_refined) - left_borders_refined(Il_refined);
    end
end

    % Finally, construct a series of gauges "rectangles"
for Ig = 1:length(lengthes_refined)

    gauge{Ig}.Iang_lb =  left_borders_refined(Ig);
    gauge{Ig}.Iang_rb = right_borders_refined(Ig);

    gauge{Ig}.ang_lb =  angle(gauge{Ig}.Iang_lb);
    gauge{Ig}.ang_rb =  angle(gauge{Ig}.Iang_rb);
    
    R_temp = endp_b_R(gauge{Ig}.Iang_lb:gauge{Ig}.Iang_rb);
    gauge{Ig}.b_R = mean(R_temp(R_temp~=0));

    R_temp = endp_e_R(gauge{Ig}.Iang_lb:gauge{Ig}.Iang_rb);
    gauge{Ig}.e_R = mean(R_temp(R_temp~=0));
end

%visualize the process of plates/gauge seeking - uncomment to watch   
% fig = figure;
% imagesc(data_slice(:,:,1));
% set(gca, 'YDir','normal');
% axis square
% hold on;
% 
% for Ig = 1:length(gauge)
%     rect_x = c_x + [gauge{Ig}.b_R * cos((gauge{Ig}.ang_lb/180)*pi),...
%                     gauge{Ig}.e_R * cos((gauge{Ig}.ang_lb/180)*pi),...
%                     gauge{Ig}.e_R * cos((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * cos((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * cos((gauge{Ig}.ang_lb/180)*pi)];
%     rect_y = c_y + [gauge{Ig}.b_R * sin((gauge{Ig}.ang_lb/180)*pi),... 
%                     gauge{Ig}.e_R * sin((gauge{Ig}.ang_lb/180)*pi),...
%                     gauge{Ig}.e_R * sin((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * sin((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * sin((gauge{Ig}.ang_lb/180)*pi)];
% 
%     line(rect_x, rect_y, 'Color','r');
% 
%     text(c_x+ gauge{Ig}.b_R*0.8 * cos(( ((gauge{Ig}.ang_lb+gauge{Ig}.ang_rb)/2) /180)*pi),...
%          c_y+ gauge{Ig}.b_R*0.8 * sin(( ((gauge{Ig}.ang_lb+gauge{Ig}.ang_rb)/2) /180)*pi),...
%           num2str(Ig), 'FontSize', 6);
% end

    % Having boundaries - find all image points that lies in appropriate
    % boundaries (find indices of all points inside boundaries to compute
    % std/mean value precisely)
for Ig = 1:length(gauge)

    gauge{Ig}.points = [];
    
    alpha1 = gauge{Ig}.ang_lb;
    alpha2 = gauge{Ig}.ang_rb;
    R1 = gauge{Ig}.b_R; 
    R2 = gauge{Ig}.e_R;
    p_x = c_x + [ R1 * cos((alpha1/180)*pi),...
                  R2 * cos((alpha1/180)*pi),...
                  R2 * cos((alpha2/180)*pi),...
                  R1 * cos((alpha2/180)*pi)];
    p_y = c_y + [ R1 * sin((alpha1/180)*pi),...
                  R2 * sin((alpha1/180)*pi),...
                  R2 * sin((alpha2/180)*pi),...
                  R1 * sin((alpha2/180)*pi)];  
    
    b1 = (p_y(2)-p_y(1)) / (p_x(2)-p_x(1));
    a1 = p_y(1) - b1 * p_x(1);
       
    b2 = (p_y(3)-p_y(4)) / (p_x(3)-p_x(4));
    a2 = p_y(4) - b2 * p_x(4);
                  
    if( (alpha1 >= 0) && (alpha2 <= 90) )

        x_min = floor(min(p_x)); x_max = ceil(max(p_x));
        y_min = floor(min(p_y)); y_max = ceil(max(p_y));
        
        for x = x_min:1:x_max
            for y = y_min:1:y_max
                if((y >= a1+b1*x) && (y <= a2+b2*x)&&...
                   ( ((x-c_x)^2 + (y-c_y)^2) >= R1^2) &&... 
                   ( ((x-c_x)^2 + (y-c_y)^2) <= R2^2) )
                    gauge{Ig}.points = [gauge{Ig}.points; [x, y] ];
                end
            end
        end
        
    end
    
    if( (alpha1 >= 90) && (alpha2 <= 180) )

        x_min = floor(min(p_x)); x_max = ceil(max(p_x));
        y_min = floor(min(p_y)); y_max = ceil(max(p_y));
        
        for x = x_min:1:x_max
            for y = y_min:1:y_max
                if((y <= a1+b1*x) && (y >= a2+b2*x)&&...
                   ( ((x-c_x)^2 + (y-c_y)^2) >= R1^2) &&... 
                   ( ((x-c_x)^2 + (y-c_y)^2) <= R2^2) )
                    gauge{Ig}.points = [gauge{Ig}.points; [x, y] ];
                end
            end
        end
        
    end

    if( (alpha1 >= 180) && (alpha2 <= 270) )

        x_min = floor(min(p_x)); x_max = ceil(max(p_x));
        y_min = floor(min(p_y)); y_max = ceil(max(p_y));
        
        for x = x_min:1:x_max
            for y = y_min:1:y_max
                if((y <= a1+b1*x) && (y >= a2+b2*x)&&...
                   ( ((x-c_x)^2 + (y-c_y)^2) >= R1^2) &&... 
                   ( ((x-c_x)^2 + (y-c_y)^2) <= R2^2) )
                    gauge{Ig}.points = [gauge{Ig}.points; [x, y] ];
                end
            end
        end
        
    end
    
    if( ((alpha1 >= 270) && (alpha2 <= 360)) ||...
        ((alpha1 >= -90) && (alpha2 <= 0)))

        x_min = floor(min(p_x)); x_max = ceil(max(p_x));
        y_min = floor(min(p_y)); y_max = ceil(max(p_y));
        
        for x = x_min:1:x_max
            for y = y_min:1:y_max
                if((y >= a1+b1*x) && (y <= a2+b2*x)&&...
                   ( ((x-c_x)^2 + (y-c_y)^2) >= R1^2) &&... 
                   ( ((x-c_x)^2 + (y-c_y)^2) <= R2^2) )
                    gauge{Ig}.points = [gauge{Ig}.points; [x, y] ];
                end
            end
        end
        
    end

    
    if( ((alpha1 < 0)&&(alpha1 > -90)) && ((alpha2 > 0)&&(alpha2 < 90)) )

        x_min = floor(min(p_x)); x_max = ceil(c_x + R2);
        y_min = floor(min(p_y)); y_max = ceil(max(p_y));
        
        for x = x_min:1:x_max
            for y = y_min:1:y_max
                if((y >= a1+b1*x) && (y <= a2+b2*x)&&...
                   ( ((x-c_x)^2 + (y-c_y)^2) >= R1^2) &&... 
                   ( ((x-c_x)^2 + (y-c_y)^2) <= R2^2) )
                    gauge{Ig}.points = [gauge{Ig}.points; [x, y] ];
                end
            end
        end
        
    end

    if( ((alpha1 < 90)&&(alpha1 > 0)) && ((alpha2 > 90)&&(alpha2 < 180)) )

        x_min = floor(min(p_x)); x_max = ceil(max(p_x));
        y_min = floor(min(p_y)); y_max = ceil(c_y + R2);
        
        for x = x_min:1:x_max
            for y = y_min:1:y_max
                if((y >= a1+b1*x) && (y >= a2+b2*x)&&...
                   ( ((x-c_x)^2 + (y-c_y)^2) >= R1^2) &&... 
                   ( ((x-c_x)^2 + (y-c_y)^2) <= R2^2) )
                    gauge{Ig}.points = [gauge{Ig}.points; [x, y] ];
                end
            end
        end
        
    end
    
end

%visualize the process of plates/gauge seeking - uncomment to watch   
% fig = figure;
% imagesc(data_slice(:,:,1));
% set(gca, 'YDir','normal');
% axis square
% hold on;
% 
% for Ig = 1:length(gauge)
%     rect_x = c_x + [gauge{Ig}.b_R * cos((gauge{Ig}.ang_lb/180)*pi),...
%                     gauge{Ig}.e_R * cos((gauge{Ig}.ang_lb/180)*pi),...
%                     gauge{Ig}.e_R * cos((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * cos((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * cos((gauge{Ig}.ang_lb/180)*pi)];
%     rect_y = c_y + [gauge{Ig}.b_R * sin((gauge{Ig}.ang_lb/180)*pi),... 
%                     gauge{Ig}.e_R * sin((gauge{Ig}.ang_lb/180)*pi),...
%                     gauge{Ig}.e_R * sin((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * sin((gauge{Ig}.ang_rb/180)*pi),...
%                     gauge{Ig}.b_R * sin((gauge{Ig}.ang_lb/180)*pi)];
% 
%     line(rect_x, rect_y, 'Color','r');
% 
%     text(c_x+ gauge{Ig}.b_R*0.8 * cos(( ((gauge{Ig}.ang_lb+gauge{Ig}.ang_rb)/2) /180)*pi),...
%          c_y+ gauge{Ig}.b_R*0.8 * sin(( ((gauge{Ig}.ang_lb+gauge{Ig}.ang_rb)/2) /180)*pi),...
%           num2str(Ig), 'FontSize', 6);
% 
%    scatter(gauge{Ig}.points(:,1),gauge{Ig}.points(:,2))
% end

%% Computing MTF
    % Assume amount of plates for each gauge are known
gauge_plates_Number = [5*ones(1,10), 4*ones(1,3), 3, 2];
    % Compute gauge lenghtes: from the first to the last plate
for Ig = 1:length(gauge)
    gauge_Lenghtes(Ig) = pi*((gauge{Ig}.ang_rb - gauge{Ig}.ang_lb)/180) *...
                        ( (gauge{Ig}.b_R + gauge{Ig}.e_R)/2 )*...
                         dcm_info.PixelSpacing(1);
end
    % Compute gauge "frequency"
gauge_F = gauge_plates_Number ./ gauge_Lenghtes;

    % Finally, compute MTF    
            % Get the mean and standard deviation on the plates for the
            % largest plates gauge (the last on is the largest)     
    full_region_data = permute(data_slice(bsxfun(@plus,(gauge{end}.points(:,2)+size(data_slice,1)*(gauge{end}.points(:,1)-1)).',(size(data_slice,1)*size(data_slice,2))*(0:size(data_slice,3)-1).')),[2 1]);
                % Provide 2-class clustering (to get plates points and "background
                % points")
    [idx_larg, C_larg] = kmeans(full_region_data,2);    
                % Get the higher centroid as a plate - use it further as SD/CT of
                % the plate
    SD_plate = std(full_region_data(idx_larg==find(C_larg==max(C_larg))));
    MEAN_plate = mean(full_region_data(idx_larg==find(C_larg==max(C_larg))));
                % Use central free-of-objects (free of gauge) region std and mean
                % vaues as SD/CT for the free-region or background
    SD_background   = c_std;
    MEAN_background = c_mean;
    
    M_0 = (MEAN_plate - MEAN_background)/2;

for Ig = 1:length(gauge)    
        % Get all gauge points values
    full_region_data = permute(data_slice(bsxfun(@plus,(gauge{Ig}.points(:,2)+size(data_slice,1)*(gauge{Ig}.points(:,1)-1)).',(size(data_slice,1)*size(data_slice,2))*(0:size(data_slice,3)-1).')),[2 1]);
        % Compute Standard Deviation for all gauge region
    SD_full_region = std(full_region_data);
    
        % Provide 2-class clustering (to get plates points and "background
        % points")
    [idx, C] = kmeans(full_region_data,2);
    
        % Store the ratio between "background" points and "plate" points
        % (some region\gauges image could contain more "plates" points
        % that means that plates there are "not-recognizable") - ratio
        % could help to identify such situations  OR could be 
    Points_Ratio_2(Ig) = sum( idx==find(C==max(C))) / sum( idx==find(C==min(C))); 
    
    M_f = sqrt(SD_full_region^2 - (SD_plate^2+SD_background^2)/2);
    
    MTF_2(Ig) = (pi*sqrt(2)/4) * (M_f/M_0);
end

%% Plot the resulting MTF (normalized) along with "points" ration 
%  described above

fig = figure;
yyaxis left;
plot(gauge_F ,MTF_2 / max(MTF_2),'--o');
ylabel('MTF (2) (normalized)');

yyaxis right;
plot(gauge_F ,Points_Ratio_2,'--*');
ylabel('Plates-to-bckground points ratio');

xlabel('Spatial frequency, 1/mm');

grid on;

saveas(fig, 'MTF_2.png');
