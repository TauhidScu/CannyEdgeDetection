 close all;  % Close figures
    sigma = 1.4; % Gaussian filter sigma
    highThresholdRatio = 0.275; % High threshold ratio
    lowThresholdRatio = 0.25; % Low threshold ratio OF THE high threshold
     step = 0;   
    im = imread('Test_Photos\canny.bmp');
    figure(1); imshow(im);title('Original Image');
    
    % Smooth with Gaussian 5x5 filter to reduce noise
    im = rgb2gray(im);
    figure(2); imshow(im);title('B/W Image');
    
    im = double(imgaussfilt(im,sigma));
    figure(3); imshow(im,[]);title('Gaussian Filter');
    
    % Find the intensity gradient of the image

    Gx =[-1 0 1;-2 0 2;-1 0 1];
    Gx =imfilter(im,Gx,sigma);
    
    Gy =[1 2 1;0 0 0;-1 -2 -1];
    Gy =imfilter(im,Gy,sigma);
   
    figure(4); imshow(Gx, []);title('Gx Sobel Filter');
   
    figure(5); imshow(Gy, []);title('Gy Sobel Filter');
    
    % Find the magnitude of the gradient
    Gmag = sqrt(Gx.^2 + Gy.^2);
    angle = atan2(Gy,Gx)*180/pi;
    figure(6); imshow(Gmag,[]);title('Gmag');
   % Perform non-maximum suppression using interpolation
    [h,w] = size(im);
    X=[-1,0,+1 ;-1,0,+1 ;-1,0,+1];
	Y=[-1,-1,-1 ;0,0,0 ;+1,+1,+1];
    output = zeros(h,w);
    x = [0 1];
    for i=2:h-1 % row
        for j=2:w-1 % col
            if (angle(i,j)>=0 && angle(i,j)<=45) || ...
                    (angle(i,j)<-135 && angle(i,j)>=-180)
                yBot = [Gmag(i,j+1) Gmag(i+1,j+1)];
                yTop = [Gmag(i,j-1) Gmag(i-1,j-1)];
                x_est = abs(Gy(i,j)/Gmag(i,j)); % y
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1))) % interpolation
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>45 && angle(i,j)<=90) || ...
                    (angle(i,j)<-90 && angle(i,j)>=-135)
                yBot = [Gmag(i+1,j) Gmag(i+1,j+1)];
                yTop = [Gmag(i-1,j) Gmag(i-1,j-1)];
                x_est = abs(Gx(i,j)/Gmag(i,j));
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1)))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>90 && angle(i,j)<=135) || ...
                    (angle(i,j)<-45 && angle(i,j)>=-90)
                yBot = [Gmag(i+1,j) Gmag(i+1,j-1)];
                yTop = [Gmag(i-1,j) Gmag(i-1,j+1)];
                x_est = abs(Gx(i,j)/Gmag(i,j));
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1)))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            elseif (angle(i,j)>135 && angle(i,j)<=180) || ...
                    (angle(i,j)<0 && angle(i,j)>=-45)
                yBot = [Gmag(i,j-1) Gmag(i+1,j-1)];
                yTop = [Gmag(i,j+1) Gmag(i-1,j+1)];
                x_est = abs(Gx(i,j)/Gmag(i,j));
                if (Gmag(i,j) >= ((yBot(2)-yBot(1))*x_est+yBot(1)) && ...
                    Gmag(i,j) >= ((yTop(2)-yTop(1))*x_est+yTop(1)))
                    output(i,j)= Gmag(i,j);
                else
                    output(i,j)=0;
                end
            end           
        end
    end
    
Gmag = output/max(max(output));
    figure(7); imshow(Gmag,[]);title('Non Maximum Suppression');
    % Perform double thresholding
    highThreshold = max(max(Gmag))*highThresholdRatio;
    lowThreshold = highThreshold*lowThresholdRatio;
    strongEdgesRow = zeros(1,h*w); % Keep track of the strong edge row index
    strongEdgesCol = zeros(1,h*w); % Keep track of the strong edge col index
    weakEdgesRow = zeros(1,h*w);  % Keep track of the weak edge row index
    weakEdgesCol = zeros(1,h*w);  % Keep track of the weak edge col index
    strongIndex = 1;
    weakIndex = 1;
    for i=2:h-1 % row
        for j=2:w-1 % col
            if Gmag(i,j) > highThreshold    % Strong edge
                Gmag(i,j) = 1;
                strongEdgesRow(strongIndex) = i;
                strongEdgesCol(strongIndex) = j;
                strongIndex = strongIndex + 1;
            elseif Gmag(i,j) < lowThreshold % No edge
                Gmag(i,j) = 0;
            else                            % Weak edge
                weakEdgesRow(weakIndex) = i;
                weakEdgesCol(weakIndex) = j;
                weakIndex = weakIndex + 1;
            end
        end
    end
    figure(8); imshow(Gmag,[]);title('Double Threshold'); 
    
    % Perform edge tracking by hysteresis
    set(0,'RecursionLimit',10000)
    for i=1:strongIndex-1
        % Find the weak edges that are connected to strong edges and set 
        % them to 1
        Gmag = FindConnectedWeakEdges(Gmag, strongEdgesRow(i),...
            strongEdgesCol(i));
    end
    figure(9); imshow(Gmag,[]);title('Edge Tracking Before Clean Up'); 
 % Remove the remaining weak edges that are not actually edges
    % and is noise instead
    for i=1:weakIndex-1
        if Gmag(weakEdgesRow(i),weakEdgesCol(i)) ~= 1
            Gmag(weakEdgesRow(i),weakEdgesCol(i)) = 0;
        end
    end
    figure(10); imshow(Gmag,[]);title('Edge Tracking After Clean Up'); 
    