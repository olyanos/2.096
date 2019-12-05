% System
n      = Params.n;
E      = Params.E;
A      = Params.A;
Lambda = Params.Lambda;
Gamma  = Params.Gamma;
y      = Params.y;
P      = Params.P;
d      = Params.d;
yw     = Params.yw;
gam1   = Params.gam1;
gam2   = Params.gam2;

%% Plot Image
load('Venezuela_sim_data')
load('Color_Warm')
mymap = Color_Warm2;
%mymap = colormap(hsv);
r_x  = r.x_mat;
r_Yt = r.Yt;

for hh = 1
    
    % Call data
    x   = r_x(:,hh);
    Ytv = r_Yt(:,hh);
    
    % Update figure
    clf
    imshow('world_map.jpg')
    axis off
    hold on
    
    % Load nice color map
    ncs   = length(mymap);
    cv    = linspace(0.95*min(Ytv),1.05*max(Ytv),ncs);
    
    % Define coordinates for Australia, Brazil, Canada, EU, Ghana,
    % Israel, Malysia, Singapore, Switzerland, Turkey
    
    % Put Switzerland in the Ocean, lift EU, etc
    xy_vec = [840 480;
        270 420;
        150 200;
        480 160;
        418 365;
        555 305;
        810 380;
        750 385;
        350 230;
        525 255];
    
    % Plot
    xup   = sqrt(0:0.02:1);
    xdown = (0:0.02:1).^2;
    kk    = 1;
    for ii = 1:size(xy_vec,1)
        % Plot flows
        for jj = 1:n
            if jj == ii
                % skip
            else
                if ii < jj
                    
                    % Set x
                    if xy_vec(ii,1) < xy_vec(jj,1)
                        x1 = xy_vec(ii,1);
                        x2 = xy_vec(jj,1);
                        y1 = xy_vec(ii,2);
                        y2 = xy_vec(jj,2);
                    else
                        x2 = xy_vec(ii,1);
                        x1 = xy_vec(jj,1);
                        y2 = xy_vec(ii,2);
                        y1 = xy_vec(jj,2);
                    end
                    
                    % Define quadratic function parameters
                    c1 = y1;
                    c2 = x1;
                    c3 = (y2 - c1)/((x2-c2)^2);
                    xvals = x1:0.1:x2;
                    yvals = c3*(xvals - c2).^2+c1;
                    
                else
                    % Set x
                    if xy_vec(ii,1) < xy_vec(jj,1)
                        x1 = xy_vec(ii,1);
                        x2 = xy_vec(jj,1);
                        y1 = xy_vec(ii,2);
                        y2 = xy_vec(jj,2);
                    else
                        x2 = xy_vec(ii,1);
                        x1 = xy_vec(jj,1);
                        y2 = xy_vec(ii,2);
                        y1 = xy_vec(jj,2);
                    end
                    
                    % Define sqrt function parameters
                    c1 = y1;
                    c2 = x1;
                    c3 = (y2 - c1)/(sqrt(x2-c2));
                    xvals = x1:0.1:x2;
                    yvals = c3*sqrt(xvals - c2)+c1;
                    
                end
                flow_width = min(x(kk),20);
                flow_width = max(flow_width,0.2);
                plot(xvals,yvals,'linewidth',5*flow_width);
                kk = kk+1;
            end
        end
        
        % Make cirlce size proportional to GDP
        GDP_size = 20*log10(20*y(ii));
        [~,ind] = min(abs(Ytv(ii) - cv));
        plot(xy_vec(ii,1),xy_vec(ii,2),'o','linewidth',2,'markersize',GDP_size,'color','black','MarkerFaceColor',mymap(ind,:));
    end
    
    % Add a colorbar with tick labels
    colormap(mymap)
    drawnow
    tick = [0, 0.25, 0.5, 0.75 1];
    C  = colorbar('Location', 'east','Ticks',tick,'YTickLabel',{'$0.72','$1.37','$2.03','$2.69','$3.35'},'AxisLocation','out');
%     Frames(hh-889) = getframe(gcf);
%     hh
end

% %% Play video
% 
%     v = VideoWriter('Parameter_Convergence_T1.avi');
%     v.FrameRate = 10;
%     open(v)
%     writeVideo(v,Frames)
%     close(v)