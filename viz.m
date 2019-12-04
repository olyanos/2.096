% clean up
clear all;clc; close all;
% load the data
file_name = 'Venezuela_sim_data'
load(file_name);
Currency = r.Yt;


%% plot only the Venezula case
% figure(1);
% for i = 1:length(Currency-1)
%     
%     line=linspace(Currency(1,i),Currency(1,i+1),2);
%     timeline = linspace(i,i+1,2);
%     plot(timeline,line,'LineWidth',3,'Color','red');
%     xlim([0  length(Currency-1)]);
%     ylim([0 5]);
%     hold on;
%     pause(1e-15);
% end
% hold off;
%% plot all countries
figure(2);

% Initialize video
myVideo = VideoWriter(file_name); %open video file
myVideo.FrameRate = 50;  %can adjust this, 5 - 10 works well for me
open(myVideo)

shift = 60; % we want to plot every 60 points together to speed up the video
cmap = jet(10);
for i = 1:shift:length(Currency)-shift*2
    for j = 1:10
        line=linspace(Currency(j,i),Currency(j,i+shift),shift);
        timeline = linspace(i,i+shift,shift);
        if j==1
            plot(timeline,line,'o','Color',cmap(j, :)); 
        else
            plot(timeline,line,'LineWidth',2,'Color',cmap(j, :));
        end
        hold on;
        xlim([0  length(Currency-1)]);
        ylim([0 4]);
        set(gca,'FontName','Times','FontSize',13)
        xlabel('Time (hours)')
        ylabel('Currency ($)')
        legend('Australia', 'Brazil', 'Canada', 'EU', 'Ghana', 'Israel', 'Malysia', 'Singapore', 'Switzerland', 'Turkey')
    end
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(1e-15);
end
close(myVideo)