%% PARAMETERS 

AFP=0.2;% parameter of the alpha function (for g1)
GainFB=0.25 ; % parameter for the gain of the feedback
NL=4;% power law to account for the non-linearity (inspired from estimation going from vsd to spike) 

x=0:0.001:0.5; % temporal axis

a1=0.024; % time constant of the first area
d_tc=0.004; % increase of time constant from area to areas

colcol=parula(5);

% ---------- INTEGRATIVE TIME CONSTANTS IN ALL AREAS

tc_a=a1:d_tc:(a1+5*d_tc);

% g1 is the first input (alpha function) to the first area 
g1= x.^AFP .* exp (-x./tc_a(1)) ./ tc_a(1); % alpha function of the first area (FF driven)
g1=g1./sum(g1); % normalization

clear TC
for i=1:5
    gtemp= exp (-x./tc_a(i)) ./ tc_a(i); % integration time constant of the i-st area
    TC(i,:)=gtemp./sum(gtemp);  % normalization
end

% ---------- DELAYS BETWEEN ALL AREAS

clear del
del(1,2)=0.002; % delays from area 1 to area 2
del(2,1)=0.004; % delays from area 2 to area 1

for i=1:5
    del(i+1, i+2)=del(1,2); % delays from area i+1 to area i+2
    del(i+2, i+1)=del(2,1); % delays from area i+2 to area i+1
end



%% FIGURE 2A
% In this figure, we take a simple response (yy) and add to it the effect
% of delay yy_DL, then integration time constant (yy_DL_TC) and then
% nonlinearity (yy_DL_TC_NL)

figure(21)
clf, set(gcf,'PaperPosition',[0.25 0.25 50 10]); 

  % CortConv is a function simulating the cortical integration (see below)
  % the first step is thus to integrate an alpha function (g1)
  % depending on the parameter of the CortConv function :
  % (input, time constant, delay, time vector, power law exponent)

yy           = CortConv( g1 ,squeeze(TC(1,:)),0        ,x,NL);
yy_DL        = CortConv( yy ,zeros(1,501)    , del(1,2),x,1 );
yy_DL_TC     = CortConv( yy ,squeeze(TC(1,:)), del(1,2),x,1 );
yy_DL_TC_NL  = CortConv( yy ,squeeze(TC(1,:)), del(1,2),x,NL);

subplot(1,5,1) % ---------- plot of yy
plot(x,yy,'color',colcol(1,:),'linewidth',4);
hold on
[HWg(1)]=find(yy>=max(yy/2),1);
plot(x(HWg(1)),yy(HWg(1)),'o','color',colcol(1,:),'linewidth',2,'MarkerSize',15,'MarkerFaceColor','w');
line([ 0 0],[0 0.005],'color',colcol(1,:),'linewidth',2)

xlim([-0.005 0.10])
ylim([-0.001 0.025])
xlabel('Time ms')
ylabel('Activity (au)')
grid;
Figure_appearance


subplot(1,5,2) % ---------- plot of y_DL
plot(x,yy_DL,'color',colcol(2,:),'linewidth',4);
hold on
[HWg(2)]=find(yy_DL>=max(yy_DL/2),1);
plot(x(HWg(2)),yy_DL(HWg(2)),'o','color',colcol(2,:),'linewidth',2,'MarkerSize',15,'MarkerFaceColor','w');
line([ (del(1,2)) (del(1,2))],[0 0.005],'color',colcol(1,:),'linewidth',2)

xlim([-0.005 0.10])
ylim([-0.001 0.025])
title('+ D')
grid;
Figure_appearance

subplot(1,5,3) % ---------- plot of yy_DL_TC
plot(x,yy_DL_TC,'color',colcol(2,:),'linewidth',4);
hold on
[HWg(3)]=find(yy_DL_TC>=max(yy_DL_TC/2),1);
plot(x(HWg(3)),yy_DL_TC(HWg(3)),'o','color',colcol(2,:),'linewidth',2,'MarkerSize',15,'MarkerFaceColor','w');
line([ (del(1,2)) (del(1,2))],[0 0.005],'color',colcol(1,:),'linewidth',2)
xlim([-0.005 0.10])
ylim([-0.001 0.025])
title('+ D+TC')
grid;
Figure_appearance

subplot(1,5,4) % ---------- plot of yy_DL_TC_NL
plot(x,yy_DL_TC_NL,'color',colcol(2,:),'linewidth',4);
hold on
[HWg(4)]=find(yy_DL_TC_NL>=max(yy_DL_TC_NL/2),1);
plot(x(HWg(4)),yy_DL_TC_NL(HWg(4)),'o','color',colcol(2,:),'linewidth',2,'MarkerSize',15,'MarkerFaceColor','w');
line([ (del(1,2)) (del(1,2))],[0 0.005],'color',colcol(1,:),'linewidth',2)

xlim([-0.005 0.10])
ylim([-0.001 0.025])
title('+ NL + D + TC')
grid;
Figure_appearance


subplot(1,5,5) % ---------- plot of phase estimations across all 4 conditions
hold on
plot(1000*del(1,2)*0,1000*x(HWg(1)),'o','color',colcol(1,:),'linewidth',2,'MarkerSize',10,'MarkerFaceColor','w');

for z=2:4 % times steps
    plot(1000*del(1,2),1000*x(HWg(z)),'o','color',colcol(2,:),'linewidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end

grid; Figure_appearance
xlim([0 3])
ylim([0 20])

footer(['NL =' num2str(NL) ' | DEL12 =' num2str(del(1,2)) ' | TC1 =' num2str(a1(1))])

%% FIG 2B
% figure illustrating the case of horizontal propagation within one area

delH=0.01; % delay from each node sample along the horizontal network
tab=[1 1; 2 2;  3 3;  4 4; 5 5 ]; % receiving layer / time step
ER(1,1,:)=CortConv( g1  ,squeeze(TC(1,:)),0,x,NL);      % CortConv is a function for the cortical integration 

clear HPSI % horizontal post-synaptic integration 
for z=1:size(tab,1) % times steps
    j=tab(z,1); k=tab(z,2);
    if j>1, 
        HPSI(j,k,:)=CortConv( squeeze(ER(1,1,:))  ,squeeze(TC(j,:)), delH*(j-1),x,NL);    % Horizontal PostSynaptic Integration (PSI) of layer j (Time constant j) of FF input (coming from layer j-1)  , step k
        ER(j,k,:)=squeeze(HPSI(j,k,:));
    end
end


figure(22)
clf, set(gcf,'PaperPosition',[0.25 0.25 50 10]); 


for z=2:size(tab,1) % t---------- Times steps
    j=tab(z,1);

    subplot(1,size(tab,1),z-1) % plot the hPSI at each step
    clcl=1;

    plot(x,squeeze(ER(j,j,:)),'color',colcol(clcl,:),'linewidth',4);
    hold on

    [HWz(j,j)]=find(ER(j,j,:)>=max(squeeze(ER(j,j,:)))/2,1); % extract and plot an approximation of the phase (half height) HWz
    plot(x(HWz(j,j)),squeeze(ER(j,j,HWz(j,j))),'o','color',colcol(clcl,:),'linewidth',2,'MarkerSize',15,'MarkerFaceColor','w');

    line([ delH*(j-1) delH*(j-1)],[0 0.005],'color',colcol(clcl,:),'linewidth',2)

    xlim([-0.005 0.25])
    ylim([-0.001 0.025])
    grid;
    Figure_appearance

end

[HWz(1,1)]=find(ER(1,1,:)>=max(squeeze(ER(1,1,:)))/2,1); % extract an approximation of the phase (half height) HWz


subplot(1,5,5) % ---------- final plot of the evolution of the phase (HWz) across all steps
hold on
line([0 40],1000*x(HWz(1,1))+[0 40],'color','k','linewidth',3,'linestyle',':')
line([0 40],1000*x(HWz(1,1))+[0 80],'color','k','linewidth',2,'linestyle',':')
line([0 40],1000*x(HWz(1,1))+[0 120],'color','k','linewidth',1,'linestyle',':')
plot(0,1000*x(HWz(1,1)),'o','color',colcol(1,:),'linewidth',2,'MarkerSize',10,'MarkerFaceColor','w');
for z=2:5 % times steps
    j=z+(k-1)*2;
    plot(1000*(z-1).*delH,1000*x(HWz(z,z)),'o','color',colcol(1,:),'linewidth',2,'MarkerSize',10,'MarkerFaceColor','w');
end

ylim([0 90])
grid; Figure_appearance

footer(['NL =' num2str(NL) ' | GainFB =' num2str(GainFB) ' | DEL =' num2str(delH) ' | TC1 =' num2str(a1(1))])



%% INTER-CORTICAL SIMULATION OF EVOKED RESPONSE (ER) AS SUM OF POSTSYNAPTIC INTEGRATION (PSI) FROM Feedforward (FF) & Feedback (FB)

clear ER % ER(i,j,:) is Evoked Response in time for hierarchical level i, and cycle step k

ER(1,1,:)=CortConv( g1  ,squeeze(TC(1,:)),0,x,NL);      % CortConv is a function for the cortical integration (see below)
                                                        % the first step is thus to integrate an alpha function (g1) 
ER(2,1,:)=x.*0;
ER(3,1,:)=x.*0;
ER(4,2,:)=x.*0;
ER(5,3,:)=x.*0;
ER(6,4,:)=x.*0;

clear tab PSI
% tab is attributing, for each time step, the corresponding hierarchical level (tab(z,1)) and the cycle step (tab(z,2))
tab=[2 2;  1 3; 3 3; 2 4; 4 4;  1 5 ;3 5  ;5 5;  2 6 ; 4 6 ; 1 7 ; 3 7 ; 5 7; 2 8 ; 4 8 ; 1 9 ; 3 9 ; 5 9 ; 4 10 ]; % receiving layer / time step


for z=1:size(tab,1) % times steps
    j=tab(z,1); k=tab(z,2);
    if j>1, 
        PSI(j-1,j,k,:)=CortConv( squeeze(ER(j-1,k-1,:))  ,squeeze(TC(j,:)), del(j-1,j),x,NL);    % PostSynaptic Integration (PSI) of layer j (Time constant j) of FF input (coming from layer j-1)  , step k
        PSI(j+1,j,k,:)=CortConv( squeeze(ER(j+1,k-1,:))  ,squeeze(TC(j,:)), del(j+1,j),x,NL);    % PostSynaptic Integration (PSI) of layer j (Time constant j) of FB input (coming from layer j+1)  , step k
        ER(j,k,:)=squeeze(PSI(j-1,j,k,:))+GainFB.*squeeze(PSI(j+1,j,k,:));
    else
        PSI(j+1,j,k,:)=CortConv( squeeze(ER(j+1,k-1,:))  ,squeeze(TC(j,:)), del(j+1,j),x,NL);    % PostSynaptic Integration (PSI) of layer j (Time constant j) of FB input (coming from layer j+1)  , step k
        ER(j,k,:)=squeeze(ER(j,k-2,:))+GainFB.*squeeze(PSI(j+1,j,k,:));
    end
end




%% SUPPL FIG
% Not in the paper, but shows all the cycles (columns)
% iterations across all hierarchical levels (color-coded, row), 
% to calculate HWx, estimation of phase in each steps
% feedforward is dotted, feedback dashed

figure(24)
clf, set(gcf,'PaperPosition',[0.25 0.25 40 25]); 

tab=[1 1; 2 2;  1 3; 3 3; 2 4; 4 4;  1 5 ;3 5  ;  2 6 ; 4 6 ; 1 7 ; 3 7 ; 2 8 ; 4 8 ; 3 9 ; 4 10   ]; % receiving layer / time step


for z=1:size(tab,1) % ---------- times steps
    j=tab(z,1); k=tab(z,2);

    subplot('Position',[k/11-0.02 j/5-0.1 0.1 0.2])

    plot(x,squeeze(ER(j,k,:)),'color',colcol(j,:),'linewidth',4);
    hold on
    if j>1
        plot(x,squeeze(PSI(j-1,j,k,:)),'color',colcol(j-1,:),'linestyle',':','linewidth',2);
        plot(x,GainFB.*squeeze(PSI(j+1,j,k,:)),'color',colcol(j+1,:),'linestyle','--','linewidth',2);
    else
        if k>2, plot(x,squeeze(ER(j,k-2,:)),'k:','linewidth',2); end
        plot(x,GainFB.*squeeze(PSI(j+1,j,k,:)),'k--','linewidth',2);
    end

    % GENERATING FF (dotted) and FB (dashed) arrows
    gcap=get(gca,'Position');
    if j<4
        arh=annotation('arrow', [gcap(1)+gcap(3)/2 gcap(1)+gcap(3)/2+0.03 ],[gcap(2)+gcap(4)/1.5 gcap(2)+gcap(4)/1.5+0.05] );
        arh.LineWidth = 2; arh.LineStyle = ':';

    end
    if j>1
        arh=annotation('arrow', [gcap(1)+gcap(3)/4 gcap(1)+gcap(3)/4+0.03 ],[gcap(2)-0.02 gcap(2)-0.05] );
        arh.LineWidth = 2; arh.LineStyle = '--';
    end

    [HWx(j,k)]=find(ER(j,k,:)>=max(squeeze(ER(j,k,:)))/2,1); % extract phase as half-width 
    plot(x(HWx(j,k)),squeeze(ER(j,k,HWx(j,k))),'o','color',colcol(j,:),'linewidth',2,'MarkerSize',15,'MarkerFaceColor','w');

    line([ (del(1,2))*(j-1) (del(1,2))*(j-1)],[0 0.005],'color',colcol(j,:),'linewidth',2)

    xlim([-0.005 0.25])
    ylim([-0.001 0.025])
    grid;
    Figure_appearance

end
footer(['NL =' num2str(NL) ' | GainFB =' num2str(GainFB) ' | DEL12 =' num2str(del(1,2)) ' | DEL21 =' num2str(del(2,1)) ' | TC1 =' num2str(a1(1))])


%% FIG 2C
% Figure compiling inter-cortical interactions across hierarchy
% (color-coded) and cycle steps (rows)

figure(23)
clf,set(gcf,'PaperPosition',[0.25 0.25 50 10]); 

mxk=4 ;
tab=[1 1; 2 2; 3 3; 4 4];

for k=1:mxk
    subplot(1,mxk+1,k)
    hold on
    for z=1:4 % ---------- times steps
        j=z+(k-1)*2;
        plot(x, squeeze(ER(z,j,:)) ,'color',colcol(z,:),'linewidth',3)
        plot(x(HWx(z,j)),squeeze(ER(z,j,HWx(z,j))),'o','color',colcol(z,:),'linewidth',2,'MarkerSize',10,'MarkerFaceColor','w');
        line([ (del(1,2))*(z-1) (del(1,2))*(z-1)],[0 0.005],'color',colcol(z,:),'linewidth',2)
       
        xlim([0 0.2])
        ylim([0 0.025])
        grid; Figure_appearance
    end
end

subplot(1,mxk+1,mxk+1) % ---------- last plot recapitulating all phase (HWx)

hold on
line([0 7],1000*x(HWx(1,1))+[0 7],'color','k','linewidth',3,'linestyle',':')
line([0 7],1000*x(HWx(1,1))+[0 14],'color','k','linewidth',2,'linestyle',':')
line([0 7],1000*x(HWx(1,1))+[0 28],'color','k','linewidth',1,'linestyle',':')
line([0 7],1000*x(HWx(1,1))+[0 56],'color','k','linewidth',0.5,'linestyle',':')
line([0 7],1000*x(HWx(1,1))+[0 112],'color','k','linewidth',0.25,'linestyle',':')

sz=size(HWx);
for k=1:mxk

    zz=1:4;
    jj=zz+(k-1)*2;
    plot(1000*(zz-1)*del(1,2),1000*x(HWx(sub2ind(sz,zz,jj))),'-','color','k','linewidth',k/2);
    
    for z=1:4 % times steps
        j=z+(k-1)*2;
        plot(1000*(z-1)*del(1,2),1000*x(HWx(z,j)),'o','color',colcol(z,:),'linewidth',2,'MarkerSize',10,'MarkerFaceColor','w');
    end

end
xlim([0 6.5])
ylim([0 90])

grid; Figure_appearance

footer(['NL =' num2str(NL) ' | GainFB =' num2str(GainFB) ' | DEL12 =' num2str(del(1,2)) ' | DEL21 =' num2str(del(2,1)) ' | TC1 =' num2str(a1(1))])




%% Convolution of the output of a structure (G) with a recipient structure with time constant (g), after a propagation delay d - time is x

function Gt = CortConv(G,g,dd,x,NL)
% G : input integrated by a recipient structure with time constant (g)
% a power law exponent NL, and after a propagation delay d - time is x

if sum(G)~=0,
    Gs=G.^NL./sum(G.^NL); % Nonlinearity to go from Voltage to Spike (Chen et al 2012)
else
    Gs=G;
end

if sum(g)~=0,
    Gsg=conv(Gs,g,'full'); % convolution between the output of structure 1 and the postsynaptic integration time constant g
else
    Gsg=Gs;
end

if sum(Gsg)~=0,
    Gsg=Gsg./sum(Gsg); %  normalizing
end

d=find(x>=dd,1);
Gt=[zeros(1,d) Gsg(1:length(x)-d)]; % adding delay d


end



%% Figure appearance

function [] = Figure_appearance()
box off
set(gca,'TickDir','out')
set(gca,'XScale', 'linear')
set(gca,'fontangle','oblique')
set(gca,'fontsize',16)
set(gca,'ticklength',[0.01 0.025])
axis square
set(gcf,'Color','w')
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
set(gcf,'position',pos);
end
