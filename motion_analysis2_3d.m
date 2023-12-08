clc
clear
close all

%% load static data
msgbox('Please select and load static data!');
[file path]=uigetfile('*.xlsx');
data=xlsread([path file]);

% pelvis markers, meters
LASIS_st=0.001*mean(data(:,3:5));
RASIS_st=0.001*mean(data(:,6:8));
LPSIS_st=0.001*mean(data(:,9:11));
RPSIS_st=0.001*mean(data(:,12:14));

% knee markers
RLKNE_st=0.001*mean(data(:,36:38));
RMKNE_st=0.001*mean(data(:,111:113));

% ankle markers
RLANK_st=0.001*mean(data(:,42:44));
RMANK_st=0.001*mean(data(:,114:116));

% foot markers
RMET1_st=0.001*mean(data(:,126:128));
RMET5_st=0.001*mean(data(:,78:80));
RTO_st=0.001*mean(data(:,48:50));
RHE_st=0.001*mean(data(:,45:47));

Opelvis_st=0.5*(RASIS_st+LASIS_st);
i_pel_st=(RASIS_st-Opelvis_st)/norm(RASIS_st-Opelvis_st);
v_pel_st=(Opelvis_st-0.5*(RPSIS_st+LPSIS_st))/...
    norm(Opelvis_st-0.5*(RPSIS_st+LPSIS_st));
k_pel_st=cross(i_pel_st,v_pel_st);
j_pel_st=cross(k_pel_st,i_pel_st);  
R_pel_st=[i_pel_st;j_pel_st;k_pel_st];

Ohipr_st=[0.36 -0.19 -0.3]*norm(RASIS_st-LASIS_st);
Ohipl_st=[-0.36 -0.19 -0.3]*norm(RASIS_st-LASIS_st);
OhipR_st=R_pel_st'*Ohipr_st'+Opelvis_st';
OhipL_st=R_pel_st'*Ohipl_st'+Opelvis_st';

OhipR_st=OhipR_st';
OkneeR_st=0.5*(RMKNE_st+RLKNE_st);
OankleR_st=0.5*(RLANK_st+RMANK_st);

% segment anatomical length
L_pelvis=norm(Opelvis_st-OhipR_st);
L_thigh=norm(OhipR_st-OkneeR_st);
L_shank=norm(OkneeR_st-OankleR_st);
L_foot=norm(OankleR_st-0.5*(RMET1_st+RMET5_st));
L_heel=norm(RHE_st-OankleR_st);

% femur anatomical LCS
k_fem=(OhipR_st-0.5*(RMKNE_st+RLKNE_st))/norm(OhipR_st-0.5*(RMKNE_st+RLKNE_st));
v_fem=(RLKNE_st-RMKNE_st)/norm(RLKNE_st-RMKNE_st);
j_fem=cross(k_fem,v_fem);
i_fem=cross(j_fem,k_fem);
R_fem=[i_fem;j_fem;k_fem];
    
% tibial anatomical LCS
k_tib=(OkneeR_st-0.5*(RMANK_st+RLANK_st))/norm(OkneeR_st-0.5*(RMANK_st+RLANK_st));
V_tib=(RLANK_st-RMANK_st)/norm(RLANK_st-RMANK_st);
j_tib=cross(k_tib,V_tib);
i_tib=cross(j_tib,k_tib);
R_tib=[i_tib;j_tib;k_tib];  
    
%% load motion data
msgbox('Please select and load motion data!');
[file path]=uigetfile('*.xlsx');
data=xlsread([path file]);
frame=data(:,1)-data(1,1);
FR=100; % frame rate, Hz
t=frame/FR; %time, sec
n=length(t);
dt=1/FR;

% pelvis markers, meters
LASIS=0.001*data(:,3:5);
RASIS=0.001*data(:,6:8);
LPSIS=0.001*data(:,9:11);
RPSIS=0.001*data(:,12:14);

% knee markers
% RLKNE=0.001*data(:,36:38);
% RMKNE=0.001*data(:,111:113);

% ankle markers
RLANK=0.001*data(:,42:44);
RMANK=0.001*data(:,114:116);

% foot markers
RMET1=0.001*data(:,126:128);
RMET5=0.001*data(:,78:80);
RTO=0.001*data(:,48:50);
RHE=0.001*data(:,45:47);

% thigh markers
RTHI=0.001*data(:,33:35);
RTHIBU=0.001*data(:,51:53);
RTHIFU=0.001*data(:,54:56);
OthighR=(RTHI+RTHIBU+RTHIFU)/3;

% shank markers
RTIB=0.001*data(:,39:41);
RTIBBU=0.001*data(:,63:65);
RTIBFU=0.001*data(:,66:68);
OshankR=(RTIB+RTIBBU+RTIBFU)/3;

%% load force data
fdata=xlsread([path file],2);
% index=10:10:10*n;
% Fdata=Fdata(index,:);

for i=10:10:length(fdata) %10*n
    Fdata(0.1*i,:)=mean(fdata(i-9:i,:));
end

GRF=-Fdata(:,3:5); Fx=GRF(:,1); Fy=GRF(:,2); Fz=GRF(:,3);
Mz=0.001*Fdata(:,8); 
COP=0.001*Fdata(:,9:11); Cx=COP(:,1); Cy=COP(:,2); Cz=COP(:,3);
g=[0 0 -9.81];

%% segment local coordinate systems (LCS)

% Pelvis LCS
Opelvis=0.5*(RASIS+LASIS);
for i=1:n
    i_pel(i,:)=(RASIS(i,:)-Opelvis(i,:))/norm(RASIS(i,:)-Opelvis(i,:));
    v_pel(i,:)=(Opelvis(i,:)-0.5*(RPSIS(i,:)+LPSIS(i,:)))/...
        norm(Opelvis(i,:)-0.5*(RPSIS(i,:)+LPSIS(i,:)));
    k_pel(i,:)=cross(i_pel(i,:),v_pel(i,:));
    j_pel(i,:)=cross(k_pel(i,:),i_pel(i,:));  
    R_pel(:,:,i)=[i_pel(i,:);j_pel(i,:);k_pel(i,:)];
end

% thigh LCS
for i=1:n
    j_th(i,:)=(RTHIFU(i,:)-RTHIBU(i,:))/norm(RTHIFU(i,:)-RTHIBU(i,:));
    v_th(i,:)=(0.5*(RTHIFU(i,:)+RTHIBU(i,:))-RTHI(i,:))/norm(0.5*(RTHIFU(i,:)+RTHIBU(i,:))-RTHI(i,:));
    i_th(i,:)=cross(j_th(i,:),v_th(i,:));
    k_th(i,:)=cross(i_th(i,:),j_th(i,:));
    R_th(:,:,i)=[i_th(i,:);j_th(i,:);k_th(i,:)];
end

% shank LCS
for i=1:n
    j_sh(i,:)=(RTIBFU(i,:)-RTIBBU(i,:))/norm(RTIBFU(i,:)-RTIBBU(i,:));
    v_sh(i,:)=(0.5*(RTIBFU(i,:)+RTIBBU(i,:))-RTIB(i,:))/norm(0.5*(RTIBFU(i,:)+RTIBBU(i,:))-RTIB(i,:));
    i_sh(i,:)=cross(j_sh(i,:),v_sh(i,:));
    k_sh(i,:)=cross(i_sh(i,:),j_sh(i,:));
    R_sh(:,:,i)=[i_sh(i,:);j_sh(i,:);k_sh(i,:)];    
end

% hip center
for i=1:n
    Ohipr(i,:)=[0.36 -0.19 -0.3]*norm(RASIS(i,:)-LASIS(i,:));
%     Ohipl(i,:)=[-0.36 -0.19 -0.3]*norm(RASIS(i,:)-LASIS(i,:));
    OhipR(i,:)=R_pel(:,:,i)'*Ohipr(i,:)'+Opelvis(i,:)';
%     OhipL(i,:)=R_pel(:,:,i)'*Ohipl(i,:)'+Opelvis(i,:)';
    v_hp(i,:)=(Opelvis(i,:)-OhipR(i,:))/norm(Opelvis(i,:)-OhipR(i,:)); %unit vector from hip to pelvis origin
end

%     k_th(i,:)=(OhipR(i,:)-0.5*(RMKNE(i,:)+RLKNE(i,:)))/norm(OhipR(i,:)-0.5*(RMKNE(i,:)+RLKNE(i,:)));
%     v_th(i,:)=(RLKNE(i,:)-RMKNE(i,:))/norm(RLKNE(i,:)-RMKNE(i,:));
%     j_th(i,:)=cross(k_th(i,:),v_th(i,:));
%     i_th(i,:)=cross(j_th(i,:),k_th(i,:));
%     R_th(:,:,i)=[i_th(i,:);j_th(i,:);k_th(i,:)];    
% end

% shank LCS
% OkneeR=0.5*(RMKNE+RLKNE);
% for i=1:n
%     k_sh(i,:)=(OkneeR(i,:)-0.5*(RMANK(i,:)+RLANK(i,:)))/norm(OkneeR(i,:)-0.5*(RMANK(i,:)+RLANK(i,:)));
%     v_sh(i,:)=(RLANK(i,:)-RMANK(i,:))/norm(RLANK(i,:)-RMANK(i,:));
%     j_sh(i,:)=cross(k_sh(i,:),v_sh(i,:));
%     i_sh(i,:)=cross(j_sh(i,:),k_sh(i,:));
%     R_sh(:,:,i)=[i_sh(i,:);j_sh(i,:);k_sh(i,:)];    
% end 

% foot LCS
OankleR=0.5*(RLANK+RMANK);
for i=1:n
    k_fo(i,:)=(OankleR(i,:)-0.5*(RMET1(i,:)+RMET5(i,:)))/norm(OankleR(i,:)-0.5*(RMET1(i,:)+RMET5(i,:)));
    v_fo(i,:)=(RLANK(i,:)-RMANK(i,:))/norm(RLANK(i,:)-RMANK(i,:));
    j_fo(i,:)=cross(k_fo(i,:),v_fo(i,:));
    i_fo(i,:)=cross(j_fo(i,:),k_fo(i,:));
    R_fo(:,:,i)=[i_fo(i,:);j_fo(i,:);k_fo(i,:)];  
    v_ah(i,:)=(RHE(i,:)-OankleR(i,:))/norm(RHE(i,:)-OankleR(i,:));
end 
    
%% correcting anatomical joint positions
% OankleR=0.5*(RMET1+RMET5)+L_foot*k_fo;
% OkneeR=OankleR+L_shank*k_sh;
% OhipR=OkneeR+L_thigh*k_th;
for i=1:n
    OankleR(i,:)=R_fo(:,:,i)'*[0;0;L_foot]+0.5*(RMET1(i,:)'+RMET5(i,:)');
    OkneeR(i,:)=R_sh(:,:,i)'*[0;0;L_shank]+OankleR(i,:)';
    OhipR(i,:)=R_th(:,:,i)'*[0;0;L_thigh]+OkneeR(i,:)';
end
Opel=OhipR+L_pelvis*v_hp;
rpn=Opel-Opelvis; %different between old and new pelvis origin 
RASIS=RASIS+rpn; RPSIS=RPSIS+rpn;
LASIS=LASIS+rpn; LPSIS=LPSIS+rpn;
RHE=OankleR+L_heel*v_ah;

%% COM calculation
mass=102.5; %body mass, Kg
m_th=0.1*mass;
m_sh=0.0465*mass;
m_fo=0.0165*mass;

Rp_th=norm(RASIS_st-OhipR_st); % thigh proximal radius
Rd_th=0.5*norm(RLKNE_st-RMKNE_st); % thigh distal radius
Rp_sh=Rd_th;
Rd_sh=0.5*norm(RLANK_st-RMANK_st);
Rp_fo=Rd_sh;
Rd_fo=0.5*norm(RMET1_st-RMET5_st);

x_th=Rd_th/Rp_th; x_sh=Rd_sh/Rp_sh; x_fo=Rd_fo/Rp_fo; 

C=@(x) (1+2*x+3*x^2)/(4*(1+x+x^2));
c=@(x) C(x)*(x<=1)+(1-C(x))*(x>1);

r_th=[0;0;-c(x_th)*L_thigh]; %com in LCS
r_sh=[0;0;-c(x_sh)*L_shank];
r_fo=[0;0;-c(x_fo)*L_foot];

for i=1:n
    rcm_th(i,:)=R_th(:,:,i)'*r_th+OhipR(i,:)';
    rcm_sh(i,:)=R_sh(:,:,i)'*r_sh+OkneeR(i,:)';
    rcm_fo(i,:)=R_fo(:,:,i)'*r_fo+OankleR(i,:)';
end

% COM velocities and accelerations
[vcm_th acm_th]=num_diff(t,rcm_th);
[vcm_sh acm_sh]=num_diff(t,rcm_sh);
[vcm_fo acm_fo]=num_diff(t,rcm_fo);

%% moment of inertia calculation
I_th=MOI(m_th,L_thigh,Rp_th,Rd_th);
I_sh=MOI(m_sh,L_shank,Rp_sh,Rd_sh);
I_fo=MOI(m_fo,L_foot,Rp_fo,Rd_fo);

%% muscle attachments
% tibialis anterior
tib_ins=OankleR-0.3*L_foot*k_fo;
tib_org=OkneeR-0.35*L_shank*k_sh+0.02*L_shank*j_sh;

% soleus
sol_ins=RHE;
sol_org=OkneeR-0.25*L_shank*k_sh-0.02*L_shank*j_sh;

% gastrocnemeus
gas_ins=RHE;
gas_org1=OhipR-0.84*L_thigh*k_th-0.02*L_thigh*j_th-0.07*L_thigh*i_th;
gas_org2=OhipR-0.84*L_thigh*k_th-0.02*L_thigh*j_th+0.07*L_thigh*i_th;

% quadriceps
vi_ins=OkneeR+0.14*L_thigh*j_th;
vi_org=OhipR-0.3*L_thigh*k_th+0.03*L_thigh*j_th;

vm_ins=OkneeR+0.14*L_thigh*j_th-0.01*L_thigh*i_th;
vm_org=OhipR-0.4*L_thigh*k_th+0.03*L_thigh*j_th-0.05*L_thigh*i_th;
vmo_org=OhipR-0.8*L_thigh*k_th+0.03*L_thigh*j_th-0.08*L_thigh*i_th;

vl_ins=OkneeR+0.14*L_thigh*j_th+0.01*L_thigh*i_th;
vl_org=OhipR-0.4*L_thigh*k_th+0.03*L_thigh*j_th+0.05*L_thigh*i_th;
vlo_org=OhipR-0.65*L_thigh*k_th+0.03*L_thigh*j_th+0.07*L_thigh*i_th;

rf_ins=OkneeR+0.16*L_thigh*j_th;
rf_org=0.5*(OhipR+RASIS);

% hamstrings
sm_ins=OkneeR-0.25*L_shank*k_sh-0.07*L_shank*i_sh;
sm_org=OhipR-0.4*L_pelvis*k_pel-0.4*L_pelvis*j_pel;

st_ins=OkneeR-0.25*L_shank*k_sh-0.08*L_shank*i_sh;
st_org=OhipR-0.4*L_pelvis*k_pel-0.4*L_pelvis*j_pel-0.01*L_pelvis*i_pel;

bf_ins=OkneeR-0.25*L_shank*k_sh+0.07*L_shank*i_sh;
bf_org1=OhipR-0.4*L_pelvis*k_pel-0.4*L_pelvis*j_pel+0.005*L_pelvis*i_pel;
bf_org2=OhipR-0.5*L_thigh*k_th-0.03*L_thigh*j_th;

%% showing the model 
figure('WindowState','Maximized')
plotx=[RASIS(:,1) LASIS(:,1) LPSIS(:,1) RPSIS(:,1) RASIS(:,1) OhipR(:,1) OkneeR(:,1) OankleR(:,1) RMET1(:,1) RMET5(:,1) OankleR(:,1) RHE(:,1)];
ploty=[RASIS(:,2) LASIS(:,2) LPSIS(:,2) RPSIS(:,2) RASIS(:,2) OhipR(:,2) OkneeR(:,2) OankleR(:,2) RMET1(:,2) RMET5(:,2) OankleR(:,2) RHE(:,2)];
plotz=[RASIS(:,3) LASIS(:,3) LPSIS(:,3) RPSIS(:,3) RASIS(:,3) OhipR(:,3) OkneeR(:,3) OankleR(:,3) RMET1(:,3) RMET5(:,3) OankleR(:,3) RHE(:,3)];

plotcmx=[rcm_th(:,1) rcm_sh(:,1) rcm_fo(:,1)];
plotcmy=[rcm_th(:,2) rcm_sh(:,2) rcm_fo(:,2)];
plotcmz=[rcm_th(:,3) rcm_sh(:,3) rcm_fo(:,3)];

for i=1:5:n
    plot3(plotx(i,:),ploty(i,:),plotz(i,:),'*-');
    hold on
    plot3(plotcmx(i,:),plotcmy(i,:),plotcmz(i,:),'blacko');
    plot3([tib_org(i,1) tib_ins(i,1)],[tib_org(i,2) tib_ins(i,2)],[tib_org(i,3) tib_ins(i,3)],'r','LineWidth',2.5);
    plot3([sol_org(i,1) sol_ins(i,1)],[sol_org(i,2) sol_ins(i,2)],[sol_org(i,3) sol_ins(i,3)],'r','LineWidth',2.5);
    plot3([gas_org1(i,1) gas_ins(i,1) gas_org2(i,1)],[gas_org1(i,2) gas_ins(i,2) gas_org2(i,2)],[gas_org1(i,3) gas_ins(i,3) gas_org2(i,3)],'r','LineWidth',2.5);
    plot3([vi_org(i,1) vi_ins(i,1)],[vi_org(i,2) vi_ins(i,2)],[vi_org(i,3) vi_ins(i,3)],'r','LineWidth',2.5);
    plot3([vm_org(i,1) vm_ins(i,1) vmo_org(i,1)],[vm_org(i,2) vm_ins(i,2) vmo_org(i,2)],[vm_org(i,3) vm_ins(i,3) vmo_org(i,3)],'r','LineWidth',2.5);
    plot3([vl_org(i,1) vl_ins(i,1) vlo_org(i,1)],[vl_org(i,2) vl_ins(i,2) vlo_org(i,2)],[vl_org(i,3) vl_ins(i,3) vlo_org(i,3)],'r','LineWidth',2.5);
    plot3([rf_org(i,1) rf_ins(i,1)],[rf_org(i,2) rf_ins(i,2)],[rf_org(i,3) rf_ins(i,3)],'r','LineWidth',2.5);
    plot3([sm_org(i,1) sm_ins(i,1)],[sm_org(i,2) sm_ins(i,2)],[sm_org(i,3) sm_ins(i,3)],'r','LineWidth',2.5);
    plot3([st_org(i,1) st_ins(i,1)],[st_org(i,2) st_ins(i,2)],[st_org(i,3) st_ins(i,3)],'r','LineWidth',2.5);
    plot3([bf_org1(i,1) bf_ins(i,1) bf_org2(i,1)],[bf_org1(i,2) bf_ins(i,2) bf_org2(i,2)],[bf_org1(i,3) bf_ins(i,3) bf_org2(i,3)],'r','LineWidth',2.5);
    
%     plot3([OhipR(i,1) OhipL(i,1)],[OhipR(i,2) OhipL(i,2)],[OhipR(i,3) OhipL(i,3)],'bo');
    plot_axis(Opel(i,:),i_pel(i,:),j_pel(i,:),k_pel(i,:),RASIS(i,:),LASIS(i,:),0.5*(RPSIS(i,:)+LPSIS(i,:)));
%     plot_axis(OhipR(i,:),i_th(i,:),j_th(i,:),k_th(i,:),OhipR(i,:),RLKNE(i,:),RMKNE(i,:));
    plot_axis(OthighR(i,:),i_th(i,:),j_th(i,:),k_th(i,:),RTHI(i,:),RTHIBU(i,:),RTHIFU(i,:));
%     plot_axis(OkneeR(i,:),i_sh(i,:),j_sh(i,:),k_sh(i,:),OkneeR(i,:),RLANK(i,:),RMANK(i,:));
    plot_axis(OshankR(i,:),i_sh(i,:),j_sh(i,:),k_sh(i,:),RTIB(i,:),RTIBBU(i,:),RTIBFU(i,:));
    plot_axis(OankleR(i,:),i_fo(i,:),j_fo(i,:),k_fo(i,:),OankleR(i,:),RMET1(i,:),RMET5(i,:));

    quiver3(COP(i,1),COP(i,2),COP(i,3),0.001*Fx(i),0.001*Fy(i),0.001*Fz(i),'m','LineWidth',2);
    view([1 0.5 1]); title(sprintf('motion time=%.2f sec',t(i)));
    axis equal; 
    xlim([min(plotx(:))-0.1 max(plotx(:))+0.1]);
    ylim([min(ploty(:))-0.1 max(ploty(:))+0.1]);
    zlim([min(plotz(:))-0.1 max(plotz(:))+0.1]);    
    pause(dt);
    hold off
end

%% inverse kinematics

% segment angles
for i=1:n
   [a_pel(i,:) b_pel(i,:) g_pel(i,:)]=angles(R_pel(:,:,i));
   [a_th(i,:) b_th(i,:) g_th(i,:)]=angles(R_th(:,:,i));
   [a_sh(i,:) b_sh(i,:) g_sh(i,:)]=angles(R_sh(:,:,i));
   [a_fo(i,:) b_fo(i,:) g_fo(i,:)]=angles(R_fo(:,:,i));
end

% joint angles
for i=1:n
   [a_h(i,:) b_h(i,:) g_h(i,:)]=angles(R_pel(:,:,i)'*R_th(:,:,i));
   [a_k(i,:) b_k(i,:) g_k(i,:)]=angles(R_th(:,:,i)'*R_sh(:,:,i));
   [a_a(i,:) b_a(i,:) g_a(i,:)]=angles(R_sh(:,:,i)'*R_fo(:,:,i));
end

% velocities and accelerations
[omga_pel alfa_pel]=num_diff(t,pi/180*[a_pel b_pel g_pel]);
[omga_th alfa_th]=num_diff(t,pi/180*[a_th b_th g_th]);
[omga_sh alfa_sh]=num_diff(t,pi/180*[a_sh b_sh g_sh]);
[omga_fo alfa_fo]=num_diff(t,pi/180*[a_fo b_fo g_fo]);

[omga_h alfa_h]=num_diff(t,pi/180*[a_h b_h g_h]);
[omga_k alfa_k]=num_diff(t,pi/180*[a_k b_k g_k]);
[omga_a alfa_a]=num_diff(t,pi/180*[a_a b_a g_a]);

%% kinetic calculations
% net joint forces
Fa=m_fo*(acm_fo-g)-GRF;
Fk=m_sh*(acm_sh-g)+Fa;
Fh=m_th*(acm_th-g)+Fk;

% net joint moments
for i=1:n
    %ankle joint
    omgfo(i,:)=R_fo(:,:,i)*omga_fo(i,:)';
    alffo(i,:)=R_fo(:,:,i)*alfa_fo(i,:)';
    Tin_fo(i,:)=I_fo*alffo(i,:)'+cross(omgfo(i,:)',I_fo*omgfo(i,:)');
    Tin_fo(i,:)=R_fo(:,:,i)'*Tin_fo(i,:)';
    Tz(i,:)=Mz(i)-COP(i,1)*Fy(i)-COP(i,2)*Fx(i);
    Tgrf(i,:)=[0 0 Tz(i)]; 
    ra_f(i,:)=rcm_fo(i,:)-OankleR(i,:);
    ra_grf(i,:)=COP(i,:)-OankleR(i,:);
    Ta(i,:)=Tin_fo(i,:)'-Tgrf(i,:)'-cross((ra_grf(i,:)'-ra_f(i,:)'),GRF(i,:)')...
        +cross(ra_f(i,:)',Fa(i,:)');
    %knee joint
    omgsh(i,:)=R_sh(:,:,i)*omga_sh(i,:)';
    alfsh(i,:)=R_sh(:,:,i)*alfa_sh(i,:)';
    Tin_sh(i,:)=I_sh*alfsh(i,:)'+cross(omgsh(i,:)',I_sh*omgsh(i,:)');
    Tin_sh(i,:)=R_sh(:,:,i)'*Tin_sh(i,:)';
    rk_l(i,:)=rcm_sh(i,:)-OkneeR(i,:);
    rk_a(i,:)=OankleR(i,:)-OkneeR(i,:);
    Tk(i,:)=Tin_sh(i,:)'+Ta(i,:)'+cross(rk_l(i,:)',Fk(i,:)')...
        +cross((rk_a(i,:)'-rk_l(i,:)'),Fa(i,:)');
    %hip joint
    omgth(i,:)=R_th(:,:,i)*omga_th(i,:)';
    alfth(i,:)=R_th(:,:,i)*alfa_th(i,:)';
    Tin_th(i,:)=I_th*alfth(i,:)'+cross(omgth(i,:)',I_th*omgsh(i,:)');
    Tin_th(i,:)=R_th(:,:,i)'*Tin_th(i,:)';
    rh_t(i,:)=rcm_th(i,:)-OhipR(i,:);
    rh_k(i,:)=OkneeR(i,:)-OhipR(i,:);
    Th(i,:)=Tin_th(i,:)'+Tk(i,:)'+cross(rh_t(i,:)',Fh(i,:)')...
        +cross((rh_k(i,:)'-rh_t(i,:)'),Fk(i,:)');    
end

% joint powers
Pa=Ta.*omga_a;
Pk=Tk.*omga_k;
Ph=Th.*omga_h;

%% muscle & joint force calculation

% tibialis anterior
v_tib=tib_org-tib_ins;
L_tib=sqrt(v_tib(:,1).^2+v_tib(:,2).^2+v_tib(:,3).^2);
e_tib=v_tib./L_tib;
r_tib=tib_ins-OankleR;
for i=1:n
    d_tib(i,:)=cross(r_tib(i,:),v_tib(i,:))/L_tib(i);
end

% soleus
v_sol=sol_org-sol_ins;
L_sol=sqrt(v_sol(:,1).^2+v_sol(:,2).^2+v_sol(:,3).^2);
e_sol=v_sol./L_sol;
r_sol=sol_ins-OankleR;
for i=1:n
    d_sol(i,:)=cross(r_sol(i,:),v_sol(i,:))/L_sol(i);
end

% gastrocnemius
v_gas1=gas_org1-gas_ins; % medial head
v_gas2=gas_org2-gas_ins; % lateral head
L_gas1=sqrt(v_gas1(:,1).^2+v_gas1(:,2).^2+v_gas1(:,3).^2);
L_gas2=sqrt(v_gas2(:,1).^2+v_gas2(:,2).^2+v_gas2(:,3).^2);
e_gas1=v_gas1./L_gas1; e_gas2=v_gas2./L_gas2;
r_gas1=gas_ins-OankleR; r_gas2=gas_ins-OankleR;
for i=1:n
    d_gas1(i,:)=cross(r_gas1(i,:),v_gas1(i,:))/L_gas1(i);
    d_gas2(i,:)=cross(r_gas2(i,:),v_gas2(i,:))/L_gas2(i);    
end

%%% quadriceps
% vastus intermedius
v_vi=vi_org-vi_ins;
L_vi=sqrt(v_vi(:,1).^2+v_vi(:,2).^2+v_vi(:,3).^2);
e_vi=v_vi./L_vi;
r_vi=vi_ins-OkneeR;
for i=1:n
    d_vi(i,:)=cross(r_vi(i,:),v_vi(i,:))/L_vi(i);
end

% vastus medialis 
v_vml=vm_org-vm_ins; % longus head
v_vmo=vmo_org-vm_ins; % oblique head
L_vml=sqrt(v_vml(:,1).^2+v_vml(:,2).^2+v_vml(:,3).^2);
L_vmo=sqrt(v_vmo(:,1).^2+v_vmo(:,2).^2+v_vmo(:,3).^2);
e_vml=v_vml./L_vml; e_vmo=v_vmo./L_vmo;
r_vml=vm_ins-OkneeR; r_vmo=vm_ins-OkneeR;
for i=1:n
    d_vml(i,:)=cross(r_vml(i,:),v_vml(i,:))/L_vml(i);
    d_vmo(i,:)=cross(r_vmo(i,:),v_vmo(i,:))/L_vmo(i);
end

% vastus lateralis 
v_vll=vl_org-vl_ins; % longus head
v_vlo=vlo_org-vl_ins; % oblique head
L_vll=sqrt(v_vll(:,1).^2+v_vll(:,2).^2+v_vll(:,3).^2);
L_vlo=sqrt(v_vlo(:,1).^2+v_vlo(:,2).^2+v_vlo(:,3).^2);
e_vll=v_vll./L_vll; e_vlo=v_vlo./L_vlo;
r_vll=vl_ins-OkneeR; r_vlo=vl_ins-OkneeR;
for i=1:n
    d_vll(i,:)=cross(r_vll(i,:),v_vll(i,:))/L_vll(i);
    d_vlo(i,:)=cross(r_vlo(i,:),v_vlo(i,:))/L_vlo(i);
end

% rectus femoris
v_rf=rf_org-rf_ins;
L_rf=sqrt(v_rf(:,1).^2+v_rf(:,2).^2+v_rf(:,3).^2);
e_rf=v_rf./L_rf;
r_rf=rf_ins-OkneeR;
for i=1:n
    d_rf(i,:)=cross(r_rf(i,:),v_rf(i,:))/L_rf(i);
end

%%% hamstrings
% semimembranosus
v_sm=sm_org-sm_ins;
L_sm=sqrt(v_sm(:,1).^2+v_sm(:,2).^2+v_sm(:,3).^2);
e_sm=v_sm./L_sm;
r_sm=sm_ins-OkneeR;
for i=1:n
    d_sm(i,:)=cross(r_sm(i,:),v_sm(i,:))/L_sm(i);
end

% semitendinosus
v_st=st_org-st_ins;
L_st=sqrt(v_st(:,1).^2+v_st(:,2).^2+v_st(:,3).^2);
e_st=v_st./L_st;
r_st=st_ins-OkneeR;
for i=1:n
    d_st(i,:)=cross(r_st(i,:),v_st(i,:))/L_st(i);
end

% biceps femoris
v_bf1=bf_org1-bf_ins; % long head
v_bf2=bf_org2-bf_ins; % short head
L_bf1=sqrt(v_bf1(:,1).^2+v_bf1(:,2).^2+v_bf1(:,3).^2);
L_bf2=sqrt(v_bf2(:,1).^2+v_bf2(:,2).^2+v_bf2(:,3).^2);
e_bf1=v_bf1./L_bf1; e_bf2=v_bf2./L_bf2;
r_bf1=bf_ins-OkneeR; r_bf2=bf_ins-OkneeR;
for i=1:n
    d_bf1(i,:)=cross(r_bf1(i,:),v_bf1(i,:))/L_bf1(i);
    d_bf2(i,:)=cross(r_bf2(i,:),v_bf2(i,:))/L_bf2(i);    
end

% PCSA of muscles
PCSA_tib=955e-6; % physiological cross sectional area, m^2 
PCSA_sol=3226e-6; 
PCSA_gas1=2371e-6;
PCSA_gas2=1159e-6;
PCSA_vi=2938e6;
PCSA_vml=0.5*2707e-6;
PCSA_vmo=0.5*2707e-6;
PCSA_vll=0.5*3206e-6
PCSA_vlo=0.5*3206e-6;
PCSA_rf=1853e-6;
PCSA_sm=1561e-6;
PCSA_st=1073e-6;
PCSA_bf1=998e-6;
PCSA_bf2=849e-6;

PCSA_a=[PCSA_tib PCSA_sol PCSA_gas1 PCSA_gas2]; 
PCSA_k=[PCSA_vi PCSA_vml PCSA_vmo PCSA_vll PCSA_vlo PCSA_rf PCSA_sm PCSA_st PCSA_bf1 PCSA_bf2]; 
Smax=0.6e6; %Maximum muscle stress: 0.1e6-1e6 Pa=N/m^2 

% Static Optimization 
p=3; % power

% ankle optimization
lba=zeros(1,4); uba=Smax*PCSA_a;
f0a=10*ones(1,4);
for j=1:3
    Aeqa(:,:,j)=[d_tib(:,j) d_sol(:,j) d_gas1(:,j) d_gas2(:,j)];
    beqa(:,j)=Ta(:,j);
end

% Nm=4;
% for i=2:Nm
%     OFa=@(f) (f(i)/PCSA(i))^p;
% end
OFa=@(fa) ((fa(1)/PCSA_a(1))^p+(fa(2)/PCSA_a(2))^p+(fa(3)/PCSA_a(3))^p+(fa(4)/PCSA_a(4))^p);

for i=1:n
    for j=1:3
    Fma(i,:,j)=fmincon(OFa,f0a,[],[],Aeqa(i,:,j),beqa(i,j),lba,uba);
    end
end

% ankle muscle forces
v_Ftib=[Fma(:,1,1) Fma(:,1,2) Fma(:,1,3)];
v_Fsol=[Fma(:,2,1) Fma(:,2,2) Fma(:,2,3)];
v_Fgas1=[Fma(:,3,1) Fma(:,3,2) Fma(:,3,3)];
v_Fgas2=[Fma(:,4,1) Fma(:,4,2) Fma(:,4,3)];
v_Fgas=v_Fgas1+v_Fgas2;

Ftib=sqrt(v_Ftib(:,1).^2+v_Ftib(:,2).^2+v_Ftib(:,3).^2);
Fsol=sqrt(v_Fsol(:,1).^2+v_Fsol(:,2).^2+v_Fsol(:,3).^2);
Fgas=sqrt(v_Fgas(:,1).^2+v_Fgas(:,2).^2+v_Fgas(:,3).^2);

% ankle total joint force 
Ra=m_fo*(acm_fo-g)-GRF-v_Ftib-v_Fsol-v_Fgas;
Fta=sqrt(Ra(:,1).^2+Ra(:,2).^2+Ra(:,3).^2);

% knee optimization
Tgas1=cross(r_gas1,v_Fgas1);
Tgas2=cross(r_gas2,v_Fgas2);
lbk=zeros(1,10); ubk=Smax*PCSA_k;
f0k=10*ones(1,10);
for j=1:3
    Aeqk(:,:,j)=[d_vi(:,j) d_vml(:,j) d_vmo(:,j) d_vll(:,j) d_vlo(:,j) d_rf(:,j) d_sm(:,j) d_st(:,j) d_bf1(:,j) d_bf2(:,j)];
    beqk(:,j)=Tk(:,j)-Tgas1(:,j)-Tgas2(:,j); 
end

OFk=@(fk) ((fk(1)/PCSA_k(1))^p+(fk(2)/PCSA_k(2))^p+(fk(3)/PCSA_k(3))^p+(fk(4)/PCSA_k(4))^p);

for i=1:n
    for j=1:3
    Fmk(i,:,j)=fmincon(OFk,f0k,[],[],Aeqk(i,:,j),beqk(i,j),lbk,ubk);
    end
end

% knee muscle forces
v_Fvi=[Fmk(:,1,1) Fmk(:,1,2) Fmk(:,1,3)];
v_Fvml=[Fmk(:,2,1) Fmk(:,2,2) Fmk(:,2,3)];
v_Fvmo=[Fmk(:,3,1) Fmk(:,3,2) Fmk(:,3,3)];
v_Fvm=v_Fvml+v_Fvmo;
v_Fvll=[Fmk(:,4,1) Fmk(:,4,2) Fmk(:,4,3)];
v_Fvlo=[Fmk(:,5,1) Fmk(:,5,2) Fmk(:,5,3)];
v_Fvl=v_Fvll+v_Fvlo;
v_Frf=[Fmk(:,6,1) Fmk(:,6,2) Fmk(:,6,3)];
v_Fsm=[Fmk(:,7,1) Fmk(:,7,2) Fmk(:,7,3)];
v_Fst=[Fmk(:,8,1) Fmk(:,8,2) Fmk(:,8,3)];
v_Fbf1=[Fmk(:,9,1) Fmk(:,9,2) Fmk(:,9,3)];
v_Fbf2=[Fmk(:,10,1) Fmk(:,10,2) Fmk(:,10,3)];
v_Fbf=v_Fbf1+v_Fbf2;

Fvi=sqrt(v_Fvi(:,1).^2+v_Fvi(:,2).^2+v_Fvi(:,3).^2);
Fvm=sqrt(v_Fvm(:,1).^2+v_Fvm(:,2).^2+v_Fvm(:,3).^2);
Fvl=sqrt(v_Fvl(:,1).^2+v_Fvl(:,2).^2+v_Fvl(:,3).^2);
Frf=sqrt(v_Frf(:,1).^2+v_Frf(:,2).^2+v_Frf(:,3).^2);
Fsm=sqrt(v_Fsm(:,1).^2+v_Fsm(:,2).^2+v_Fsm(:,3).^2);
Fst=sqrt(v_Fst(:,1).^2+v_Fst(:,2).^2+v_Fst(:,3).^2);
Fbf=sqrt(v_Fbf(:,1).^2+v_Fbf(:,2).^2+v_Fbf(:,3).^2);

% knee total joint force 
Rk=m_sh*(acm_sh-g)+Fa-v_Fvi-v_Fvm-v_Fvl-v_Frf-v_Fsm-v_Fst-v_Fbf;
Ftk=sqrt(Rk(:,1).^2+Rk(:,2).^2+Rk(:,3).^2);

%% Plot results 

% kinematic results
figure('WindowState','Maximized');
subplot(2,2,1)
plot(t,[a_pel b_pel g_pel]);
title('pelvis angles, deg');
legend('\alpha','\beta','\gamma');
subplot(2,2,2)
plot(t,[a_th b_th g_th]);
title('thigh angles, deg');
legend('\alpha','\beta','\gamma');
subplot(2,2,3)
plot(t,[a_sh b_sh g_sh]);
title('shank angles, deg');
legend('\alpha','\beta','\gamma');
subplot(2,2,4)
plot(t,[a_fo b_fo g_fo]);
title('foot angles, deg');
legend('\alpha','\beta','\gamma');

figure('WindowState','Maximized');
subplot(3,3,1)
plot(t,[a_h b_h g_h]);
title('hip joint angles, deg');
legend('\alpha','\beta','\gamma');
subplot(3,3,2)
plot(t,omga_h);
title('hip joint angular velocity, rad/s');
legend('\omega_x','\omega_y','\omega_z');
subplot(3,3,3)
plot(t,alfa_h);
title('hip joint angular acceleration, rad/s^2');
legend('\alpha_x','\alpha_y','\alpha_z');

subplot(3,3,4)
plot(t,[a_k b_k g_k]);
title('knee joint angles, deg');
legend('\alpha','\beta','\gamma');
subplot(3,3,5)
plot(t,omga_k);
title('knee joint angular velocity, rad/s');
legend('\omega_x','\omega_y','\omega_z');
subplot(3,3,6)
plot(t,alfa_k);
title('knee joint angular acceleration, rad/s^2');
legend('\alpha_x','\alpha_y','\alpha_z');

subplot(3,3,7)
plot(t,[a_a b_a g_a]);
title('ankle joint angles, deg');
legend('\alpha','\beta','\gamma');
subplot(3,3,8)
plot(t,omga_a);
title('ankle joint angular velocity, rad/s');
legend('\omega_x','\omega_y','\omega_z');
subplot(3,3,9)
plot(t,alfa_a);
title('ankle joint angular acceleration, rad/s^2');
legend('\alpha_x','\alpha_y','\alpha_z');

% kinetic results
figure('WindowState','Maximized')
subplot(3,3,1), plot(t,Fh);
title('hip net joint force, N');
legend('Fx','Fy','Fz');
subplot(3,3,2), plot(t,Th);
title('hip net joint moment, N.m');
legend('Mx','My','Mz');
subplot(3,3,3), plot(t,Ph);
title('hip joint power, Watts');
legend('Px','Py','Pz');
subplot(3,3,4), plot(t,Fk);
title('knee net joint force, N');
legend('Fx','Fy','Fz');
subplot(3,3,5), plot(t,Tk);
title('knee net joint moment, N.m');
legend('Mx','My','Mz');
subplot(3,3,6), plot(t,Pk);
title('knee joint power, Watts');
legend('Px','Py','Pz');
subplot(3,3,7), plot(t,Fa);
title('ankle net joint force, N');
legend('Fx','Fy','Fz');
subplot(3,3,8), plot(t,Ta);
title('ankle net joint moment, N.m');
legend('Mx','My','Mz');
subplot(3,3,9), plot(t,Pa);
title('ankle joint power, Watts');
legend('Px','Py','Pz');

% muscle force results
% ankle muscles
figure('WindowState','Maximized');
subplot(3,1,1)
plot(t,Ftib); title('Tibialis anterior');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(3,1,2)
plot(t,Fsol); title('Soleus');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(3,1,3)
plot(t,Fgas); title('Gastrocnemius');
xlabel('Muscle force, N'); ylabel('Time, sec');
% quadriceps 
figure('WindowState','Maximized');
subplot(2,2,1)
plot(t,Fvi); title('Vastus intermedius');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(2,2,2)
plot(t,Fvm); title('Vastus medialis');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(2,2,3)
plot(t,Fvl); title('Vastus lateralis');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(2,2,4)
plot(t,Frf); title('Rectus femoris');
xlabel('Muscle force, N'); ylabel('Time, sec');
% hamstrings 
figure('WindowState','Maximized');
subplot(3,1,1)
plot(t,Fsm); title('Semimembranosus');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(3,1,2)
plot(t,Fst); title('Semitendinosus');
xlabel('Muscle force, N'); ylabel('Time, sec');
subplot(3,1,3)
plot(t,Fbf); title('Biceps Femoris');
xlabel('Muscle force, N'); ylabel('Time, sec');

% total joint forces 
figure('WindowState','Maximized');
subplot(2,2,1)
plot(t,Ra); title('ankle joint force');
xlabel('Joint force, N'); ylabel('Time, sec');
legend('Fx','Fy','Fz');
subplot(2,2,2)
plot(t,Fta); title('ankle total joint force');
xlabel('Joint force, N'); ylabel('Time, sec');
subplot(2,2,3)
plot(t,Rk); title('knee total joint force');
xlabel('Joint force, N'); ylabel('Time, sec');
legend('Fx','Fy','Fz');
subplot(2,2,4)
plot(t,Ftk); title('knee total joint force');
xlabel('Joint force, N'); ylabel('Time, sec');