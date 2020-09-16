clc;
close all;
clear all;

% Data is taken from as directed by the following manuscript
% F. Zhu, Y. Wang, S. Xiang, B. Fan, and C. Pan, “Structured sparse method for hyperspectral unmixing,” ISPRS J. Photogrammetry Remote Sens., vol. 88, pp. 101–118, 2014.
%% Image Read
s=load('Cuprite.mat');
p=s.nRow;
q=s.nCol;
Bands=188;
Y=s.Y;
x=hyperConvert3d(Y,p,q,Bands);

%% Hysime
n_Hysime=12;

%VCA
[U, endmemberindex_2] = hyperVca(Y,n_Hysime);
endmemberindex_VCA=change_index(endmemberindex_2,p,q);

%Energy
[endmemberindex_1] = Energy(Y,n_Hysime);
endmemberindex_Energy=change_index(endmemberindex_1,p,q);

%gt compare
t1=load('groundTruth_Cuprite_nEnd12.mat');
gt=t1.M;
tit=t1.cood;
n1=gt(3:103,:);
n2=gt(114:147,:);
n3=gt(168:220,:);
gt=[n1;n2;n3];

[gt_m,gt_n]=size(gt);

for i=1:gt_n
    for j=1:Bands
        extracted_Energy(j,i)=x(endmemberindex_Energy(i,1),endmemberindex_Energy(i,2),j);
        extracted_VCA(j,i)=x(endmemberindex_VCA(i,1),endmemberindex_VCA(i,2),j);
    end
end

%% SAM Calculation
ex_n=gt_n;

for i=1:gt_n
    for j=1:ex_n
        %gt vs extraceted
        Mat_SAM_Energy(i,j)=real(acos(dot(gt(:,i),extracted_Energy(:,j))/(norm(gt(:,i)*norm(extracted_Energy(:,j))))));
        Mat_SAM_VCA(i,j)=real(acos(dot(gt(:,i),extracted_VCA(:,j))/(norm(gt(:,i)*norm(extracted_VCA(:,j))))));
    end
end

store_Energy=[0,0];
store_VCA=[0,0];

sam_Energy=0;
sam_VCA=0;

sam_total_Energy=0;
sam_total_VCA=0;

for i=1:gt_n
    %Energy
    [max_value1,mrow]=min(Mat_SAM_Energy);
    [max_value,col_Energy]=min(max_value1);
    sam_total_Energy=sam_total_Energy+max_value;
    sam_Energy=[sam_Energy;max_value];
    row_Energy=mrow(col_Energy);
    s1=[row_Energy,col_Energy];
    store_Energy=[store_Energy;s1];
    save_Energy(row_Energy)=max_value;
    Mat_SAM_Energy(row_Energy,:)=[100*ones];
    Mat_SAM_Energy(:,col_Energy)=[100*ones];
    e_Energy(:,row_Energy)=extracted_Energy(:,col_Energy);
    %VCA
    [max_value1,mrow]=min(Mat_SAM_VCA);
    [max_value,col_VCA]=min(max_value1);
    sam_total_VCA=sam_total_VCA+max_value;
    sam_VCA=[sam_VCA;max_value];
    row_VCA=mrow(col_VCA);
    s1=[row_VCA,col_VCA];
    store_VCA=[store_VCA;s1];
    save_VCA(row_VCA)=max_value;
    Mat_SAM_VCA(row_VCA,:)=[100*ones];
    Mat_SAM_VCA(:,col_VCA)=[100*ones];
    e_VCA(:,row_VCA)=extracted_VCA(:,col_VCA);
end

s_E=store_Energy(2:end,:);
s_CM=store_VCA(2:end,:);

for i=1:gt_n
    s_E(s_E(i,1),3)=s_E(i,2);
    s_CM(s_CM(i,1),3)=s_CM(i,2);
end

rms_sae=[rms(save_Energy);
    rms(save_VCA)];

rms_sae = radtodeg(rms_sae);
disp('RMSSAE of VCA');
disp(radtodeg(rms(save_VCA)));

disp('RMSSAE of Energy');
disp(radtodeg(rms(save_Energy)));