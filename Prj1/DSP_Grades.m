close all;clear all;clc;format bank;
% TA    PROJECT MId1    Mid2    Final   
a=ones(41,5); % do not change this....
% ------- @@@@@@@@@@@ Modify these two lines @@@@@@@@@@@@@@@ -------------
List=3; % Put your list number here , NEW LIST NUMBER 1 to 41 ...
% PUT your 5 marks here... 
a(List,:)=[...
% TA (from 4)   Project (from 100), Mid1 (20),  Mid2 (20),  Final (from 15). ...
3.98           79.17               10.7       16.15       11 
]; 
%-------------@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@--------------
Ta=a(:,1);Proj=a(:,2);
M1=a(:,3);M2=a(:,4);F=a(:,5)*20/15;
%-------------------------------------------------------
% Mid1 Bonuses 
Mid1_Bonus=[
        .75 % Tayebi
		.25 % Seyed Hosseini
		.50 % Tabatabian
		.50 % Karimi
		.25 % Kourehpaz
		.25 % Mohadeseh Sadeghi
		.25 % Amirkhanian
		.25 % Salehian
    	.25];% Jokaar
M1_B=zeros(1,41);
M1_B([7,17,21,29,30,37,38,39,41])=Mid1_Bonus(:); % Mid-Term1 Bonuses
M1=M1+M1_B(:); 
%  Those absent with medical excuse in Mid term 1
M1_Absent=[14,22,31,36]; 
% Mid2 Bonuses: Tabtabaian, Arabzadeh, Kourehpaz
M2([21,23,30])= M2([21,23,30]) + [.25,.5,.25]';
F(29)=F(29)+.25; % Final Bonus: Karimi
for k=1:41
    temp=[M1(k),M2(k),F(k)];
    if find(k==M1_Absent) % Those absent with excuse in Mid term 1 
       temp_Abs=[M2(k),F(k)];
    else
        temp_Abs=temp;
    end
    [Max,ix]=max(temp);
    temp(ix)=[];
    Mark1=0.4*Max+ .6*(mean(temp)); %  40% for the max. exam mark
    % ----------------------
    [Max_Abs,ix]=max(temp_Abs);
    temp_Abs(ix)=[];
    Mark_Abs=0.4*Max_Abs+ .6*(mean(temp_Abs));
    %-------------------------
    Mark(k)=max(Mark1,Mark_Abs);
    %-------------------------
    Mark(k)= Mark(k)*16/20 + Ta(k);
    if Proj(k)>= 75 % For those with high marks in Matlab Projects and apps
       temp= Mark(k)*15/20 + Ta(k)*5/4;
       if temp> Mark(k)
           disp([num2str(k),': 5 grades for Ta ',' Extra Mark= ',...
               num2str(temp-Mark(k))]);
           Mark(k)=temp;
       end
    end
    Mark(k)=round(Mark(k)*10)/10;
    if Mark(k)>20, Mark(k)=20;end
end
disp(['My Final mark= ',num2str(Mark(List))]);