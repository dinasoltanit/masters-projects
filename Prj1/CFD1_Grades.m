clc; clear; close all;
% ---------------- Input Parameters --------------------------------------
% Dear user,
% you can read your data from an excel file.
% All the best

%% Method: provide the inputs from an excel file
% Your excel file should be something like this
% Row|Name|HW1 |HW2 |HW3 |HW4 |HW5 |Prj1 |Prj2 |Prj3 |Prj4 |Mid |Final|Delays|Absents
% 1  |XXXX|18.0|19.0|19.5|20  |20  |18   |19   |19.5 |20   |20  |20   |5     |2
% 2  |YYYY|...
% read the input data
[numbers, strings, raw] = xlsread('CFD1_grades.xlsx');
HW1 = 0.02*(numbers(1:10,3));
HW2 = 0.02*(numbers(1:10,4));
HW3 = 0.02*(numbers(1:10,5));
HW4 = 0.02*(numbers(1:10,6));
HW5 = 0.02*(numbers(1:10,7));
Prj1 = 0.06*(numbers(1:10,8));
Prj2 = 0.08*(numbers(1:10,9));
Prj3 = 0.08*(numbers(1:10,10));
Prj4 = 0.08*(numbers(1:10,11));
Mid = 0.2*(numbers(1:10,12));
Final = 0.4*(numbers(1:10,13));
Delays = floor((numbers(1:10,14))./3)
Absents = (numbers(1:10,15));
% compute the grade
HW_Grades = HW1 + HW2 + HW3 + HW4 + HW5;
Prj_Grades = Prj1 + Prj2 + Prj3 + Prj4;
Grades = HW_Grades + Prj_Grades + Mid + Final + (4- Delays - Absents)*0.25;
% write the output data
filename = 'CFD1_FinalGrades.xlsx';
names = strings(2:11,2);
A = string(names);
sheet1 = 1;
xlRange1 = 'A1';
xlswrite(filename,A,sheet1,xlRange1)

B = Grades;
sheet1 = 1;
xlRange2 = 'B1';
xlswrite(filename,B,sheet1,xlRange2)





