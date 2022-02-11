% Dear user,
% you may read your data from an excel file. If you want to enter you
% desired grading weights, enter the <<mode>> to be <<manual>> and enter
% <<0>> as other inputs. Otherwise, the grading weights would be as the
% defual values, provided by the question.
% All the best
% Method of providing the inputs from an excel file
% Your excel file should be something like this
% Row|Name|HW1 |HW2 |HW3 |HW4 |HW5 |Prj1 |Prj2 |Prj3 |Prj4 |Mid |Final|Delays|Absents
% 1  |XXXX|18.0|19.0|19.5|20  |20  |18   |19   |19.5 |20   |20  |20   |5     |2
% 2  |YYYY|...
% Go ahead and good luck ;)
function [Grades, names] = CFD1_Grader(mode, NumOfStudents, w_hw1, w_hw2, w_hw3, w_hw4, w_hw5, w_prj1, w_prj2, w_prj3, w_prj4, w_mid, w_fin, inputFile, outputFile)

manual = 0;
defaul = 1;

if strcmpi (mode, 'manuall')
    chosen_mode = manual;
elseif strcmpi (mode, 'default')
    chosen_mode = defaul;
else
    disp ('ERROR: error in parsing the argument "mode".');
    return;
end

switch chosen_mode
    case manual
        %% the manual case, based on the user input data
        
        [numbers, strings, raw] = xlsread(inputFile);
        HW1 = w_hw1*(numbers(1:NumOfStudents,3));
        HW2 = w_hw2*(numbers(1:NumOfStudents,4));
        HW3 = w_hw3*(numbers(1:NumOfStudents,5));
        HW4 = w_hw4*(numbers(1:NumOfStudents,6));
        HW5 = w_hw5*(numbers(1:NumOfStudents,7));
        Prj1 = w_prj1*(numbers(1:NumOfStudents,8));
        Prj2 = w_prj2*(numbers(1:NumOfStudents,9));
        Prj3 = w_prj3*(numbers(1:NumOfStudents,10));
        Prj4 = w_prj4*(numbers(1:NumOfStudents,11));
        Mid = w_mid*(numbers(1:NumOfStudents,12));
        Final = w_fin*(numbers(1:NumOfStudents,13));
        Delays = floor((numbers(1:NumOfStudents,14))./3);
        Absents = (numbers(1:NumOfStudents,15));
        % compute the grade
        HW_Grades = HW1 + HW2 + HW3 + HW4 + HW5;
        Prj_Grades = Prj1 + Prj2 + Prj3 + Prj4;
        Grades = HW_Grades + Prj_Grades + Mid + Final + (4- Delays - Absents)*0.25;
        % write the output data
        filename = outputFile;
        names = strings((2:NumOfStudents+1),2);
        A = string(names);
        sheet1 = 1;
        xlRange1 = 'A1';
        xlswrite(filename,A,sheet1,xlRange1)

        B = Grades;
        sheet1 = 1;
        xlRange2 = 'B1';
        xlswrite(filename,B,sheet1,xlRange2)

    case defaul
        %% the default case, based on the project description
        %NumOfStudents = 10;
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
end

end





