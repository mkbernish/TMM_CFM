addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\"))
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
basepath = 'C:\Users\mathp\OneDrive\Documents\midterm_stud\';
t1 = readtable("C:\Users\mathp\OneDrive\Documents\midterm_stud\01 All_studentd_Midterm_K.xlsx");
for i = 1:height(t1)
    ln = string(t1(i,:).LastName);
    fn = append('Midterm_',ln,'.xlsx');
    ffn = append(basepath,fn);
    writetable(t1(i,:),ffn);
end