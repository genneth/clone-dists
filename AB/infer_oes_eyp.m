function infer_oes_eyp

t1 = 3/7;
data1 = [2 71; 3 8; 4 3];

t2 = 10/7;
data2 = [2 109; 3 35; 4 15; 5 4; 6 2];

t3 = 3;
data3 = [2 79; 3 38; 4 25; 5 11; 6 4; 7 1; 8 1; 11 1];

t4 = 6;
data4 = [2 60; 3 44; 4 26; 5 22; 6 10; 7 3; 8 4; 10 2; 11 4; 13 2; 14 1; 20 1; 28 1];

t5 = 12;
data5 = [2 87; 3 50; 4 52; 5 26; 6 24; 7 14; 8 8; 9 9; 10 8; 11 6; 12 7; 13 4; 14 3; 15 3; 16 2; 17 2; 20 2; 22 1; 26 1; 29 1];

t6 = 26;
data6a = [2 59; 3 43; 4 18; 5 27];
data6b = [6 23; 7 16; 8 12; 9 14; 10 10; 11 11; 12 5; 13 8; 14 7; 15 6; 16 3; 17 4; 18 3; 19 3; 20 4; 21 1; 22 2; 23 1; 24 3; 25 4; 27 1; 30 1; 33 1; 34 1; 35 5; 36 1; 37 2; 39 1; 50 2; 56 1; 85 1; 90 1];

t7 = 52;
data7a = [2 18; 3 15; 4 14; 5 6; 6 11; 7 8; 8 6; 9 9; 10 9; 11 3; 12 4];
data7b = [13 4; 14 2; 15 4; 16 4; 17 2; 18 4; 19 2; 20 2; 21 3; 22 4; 23 1; 24 5; 26 2; 27 4; 28 3; 19 1; 30 3; 31 3; 32 3; 33 2; 34 4; 35 2; 36 1; 37 1; 38 1; 40 1; 41 1; 42 1; 43 1; 44 1; 45 1; 47 2; 48 1; 49 1; 54 1; 55 3; 58 1; 60 2; 62 1; 63 1; 70 2; 75 1; 80 1; 85 1; 96 1; 110 1; 130 1; 140 1; 175 1];

rhos = linspace(0.1, 0.8, 31);
rs = linspace(0.0, 0.5, 30);
%lambdas = linspace(0.5, 1.5, 10);

%infer_eyp(rhos, rs, 1.0624, [t1], {data1}, 'oes-eyp-b-fixed-lambda-day-3');
%infer_eyp(rhos, rs, 1.0624, [t2], {data2}, 'oes-eyp-b-fixed-lambda-day-10');
%infer_eyp(rhos, rs, 1.0624, [t3], {data3}, 'oes-eyp-b-fixed-lambda-week-3');
%infer_eyp(rhos, rs, 1.0624, [t4], {data4}, 'oes-eyp-b-fixed-lambda-week-6');
%infer_eyp(rhos, rs, 1.0624, [t5], {data5}, 'oes-eyp-b-fixed-lambda-week-12');
infer_eyp(rhos, rs, 1.0624, [t6 t6], {data6a data6b}, 'oes-eyp-b-fixed-lambda-week-26');
%infer_eyp(rhos, rs, 1.0624, [t5 t1 t4 t2 t3], {data5 data1 data4 data2 data3}, 'oes-eyp-b-fixed-lambda-upto-12-weeks');
%infer_eyp(rhos, rs, lambdas, [t7 t7 t1 t6 t6 t2 t5 t3 t4], {data7a data7b data1 data6b data6b data2 data5 data3 data4}, 'oes-eyp-b');

end
