function infer_ear_eyp

% Basal layer data on ear epidermis. Should very sensitive to the scaling
% limit, but be significantly worse for distinguishing transient behaviour.
% Works almost orthogonally to the short term basal+suprabasal experiment.

t1 = 3;
data1 = [2 65; 3 19; 4 8; 5 1; 6 1];

t2 = 6;
data2 = [2 112; 3 53; 4 27; 5 17; 6 6; 7 1; 8 3];

t3 = 26;
data3 = [2 57; 3 38; 4 20; 5 24; 6 15; 7 8; 8 2; 9 8; 10 6; 11 2; 12 2; 13 5; 16 1; 20 1; 22 1];

t4 = 52;
data4 = [2 33; 3 26; 4 21; 5 15; 6 10; 7 9; 8 6; 9 8; 10 6; 11 1; 12 5; 13 6; 14 2; 15 3; 16 7; 17 2; 19 1; 20 1; 21 1; 25 2; 35 1; 44 1; 50 1];

infer_eyp(linspace(0.1, 0.8, 31), linspace(0.0, 0.5, 30), linspace(1/6, 1/3, 29), [t1 t2 t3 t4], {data1 data2 data3 data4}, 'ear-eyp-b');

end
