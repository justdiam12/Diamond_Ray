function categ = category (Klass,n) 
allcateg=[...
    0,   0,   0,   0,    0;  % Klass = 1
    
    0,   0,   n-1, n,    1;  % Klass = 2
    0,   0,   n,   n,   -1;  % Klass = 3
    0,   0,   n,   n,    1;  % Klass = 4
    0,   0,   n,   n-1, -1;  % Klass = 5

    n,   n+1, 0,   0,   -1;  % Klass = 6
    n+1, n,   0,   0,    1;  % Klass = 7
    n,   n,   0,   0,   -1;  % Klass = 8
    n,   n,   0,   0,    1;  % Klass = 9

    n,   0,   n-1, 0,    1;  % Klass = 10
    n,   0,   n,   0,   -1;  % Klass = 11
    n,   0,   n,   0,    1;  % Klass = 12
    n,   0,   n+1, 0,   -1;  % Klass = 13

    0,   n,   0,   n-1, -1;  % Klass = 14
    0,   n,   0,   n,   -1;  % Klass = 15
    0,   n,   0,   n,    1;  % Klass = 16
    0,   n,   0,   n+1,  1]; % Klass = 17
%  bot  suf  abv  blw   dir
categ=allcateg(Klass,:);