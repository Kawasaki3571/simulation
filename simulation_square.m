i = 0 : 1000;
f = i + 2;
q = [i; f];
csvwrite('myFile.csv',q)