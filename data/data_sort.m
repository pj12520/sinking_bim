function data_sort()
%data_sort Sorts data from the parameter sweep into smaller files
%   Detailed explanation goes here

%Read in data

infile = 'regime.dat';

fid = fopen(infile);
data = textscan(fid, '%f %f %f %s', 'HeaderLines', 1);
fclose(fid);

%Get size of data
data_size = size(data{1},1);

%Extract data.
for i = 1:data_size
    mdr_data(i) = data{1}(i);
    bond_data(i) = data{2}(i);
    viscos_rat_data(i) = data{3}(i);
    mode(i) = data{4}(i);
end

%Values used in sweep
mdr(1) = 1;
mdr(2) = 10;
mdr(3) = 100;
mdr(4) = 1000;

bond(1) = 0.001;
bond(2) = 0.01;
bond(3) = 0.1;
bond(4) = 1;
bond(5) = 10;
bond(6) = 100;
bond(7) = 1000;

viscos(1) = 0.001;
viscos(2) = 0.01;
viscos(3) = 0.1;
viscos(4) = 1;
viscos(5) = 10;
viscos(6) = 100;
viscos(7) = 1000;

%Output data
for i = 1:4
    mdr_file = strcat('const_mdr/mdr',num2str(mdr(i)),'.dat');
    mdr_out(i) = fopen(mdr_file,'w');
    fprintf(mdr_out(i), '%s\t %s\t %s\t %s\n', '#mdr', 'bond', 'viscos_rat', 'regime');
    
    for j = 1:data_size
        if mdr_data(j) == mdr(i)
            fprintf(mdr_out(i), '%f\t %f\t %f\t %s\n', mdr_data(j), bond_data(j), viscos_rat_data(j), mode{j});
        end
    end
    
    fclose(mdr_out(i));
    
end

for i = 1:7
   bond_file = strcat('const_bond/bond',num2str(bond(i)),'.dat');
   bond_out(i) = fopen(bond_file, 'w'); 
   fprintf(bond_out(i), '%s\t %s\t %s\t %s\n', '#mdr', 'bond', 'viscos_rat', 'regime');

   viscos_file = strcat('const_viscos/viscos',num2str(viscos(i)),'.dat');
   viscos_out(i) = fopen(viscos_file, 'w'); 
   fprintf(viscos_out(i), '%s\t %s\t %s\t %s\n', '#mdr', 'bond', 'viscos_rat', 'regime');
   
   for j = 1:data_size
       if bond_data(j) == bond(i)
           fprintf(bond_out(i), '%f\t %f\t %f\t %s\n', mdr_data(j), bond_data(j), viscos_rat_data(j), mode{j});
       end
       
       if viscos_rat_data(j) == viscos(i)
           fprintf(viscos_out(i), '%f\t %f\t %f\t %s\n', mdr_data(j), bond_data(j), viscos_rat_data(j), mode{j});
       end
   end

    fclose(bond_out(i));
    fclose(viscos_out(i));
   
end

end

