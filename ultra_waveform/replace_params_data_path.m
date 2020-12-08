% Read txt into cell A
fid = fopen('params.py','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

A1 = A{1};
A1 = A1(13:end-1);
A1 = ['"' rootZ '\' A1 '"'];
A2 = ['dat_path = r' A1];

% Change cell A
A{1} = A2;
% Write cell A into txt
fid = fopen('params.py', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end