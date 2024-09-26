function [vertices, faces] = readObj(filename)
% 打开文件
fileID = fopen(filename,'r');
if fileID == -1
    error('无法打开文件');
end

vertices = [];
faces = [];

% 逐行读取文件
line = fgetl(fileID);
while ischar(line)
    % 解析顶点坐标行
    if startsWith(line,'v ')
        parts = strsplit(line);
        vertex = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];
        vertices = [vertices; vertex];
    % 解析面索引行
    elseif startsWith(line,'f ')
        parts = strsplit(line);
        face = zeros(1, length(parts)-1);
        for i = 2:length(parts)
            faceIndexParts = strsplit(parts{i},'/');
            face(i-1) = str2double(faceIndexParts{1});
        end
        faces = [faces; face];
    end
    line = fgetl(fileID);
end

fclose(fileID);
end